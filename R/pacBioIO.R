## Copyright (c) 2010, Pacific Biosciences of California, Inc.

## All rights reserved.

## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:

##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.

##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.

##     * Neither the name of Pacific Biosciences nor the names of its
##       contributors may be used to endorse or promote products derived
##       from this software without specific prior written permission.

## THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED

## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS CONTRIBUTORS
## BE LIABLE FOR ANY

## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
## GOODS OR SERVICES;

## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
## LIABILITY, OR TORT

## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## ###################################################################
##
## R code for interfacing to *.h5 files 
##
## This file manages data access to the files and their object
## model. More complicated functionality, i.e., manipulation of
## alignments, joining data together is deal with in the *Util.R
## files.
##
## R is 1-based. h5 is zero-based. The h5r package expects 1-based
## indexing, therefore in this file 0-based values in *.h5 files need
## to be converted to 1-based values so that they are handled
## correctly by the h5r package. 
##
## ###################################################################

##
## Class definitions
##
setClassUnion("stringOrNull", c("character", "NULL"))
setClassUnion("H5ObjOrNull", c("H5Obj", "NULL"))
setClassUnion("matrixOrNull", c("matrix", "NULL"))

setClass("PacBioDataFile", contains = "H5File",
         representation = representation(
           version  = "stringOrNull"),
         prototype = prototype(version = NULL))

setClass("PacBioCmpH5", contains = "PacBioDataFile",
         representation = representation(
           AlnIndex = "data.frame",
           AlnGroup = "data.frame",
           RefGroup = "data.frame",
           MovieInfo = "data.frame",
           RefInfo = "data.frame",
           isSorted = "logical"))

setClass("PacBioBasH5", contains = "PacBioDataFile",
         representation = representation(
           baseEvents = "matrixOrNull",
           baseCallsG = "H5ObjOrNull"),
         prototype = prototype(baseEvents = NULL, baseCallsG = NULL))

setClass("PacBioPlsH5", contains = "PacBioBasH5",
         representation = representation(
           pulseEvents = "matrixOrNull",
           pulseCallsG = "H5ObjOrNull"),
         prototype = prototype(pulseEvents = NULL, pulseCallsG = NULL))
         
setClass("PacBioTrcH5", contains = "PacBioDataFile",
         representation = representation())

setClass("PacBioAlnH5", contains = "PacBioDataFile",
         representation = representation(
           ZMWs = "matrix",
           alignment = "H5Obj"))

##
## Generics - These generics represent functions which make sense for
## various types of h5 files, e.g., aln.h5 files have alignments as do
## cmp.h5 files.
##
setGeneric("getAlignments", function(h5Obj, ...) {
  standardGeneric("getAlignments")
})

setGeneric("getQualityValue", function(h5Obj, ...) {
  standardGeneric("getQualityValue")
})

setGeneric("getWidthInFrames", function(h5Obj, ...) {
  standardGeneric("getWidthInFrames")
})

setGeneric("getMovieName", function(h5Obj, ...) {
  standardGeneric("getMovieName")
})

setGeneric("getSNR", function(h5Obj, ...) {
  standardGeneric("getSNR")
})

setMethod("initialize", "PacBioDataFile", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)

  if (h5AttributeExists(.Object, "Version")) {
    version <- getH5Attribute(.Object, "Version")[]
  } else {
    version <- ""
  }
  .Object@version <- version
  return(.Object)
})

setMethod("show", "PacBioDataFile", function(object) {
  callNextMethod()
  cat("Version:", object@version, "\n")
})

getVersion <- function(pacBioDataFile) {
  return(pacBioDataFile@version)
}

## #############################################################################
##
## Compare H5 Files
##
## #############################################################################
PacBioCmpH5 <- function(fileName) {
  obj <- new("PacBioCmpH5", fileName = fileName)

  .initMap <- function(grpName, names = c("ID", "Path")) {
    group <- getH5Group(obj, grpName)
    names(names) <- names
    
    dta <- data.frame(lapply(names, function(nm) {
      if(h5DatasetExists(group, nm)) {
        getH5Dataset(group, nm)[]
      } else {
        NA
      }
    }), stringsAsFactors = FALSE)
    rownames(dta) <- dta$ID
    return(dta)
  }
  
  AlnGroup <- .initMap("AlnGroup")
  RefGroup <- .initMap("RefGroup", c("ID", "Path", "RefInfoID"))
  RefInfo <- .initMap("RefInfo", c("ID", "Name", "FullName", "Length", "MD5"))
  MovieInfo <- .initMap("MovieInfo", c("ID", "Name", "Exp", "Run"))

  ## add some things to these data structures.
  AlnGroup$RefGroupName <- sapply(strsplit(AlnGroup$Path, "/"), "[[", 2)
  
  if (h5DatasetExists(obj, "RefGroup/OffsetTable")) {
    refBlocks <- getH5Dataset(obj, "RefGroup/OffsetTable")[]
    if (is.null(dim(refBlocks))) {
      dim(refBlocks) <- c(1, 3)
    }
    refBlocks <- refBlocks[,-1, drop = FALSE]
    colnames(refBlocks) <- c("offsetBegin", "offsetEnd")
    refBlocks[,1] <- refBlocks[,1] + 1
    RefGroup <- cbind(RefGroup, refBlocks)
    isSorted <- TRUE
  } else {
    isSorted <- FALSE
  }
  
  RefGroup <- merge(RefGroup, RefInfo, by.x = "RefInfoID", by.y = "ID")
  RefGroup$Name <- gsub("^/", "", RefGroup$Path)

  obj@isSorted <- isSorted
  obj@AlnGroup <- AlnGroup
  obj@RefGroup <- RefGroup
  obj@RefInfo <- RefInfo
  obj@MovieInfo <- MovieInfo
  
  aI <- getH5Dataset(obj, "AlnInfo/AlnIndex")
  aI <- aI[]
  
  obj@AlnIndex <-
    data.frame("ID"            = aI[,1],
               "alnGroupPath"  = AlnGroup[match(aI[,2], AlnGroup$ID),  "Path"],
               "movieName"     = MovieInfo[match(aI[,3], MovieInfo$ID), "Name"],
               "refName"       = RefGroup[match(aI[,4], RefGroup$ID),  "Name"],
               "fullRefName"   = RefGroup[match(aI[,4], RefGroup$ID),  "FullName"],
               "tStart"        = aI[,5],
               "tEnd"          = aI[,6],
               "alignedStrand" = aI[,7],
               "holeNumber"    = aI[,8],
               "setNumber"     = aI[,9],
               "strobeNumber"  = aI[,10],
               "moleculeID"    = aI[,11],
               "rStart"        = aI[,12],
               "rEnd"          = aI[,13],
               "mapQV"         = aI[,14],
               "nMatches"      = aI[,15],
               "nMisMatches"   = aI[,16],
               "nInsertions"   = aI[,17],
               "nDeletions"    = aI[,18],
               "offsetBegin"   = aI[,19],
               "offsetEnd"     = aI[,20],
               "nBackRead"     = aI[,21],
               "nOverlap"      = aI[,22], stringsAsFactors = FALSE)
  
  ## adjusting for 0-v-1 based.
  obj@AlnIndex[, "tStart"] <- as.integer(obj@AlnIndex[, "tStart"] + 1)
  obj@AlnIndex[, "rStart"] <- as.integer(obj@AlnIndex[, "rStart"] + 1)
  obj@AlnIndex[, "offsetBegin"] <- as.integer(obj@AlnIndex[, "offsetBegin"] + 1)
  
  return(obj)
}

alnIndex  <- function(cmpH5) cmpH5@AlnIndex
alnGroup  <- function(cmpH5) cmpH5@AlnGroup
refGroup  <- function(cmpH5) cmpH5@RefGroup
refInfo   <- function(cmpH5) cmpH5@RefInfo
movieInfo <- function(cmpH5) cmpH5@MovieInfo
isSorted  <- function(cmpH5) cmpH5@isSorted

##
## This function gets the ref elt by Name, FullName, or ID.
##
.getRefElt <- function(cmpH5, refSeq, which) {
  index <- if (is.character(refSeq)) {
    tmp <- which(refSeq == refGroup(cmpH5)$Name)
    if (length(tmp) == 0) {
      which(refSeq == refGroup(cmpH5)$FullName)
    } else {
      tmp
    }
  } else if (is.numeric(refSeq)) {
    which(refSeq == refGroup(cmpH5)$ID)
  } else {
    stop("refSeq needs to be either integer (ID) or character (Name or FullName).")
  }
  if (length(index) != 1) {
    stop(paste("refSeq:", refSeq, "not found or not unique:", index))
  }
  refGroup(cmpH5)[index, which]
}
getRefPath <- function(cmpH5, refSeq) {
  .getRefElt(cmpH5, refSeq, "Path")
}
getRefName <- function(cmpH5, refSeq) {
  .getRefElt(cmpH5, refSeq, "Name")
}
getRefFullName <- function(cmpH5, refSeq) {
  .getRefElt(cmpH5, refSeq, "FullName")
}
getRefLength <- function(cmpH5, refSeq) {
  .getRefElt(cmpH5, refSeq, "Length")
}
.getRefBlock <- function(cmpH5, refSeq) {
  if (! isSorted(cmpH5)) {
    stop("cmpH5 file not sorted.")
  }
  as.integer(.getRefElt(cmpH5, refSeq, c("offsetBegin", "offsetEnd")))
}

setMethod("show", "PacBioCmpH5", function(object) {
  callNextMethod(object)
  cat("N Alignments:", nrow(alnIndex(object)), "\n")
  cat("N ReadGroups:", nrow(alnGroup(object)), "\n")
  cat("N RefSeqs:",    nrow(refGroup(object)), "\n")
  cat("N Read Bases:", sum(sum(getReadLength(object))), "\n")
})

setMethod("summary", "PacBioCmpH5", function(object) {
  show(object)
  print(refGroup(object))
  print(movieInfo(object))
})

setMethod("nrow", "PacBioCmpH5", function(x) {
  nrow(alnIndex(x))
})

setMethod("head", "PacBioCmpH5", function(x, ...) {
  head(alnIndex(x), ...)
})

setMethod("colnames", "PacBioCmpH5", function(x) {
  colnames(alnIndex(x))
})

setMethod("$", "PacBioCmpH5", function(x, name) {
  alnIndex(x)[,name]
})

setMethod("[", c("PacBioCmpH5", "ANY", "ANY", "ANY"), function(x, i, j) {
  alnIndex(x)[i, j, drop = FALSE]
})

setMethod("getMovieName", "PacBioCmpH5", function(h5Obj, idx = seq.int(1, nrow(h5Obj))) {
  h5Obj$movieName[idx]
})


##
## This stores the mapping from pairs to characters.
##
.NUMBERS <- c(17,18,20,24,16,31,33,34,36,40,32,47,65,66,68,72,64,79,
              129,130,132,136,128,143,1,2,4,8,0,15,241,242,244,248,240,255)
.PAIRS   <- c("AA", "AC", "AG", "AT", "A-", "AN", "CA", "CC", "CG", "CT", "C-",
              "CN", "GA", "GC", "GG", "GT", "G-", "GN", "TA", "TC", "TG", "TT",
              "T-", "TN", "-A", "-C", "-G", "-T", "--", "-N", "NA", "NC", "NG",
              "NT", "N-", "NN")
.bMap                 <- matrix(character(256*2), ncol = 2)
colnames(.bMap)       <- c("read", "reference")
.bMap[]               <- NA
.bMap[.NUMBERS + 1, ] <- do.call(rbind, strsplit(.PAIRS, ""))

.bMapC <- function(idx) .bMap[idx + 1, ]

## 
## This reads a dataset which is structured like the alignments. It
## looks overly complicated because we read the file in an organized
## fashion to avoid reading the file too much
##
.getDatasetByIdxFast <- function(cmpH5, idx = seq.int(1, nrow(cmpH5)), dsName, convert = NULL) {
  if (missing(idx))
    idx <- seq.int(1, nrow(cmpH5))
  if (length(idx) == 0)
    stop("Index length must be > 0.")
  
  aln <- alnIndex(cmpH5)
  dsetNames <- base:::paste(aln$alnGroupPath[idx], dsName, sep = '/')
  
  udsets <- dsetNames[!duplicated(dsetNames)]
  dsets <- lapply(udsets, function(a) getH5Dataset(cmpH5, a))
  names(dsets) <- udsets
  dfactor <- factor(dsetNames, udsets)
  sidxs <- split(idx, dfactor)
  iidxs <- split(seq.int(length(idx)), dfactor)
  res <- vector("list", length(idx))
  
  for ( j in 1:length(dsets) ) {
    x <- read1DSlabs(dsets[[j]], aln$offsetBegin[sidxs[[j]]], aln$offsetEnd[sidxs[[j]]] -
                     aln$offsetBegin[sidxs[[j]]] + 1)
    if (! is.null(convert)) {
      x <- lapply(x, convert)
    }
    res[iidxs[[j]]] <- x
  }
  return(res)
}


##
## This was the old function which still works for datasets with
## dimensions, hence it is preserved.
## 
.getDatasetByIdx <- function(cmpH5, idx = seq.int(1, nrow(cmpH5)), dsName, convert = NULL) {
  if (missing(idx))
    idx <- seq.int(1, nrow(cmpH5))
  if (length(idx) == 0)
    stop("Index length must be > 0.")

  
  aln <- alnIndex(cmpH5)[idx, c("alnGroupPath", "offsetBegin", "offsetEnd"), drop = FALSE]

  rgp <- aln[,1]
  sqa <- seq.int(1, nrow(aln))
  rgp <- factor(rgp, rgp[!duplicated(rgp)])
  lst <- split(sqa, rgp)

  a <- lapply(lst, function(j) {
    gg <- getH5Group(cmpH5, aln[j[1], "alnGroupPath"])

    if (! h5DatasetExists(gg, dsName)) {
      stop(sprintf("dataset: %s doesn't exist.", dsName))
    }
    aa <- getH5Dataset(gg, dsName)

    if (is.null(dim(aa))) {
      v <- read1DSlabs(aa, aln$offsetBegin[j], aln$offsetEnd[j] - aln$offsetBegin[j] + 1)
      if (! is.null(convert)) {
        lapply(v, convert)
      } else {
        return ( v )
      }
    } else {
      lapply(j, function(i) {
        hs <- hSlab(c(aln[i,"offsetBegin"], 1), end = c(aln[i,"offsetEnd"], ncol(aa)))
        if (is.null(convert))
          aa[hs]
        else 
          convert(aa[hs])
      })
    }
  })

  ## make sure that the order of the results is same as the idx.
  a <- do.call(c, a)
  b <- character(length(a))
  b[do.call(c, lst)] <- names(a)
  a <- a[b]
  names(a) <- NULL
  
  return(a)
}

getAlignmentsRaw <- function(cmpH5, idx) .getDatasetByIdxFast(cmpH5, idx, "AlnArray")
getIPD <- function(cmpH5, idx) .getDatasetByIdxFast(cmpH5, idx, "IPD")
getPulseWidth <- function(cmpH5, idx) .getDatasetByIdxFast(cmpH5, idx, "PulseWidth")
getStartTime <- function(cmpH5, idx) .getDatasetByIdxFast(cmpH5, idx, "StartTime")
getPkmid <- function(cmpH5, idx) .getDatasetByIdxFast(cmpH5, idx, "pkmid")

setMethod("getAlignments", "PacBioCmpH5", function(h5Obj, idx) {
  .getDatasetByIdxFast(h5Obj, idx, "AlnArray", convert = .bMapC)
})

setMethod("getQualityValue", "PacBioCmpH5", function(h5Obj, idx) {
  .getDatasetByIdxFast(h5Obj, idx, "QualityValue")
})

getEviconsCalls <- function(cmpH5, refSeq, collapse = FALSE) {
  refSeq <- getRefPath(cmpH5, refSeq)
  s <- paste(refSeq, "Consensus/ConsensusCalls", sep = "/")
  if (! h5DatasetExists(cmpH5, s)) {
    stop("Consensus Information not stored in file.")
  }
  a <- c("A","C","G","T")[match(getH5Dataset(cmpH5, s)[],
                                c(1, 2, 4, 8))]
  if (collapse)
    paste(a, collapse = "")
  else
    a
}

getEviconsConfidence <- function(cmpH5, refSeq) {
  refSeq <- getRefPath(cmpH5, refSeq)
  s <- paste(refSeq, "Consensus/ConsensusConfidence", sep = "/")
  if (! h5DatasetExists(cmpH5, s)) {
    stop("Consensus Information not stored in file.")
  }
  getH5Dataset(cmpH5, s)[]
}

## #############################################################################
##
## Base H5 Files
##
## #############################################################################
.ofBasecalls <- c("Basecall", "DeletionQV", "DeletionTag",
                  "InsertionQV", "PreBaseFrames", "PulseIndex",
                  "QualityValue", "SubstitutionQV", "SubstitutionTag",
                  "WidthInFrames")

setMethod("getMovieName", "PacBioBasH5", function(h5Obj) {
  getH5Attribute(getH5Group(h5Obj, "ScanData/RunInfo"), "MovieName")[]
})

setMethod("initialize", "PacBioBasH5", function(.Object, fileName = NULL) {
  .Object <- callNextMethod(.Object, fileName = fileName)

  baseCalls <- getH5Group(.Object, "PulseData/BaseCalls")
  bcZMW <- do.call(cbind, lapply(c("HoleNumber", "NumEvent", "HoleXY"), function(nm) {
    getH5Dataset(getH5Group(baseCalls, "ZMW"), nm)[]
  }))
  bcZMW <- cbind(bcZMW, "offset" = cumsum(c(1, bcZMW[-nrow(bcZMW), 2])))
  colnames(bcZMW) <- c("holeNumber", "numEvent", "x", "y", "offset")
  
  .Object@baseEvents  <- bcZMW
  .Object@baseCallsG  <- baseCalls

  return(.Object)
})

PacBioBasH5 <- function(fileName) {
  new("PacBioBasH5", fileName = fileName)
}

## #############################################################################
##
## Pulse H5 Files
##
## #############################################################################
.ofPulses <- c("Channel", "Chi2", "ClassifierQV", "IsPulse", "MaxSignal",
               "MeanSignal", "MidSignal", "MidStdDev",
               "StartFrame", "WidthInFrames")

PacBioPlsH5 <- function(fileName) {
  obj <- new("PacBioPlsH5", fileName = fileName)
  
  pulseCalls <- getH5Group(obj, "PulseData/PulseCalls")
  
  pcZMW <- do.call(cbind, lapply(c("HoleNumber", "NumEvent", "HoleXY"), function(nm) {
    getH5Dataset(getH5Group(pulseCalls, "ZMW"), nm)[]
  }))
  pcZMW <- cbind(pcZMW, "offset" = cumsum(c(1, pcZMW[-nrow(pcZMW), 2])))

  colnames(pcZMW) <- c("holeNumber", "numEvent", "x", "y", "offset")
  
  obj@pulseEvents <- pcZMW
  obj@pulseCallsG <- pulseCalls
  
  return(obj)
}

getPulseEvents <- function(plsH5) plsH5@pulseEvents
getBaseEvents  <- function(basH5) basH5@baseEvents
getHoleNumbers <- function(basH5) basH5@baseEvents[, "holeNumber"]

getBaselineSigma <- function(plsH5) {
  getH5Dataset(plsH5, "PulseData/PulseCalls/ZMW/BaselineSigma")
}

.getFromPlsH5 <- function(plsH5, plsH5GName = c("BaseCalls", "PulseCalls"), dsName,
                          convert = NULL, holeNumbers = NULL) {
  grp <- switch(match.arg(plsH5GName), "BaseCalls" = plsH5@baseCallsG,
                "PulseCalls" = plsH5@pulseCallsG)
  evt <- switch(match.arg(plsH5GName), "BaseCalls" = getBaseEvents(plsH5),
                "PulseCalls" = getPulseEvents(plsH5))
  d <- getH5Dataset(grp, dsName)
  
  if (missing(holeNumbers) || is.null(holeNumbers)) {
    rows <- seq.int(1, nrow(evt))
    
    ## cache the dataset if you are to read it all.
    d <- d[]
  }
  else {
    rows <- match(holeNumbers, evt[,"holeNumber"])
  }
  
  if (any(is.na(rows)) || length(rows) <= 0)
    stop("Invalid ZMW: ", paste(holeNumbers, collapse = ", "))
  
  isV <- if (is.null(dim(d))) TRUE else FALSE
  
  v <- lapply(rows, function(r) {
    s <- evt[r, c("offset", "numEvent")]
    if (s[2] == 0)
      return(NULL)
    
    i <- seq.int(s[1], s[1] + s[2] - 1)
    
    if (is.null(convert)) {
      if (isV) d[i] else d[i,,drop=FALSE]
    }
    else {
      convert(if (isV) d[i] else d[i,,drop=FALSE])
    }
  })
  if (is.null(holeNumbers))
    names(v) <- evt[rows, "holeNumber"]
  else
    names(v) <- holeNumbers
  
  return(v)
}

.convertToLettersFromAscii <- function(a) {
  m <- c(65, 67, 71, 84, 45, 97, 99, 103, 116)
  b <- c("A", "C", "G", "T", "-", "A", "C", "G", "T")
  b[match(a, m)]
}

getBasecalls <- function(basH5, holeNumbers = getHoleNumbers(basH5))
  .getFromPlsH5(basH5, "BaseCalls", "Basecall", convert = .convertToLettersFromAscii,
                holeNumbers = holeNumbers)

getDeletionQV <- function(basH5, holeNumbers = getHoleNumbers(basH5))
  .getFromPlsH5(basH5, "BaseCalls", "DeletionQV", holeNumbers = holeNumbers)

getDeletionTag <- function(basH5, holeNumbers = getHoleNumbers(basH5))
  .getFromPlsH5(basH5, "BaseCalls", "DeletionTag", convert = .convertToLettersFromAscii,
                holeNumbers = holeNumbers)

getInsertionQV <- function(basH5, holeNumbers = getHoleNumbers(basH5))
  .getFromPlsH5(basH5, "BaseCalls", "InsertionQV", holeNumbers = holeNumbers)

getPreBaseFrames <- function(basH5, holeNumbers = getHoleNumbers(basH5))
  .getFromPlsH5(basH5, "BaseCalls", "PreBaseFrames", holeNumbers = holeNumbers)

getPulseIndex <- function(basH5, holeNumbers = getHoleNumbers(basH5))
  .getFromPlsH5(basH5, "BaseCalls", "PulseIndex", holeNumbers = holeNumbers)

getSubstitutionQV <- function(basH5, holeNumbers = getHoleNumbers(basH5))
  .getFromPlsH5(basH5, "BaseCalls", "SubstitutionQV", holeNumbers = holeNumbers)

getSubstitutionTag <- function(basH5, holeNumbers = getHoleNumbers(basH5))
  .getFromPlsH5(basH5, "BaseCalls", "SubstitutionTag", holeNumbers = holeNumbers)

setMethod("getQualityValue", "PacBioPlsH5", function(h5Obj, ...) {
  .getFromPlsH5(h5Obj, "BaseCalls", "QualityValue", ...)
})

setMethod("getWidthInFrames", "PacBioBasH5", function(h5Obj, ...) {
  .getFromPlsH5(h5Obj, "BaseCalls", "WidthInFrames", ...)
})

getChannel <- function(plsH5, holeNumbers = getHoleNumbers(plsH5))
  .getFromPlsH5(plsH5, "PulseCalls", "Channel", holeNumbers = holeNumbers)
getChi2 <- function(plsH5, holeNumbers = getHoleNumbers(plsH5))
  .getFromPlsH5(plsH5, "PulseCalls", "Chi2", holeNumbers = holeNumbers)
getClassifierQV <- function(plsH5, holeNumbers = getHoleNumbers(plsH5))
  .getFromPlsH5(plsH5, "PulseCalls", "ClassifierQV", holeNumbers = holeNumbers)
getIsPulse <- function(plsH5, holeNumbers = getHoleNumbers(plsH5))
  .getFromPlsH5(plsH5, "PulseCalls", "IsPulse",  holeNumbers = holeNumbers)
getMaxSignal <- function(plsH5, holeNumbers = getHoleNumbers(plsH5))
  .getFromPlsH5(plsH5, "PulseCalls", "MaxSignal", holeNumbers = holeNumbers)
getMeanSignal <- function(plsH5, holeNumbers = getHoleNumbers(plsH5))
  .getFromPlsH5(plsH5, "PulseCalls", "MeanSignal", holeNumbers = holeNumbers)
getMidSignal <- function(plsH5, holeNumbers = getHoleNumbers(plsH5))
  .getFromPlsH5(plsH5, "PulseCalls", "MidSignal", holeNumbers = holeNumbers)
getMidStdDevSignal <- function(plsH5, holeNumbers = getHoleNumbers(plsH5))
  .getFromPlsH5(plsH5, "PulseCalls", "MidStdDevSignal", holeNumbers = holeNumbers)
getStartFrame <- function(plsH5, holeNumbers = getHoleNumbers(plsH5))
  .getFromPlsH5(plsH5, "PulseCalls", "StartFrame", holeNumbers = holeNumbers)

setMethod("getWidthInFrames", "PacBioPlsH5", function(h5Obj, ...) {
  .getFromPlsH5(h5Obj, "PulseCalls", "WidthInFrames", ...)
})


## ############################################################################
##
## Align H5 Files
##
## #############################################################################
PacBioAlnH5 <- function(fileName, whAlignment = c("Alignments", "CRFAlignments",
                                    "GlobalAlignments")) {
  obj <- new("PacBioAlnH5", fileName = fileName)
  whAlignment <- match.arg(whAlignment)
    
  ## initialize groups.
  alignments <- getH5Group(obj, whAlignment)
  
  ZMW <- sapply(c("HoleNumber", "ReadLength", "ReadStartBase",
                  "ReadString.Index", "TemplateLength", "TemplateStartBase",
                  "TemplateString.Index"), function(nm) {
                    getH5Dataset(alignments, nm)[]
                  })
  cn <- colnames(ZMW)
  ZMW <- cbind(ZMW, getH5Dataset(alignments, "Accuracy")[])
  colnames(ZMW) <- c(cn, "Accuracy")

  ZMW <- ZMW[order(ZMW[,"HoleNumber"]), ]
    
  obj@ZMWs <- ZMW 
  obj@alignment <- alignments
  return(obj)
}

getZMWs <- function(alnH5) {
  alnH5@ZMWs
}

setMethod("nrow", "PacBioAlnH5", function(x) {
  nrow(x@ZMWs)
})

.getAlnH5String <- function(alnH5, zmws, wh = c("ReadString", "TemplateString"),
                            convert = .convertToLettersFromAscii) {
  wh   <- match.arg(wh)
  whS  <- paste(wh, "Index", sep = ".")
  cs   <- cbind(start  = cumsum(c(1, getZMWs(alnH5)[-nrow(alnH5), whS])),
                offset = getZMWs(alnH5)[, whS])

  if (missing(zmws))
    zmws <- getZMWs(alnH5)[, "HoleNumber"]
  
  mIdx <- match(zmws, getZMWs(alnH5)[,"HoleNumber"])
  cs <- cs[mIdx,]
  rownames(cs) <- zmws
  cs <- cs[cs[,2] > 0, ]
  cs[,2] <- cs[,2]
  
  ds <- getH5Dataset(alnH5@alignment, wh)
  
  apply(cs, 1, function(slc) {
    convert(ds[hSlab(slc[1], width = slc[2])])
  })
}

setMethod("getAlignments", "PacBioAlnH5", function(h5Obj, zmws, ...) {
  mapply(function(r, t) {
    cbind("read" = r, "reference" = t)
  }, .getAlnH5String(h5Obj, zmws = zmws, wh = "ReadString", ...),
         .getAlnH5String(h5Obj, zmws = zmws, wh = "TemplateString", ...),
         SIMPLIFY = FALSE)
})

## #############################################################################
##
## Trace H5 Files
##
## #############################################################################
PacBioTrcH5 <- function(fileName) {
  obj <- new("PacBioTrcH5", fileName = fileName)  
  return(obj)
}

getTraces <- function(trcH5) {
  getH5Dataset(trcH5, "TraceData/Traces", inMemory = FALSE)
}

