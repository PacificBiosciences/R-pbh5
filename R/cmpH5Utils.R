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

#############################################################################
##
## Accessors
##
#############################################################################
getReadLength <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  alnIndex(cmpH5)[idx, "rEnd"] - alnIndex(cmpH5)[idx, "rStart"] + 1
}

getAlignedLength <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  alnIndex(cmpH5)[idx, "offsetEnd"] - alnIndex(cmpH5)[idx, "offsetBegin"] + 1
}

getZScore <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  if (h5DatasetExists(cmpH5, "AlnInfo/ZScore"))
    getH5Dataset(cmpH5, "AlnInfo/ZScore")[idx]
  else
    stop("cmpH5 does not contain ZScore dataset.")
}

getFullRefNames <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  cmpH5$fullRefName[idx]
}

getRefNames <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  cmpH5$refName[idx]
}

getAccuracy <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  1 - rowSums(alnIndex(cmpH5)[idx, c("nMisMatches", "nInsertions", "nDeletions")])/
    getReadLength(cmpH5, idx)
}

getTemplateStrand <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  cmpH5$alignedStrand[idx]
}

getTemplateStart <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  cmpH5$tStart[idx]
}

getTemplateEnd <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  cmpH5$tEnd[idx]
}

getMachineName <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  sapply(strsplit(alnIndex(cmpH5)[idx, "movieName"], "_"), "[", 3)
}

getTemplateSpan <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  cmpH5$tEnd[idx] - cmpH5$tStart[idx] + 1
}

getAlignedLength <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  cmpH5$offsetEnd[idx] - cmpH5$offsetStart[idx] + 1
}

getAdvanceTime <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  .Call("PBR_advance_time", getCumulativeAdvanceTime(cmpH5, idx), getAlignmentsRaw(cmpH5, idx), PACKAGE = "pbh5")
}

getPolymerizationRate <- function(cmpH5, idx = seq.int(1, nrow(cmpH5))) {
  totalTime <- sapply(getCumulativeAdvanceTime(cmpH5, idx = idx), function(a) {
    a[length(a)] - a[1]
  })
  getTemplateSpan(cmpH5, idx)/totalTime
}

getLocalPolymerizationRate <- function(cmpH5, idx = seq.int(1, nrow(cmpH5)),
                                       binFunction = function(aln) {
                                         cut(seq.int(1, nrow(aln)), 10)
                                       }) {
 alns <- getAlignments(cmpH5, idx = idx)
  stimes <- getCumulativeAdvanceTime(cmpH5, idx = idx)
  mapply(function(aln, stime) {
    tapply(seq.int(1, nrow(aln)), binFunction(aln), function(i) {
      sum(aln[i,2] != '-')/(max(stime[i], na.rm = TRUE) - min(stime[i], na.rm = T))
    })
  }, alns, stimes, SIMPLIFY = FALSE)
}

getNumPasses <- function(cmpH5, idx = 1:nrow(cmpH5)) {
  stopifnot(h5DatasetExists(cmpH5, "AlnInfo/NumPasses"))
  getH5Dataset(cmpH5, "AlnInfo/NumPasses")[idx]
}

#############################################################################
##
## PulseH5 + CmpH5 access
##
#############################################################################
doWithPlsAndCmp <- function(cmpH5, plsH5s, fx, ..., SIMPLIFY = TRUE) {
  names(plsH5s) <- sapply(plsH5s, function(plsH5) {
    getH5Attribute(getH5Group(plsH5, "ScanData/RunInfo"), "MovieName")[]
  })
  stopifnot(all(getMovieName(cmpH5) %in% names(plsH5s)))

  plsFactor <- factor(getMovieName(cmpH5), names(plsH5s))
  idx <- seq.int(1, nrow(cmpH5))
  idxByMovie <- split(idx, plsFactor)
  
  lst <- mapply(function(idxs, plsH5) {
    if (length(idxs) <= 0)
      return(NULL)
    fx(cmpH5, plsH5, idxs)
  }, idxByMovie, plsH5s, SIMPLIFY = FALSE)

 
  ## this is necessary to produce a list in the same order as the
  ## original AlignmentIndex.
  if (SIMPLIFY) {
    lst <- Filter(Negate(is.null), lst)
    if (! is.null(dim(lst[[1]]))) {
      lst <- do.call(rbind, lst)
      rownames(lst) <- do.call(c, idxByMovie)
      lst <- lst[as.character(idx),,drop=FALSE]
    } else {
      lst <- do.call(c, lst)
      names(lst) <- do.call(c, idxByMovie)
      lst <- lst[as.character(idx)]
    }
  }
  return(lst)
}

getChannelToBaseMap <- function(plsH5) {
  sapply(0:3, function(d) {
    getH5Attribute(getH5Group(plsH5, sprintf("ScanData/DyeSet/Analog[%d]", d)), "Base")[]
  })
}

setMethod("getSNR", "PacBioCmpH5", function(h5Obj, plsH5s) {
  doWithPlsAndCmp(h5Obj, plsH5s, fx = function(cmpH5, plsH5, idxs) {
    snr <- getSNR(plsH5)
    snr[match(cmpH5$holeNumber[idxs], getHoleNumbers(plsH5)),,drop = FALSE]
  })
})

getFrameByBaselineSNR <- function(cmpH5, plsH5s, fx = function(x) median(x, na.rm = T)) {
  doWithPlsAndCmp(cmpH5, plsH5s, fx = function(cmpH5, plsH5, idxs) {
    bs <- getH5Dataset(plsH5, "PulseData/PulseCalls/ZMW/BaselineSigma")[]
    mtch <- match(cmpH5$holeNumber[idxs], plsH5@pulseEvents[,"holeNumber"])
    channels <- sapply(0:3, function(d) {
      getH5Attribute(getH5Group(plsH5, sprintf("ScanData/DyeSet/Analog[%d]", d)), "Base")[]
    })
    mapply(function(pk, cb, m) {
      channel <- match(cb[,1], channels)
      fx(pk/bs[m, channel])
    }, getPkmid(cmpH5, idxs), getAlignments(cmpH5, idxs), as.list(mtch), SIMPLIFY = FALSE)
  })
}

.narrowRegions <- function(regTable) {
  regTable <- regTable[order(regTable$holeNumber),]
  hqregs <- regTable$type == "HQRegion"
  starts <- regTable[hqregs, "start"]
  ends <- regTable[hqregs, "end"]
  rcycle <- table(regTable$holeNumber)
  starts <- rep(starts, rcycle)
  ends <- rep(ends, rcycle)

  msk <- !((regTable$start < starts & regTable$end < starts) | (regTable$start > ends & regTable$end > ends))
  regTable <- regTable[msk,]
  starts <- starts[msk]
  ends <- ends[msk]
  
  regTable[,"start"] <- ifelse(regTable[,"start"] < starts, starts, regTable[,"start"])
  regTable[,"end"]   <- ifelse(regTable[,"end"] > ends, ends, regTable[,"end"])
  return(regTable)
}

getMoleculeIndex <- function(cmpH5) {
  factor(interaction(cmpH5$movieName, cmpH5$holeNumber, drop = T))
}

getUnrolledTemplateSpan <- function(cmpH5, idx = 1:nrow(cmpH5), unique = TRUE) {
  x <- tapply(getTemplateSpan(cmpH5)[idx], mi <- getMoleculeIndex(cmpH5)[idx], sum)
  if (! unique) {
    x <- rep(x, tapply(idx, mi, length))
  }
  return(x)
}

#############################################################################
##
## kMer Code.
##
#############################################################################
.makeKmerLabelCache <- function(maxSize = 12, alphabet = c("A", "C", "G", "T")) {
  kMerLabelCache <- new.env(hash = T, parent = emptyenv())
  lapply(1:12, function(i) assign(as.character(i), NULL, kMerLabelCache))
  function(k) {
    key <- as.character(k)
    if (is.null(get(key, kMerLabelCache))) {
      M <- as.matrix(expand.grid(lapply(1:k, function(i) alphabet), stringsAsFactors = FALSE)[,k:1])
      colnames(M) <- sprintf("P%02d", 1:k)
      assign(key, list(M = M, S = apply(M, 1, base:::paste, collapse = "")), kMerLabelCache)
    }
    get(key, kMerLabelCache)
  }
}
.kMerLabelCache <- .makeKmerLabelCache()

.kMerTable <- function(lst, k = 2) {
  res <- .Call("PBR_k_mer_tabulator", lst, "ACGT", as.integer(k), PACKAGE = "pbh5")
  res <- matrix(res, nrow = length(lst), byrow = TRUE)
  colnames(res) <- .kMerLabelCache(k)$S
  return(res)
}

.kMerLabeler <- function(lst, k = 2, collapse = FALSE) {
  stopifnot(k > 1 || k > 12)
  res <- .Call("PBR_k_mer_labeler", lst, "ACGT", as.integer(k), PACKAGE = "pbh5")
  nms <- .kMerLabelCache(k)
  if (collapse) {
    lapply(res, function(a) nms$S[a+1])
  } else {
    lapply(res, function(a) nms$M[a+1,])
  }
}

getKmerContext <- function(cmpH5, idx = 1:nrow(cmpH5), up = 2, down = 2) {
  k <- up + down + 1
  alns <- getAlignments(cmpH5, idx)
  msks <- lapply(alns, function(a) a[, 2] != "-")
  alns <- mapply(msks, alns, FUN = function(m, d) d[m, 2], SIMPLIFY = FALSE)
  labels <- .kMerLabeler(alns, k = k, collapse = TRUE)
  mapply(labels, msks, FUN = function(a, b) {
    o <- rep(NA, length(b))
    o[b] <- a
    o
  }, SIMPLIFY = FALSE)
}

associateWithContext <- function(cmpH5, idx = 1:nrow(cmpH5), f = getIPD, up = 2, down = 2,
                                 useReference = TRUE, collapse = FALSE) {    
  k <- up + down + 1
  col <- if (useReference) 2 else 1
  hasDim <- !is.null(dim(f(cmpH5, idx[1])[[1]]))
  
  alns <- getAlignments(cmpH5, idx)
  msks <- lapply(alns, function(a) a[, col] != "-")
  alns <- mapply(msks, alns, FUN = function(m, d) d[m, col], SIMPLIFY = FALSE)
  elts <- mapply(msks, f(cmpH5, idx), FUN = function(m, d) if (hasDim) d[m,] else d[m],
                 SIMPLIFY = FALSE)
  labels <- .kMerLabeler(alns, k = k, collapse = collapse)
  r <- mapply(elts, labels, alns, FUN = function(elt, lab, aa) {
    r1 <- (up+1):((if (hasDim) nrow(elt) else length(elt)) - down)
    r2 <- 1:length(r1)
    list(elt = if (hasDim) elt[r1,] else elt[r1],
         label = if (collapse) lab[r2] else lab[r2,])
  }, SIMPLIFY = FALSE)
    
  data.frame(elt = do.call(if (hasDim) rbind else c, lapply(r, "[[", 1)),
             context = do.call(if (collapse) c else rbind, lapply(r, "[[", 2)),
             stringsAsFactors = FALSE)
}

summarizeByContext <- function(cmpH5, idx = 1:nrow(cmpH5), up = 2, down = 2,
                               statF = getIPD, summaryF = function(a) median(a, na.rm = T), ...) {
  ctx <- associateWithContext(cmpH5, up = up, down = down, idx = idx, f = statF, collapse = TRUE, ...)
  levels <- .kMerLabelCache(up + down + 1)$S
  fct <- factor(ctx[,"context"], levels = levels)
  dta <- data.frame(count = as.numeric(table(fct)),
                    value = tapply(ctx[,"elt"], fct, summaryF), stringsAsFactors = FALSE)
  rownames(dta) <- levels(fct)
  return(dta)
}

makeContextDataTable <- function(cmpH5, idx = 1:nrow(cmpH5), up = 2, down = 2,
                                 fxs = list(ipd = getIPD, pw = getPulseWidth, tpos = getTemplatePosition), ...) {
  f <- function(cmpH5, idxs) {
    lst <- lapply(fxs, function(f) {
      f(cmpH5, idxs)
    })
    lapply(1:length(idxs), function(i) {
      d <- do.call(base:::cbind, lapply(1:length(fxs), function(j) {
        lst[[j]][[i]]
      }))
      colnames(d) <- names(fxs)
      return(d)
    })
  }
  associateWithContext(cmpH5, idx, f, up = up, down = down, ...)
}

#############################################################################
##
## Positional Access.
##
#############################################################################
getReadsInRange <- function(cmpH5, refSeq, refStart = 1, refEnd = NA, idx = 1:nrow(cmpH5)) {
  refBlock <- .getRefBlock(cmpH5, refSeq)
  ranges <- alnIndex(cmpH5)[seq.int(refBlock[1], refBlock[2]), c("tStart", "tEnd", "nBackRead", "nOverlap")]
  if (is.na(refEnd)) refEnd <- getRefLength(cmpH5, refSeq)
  idxs <- .Call("PBR_indices_within_range", ranges[,1], ranges[,2], ranges[,3], ranges[,4],
                as.integer(refStart), as.integer(refEnd))
  ## XXX: This shouldn't happen, but it has.
  idxs[length(idxs)] <- ifelse(idxs[length(idxs)] >= refBlock[2], refBlock[2] - 1,
                               idxs[length(idxs)])
  idxs <- idxs + 1 + ifelse(0 == refBlock[1], 0, refBlock[1] - 1)
  return(intersect(idxs, idx))
}

getCoverageInRange <- function (cmpH5, refSeq, refStart = 1, refEnd = NA, idx = NULL) {
  refBlock <- .getRefBlock(cmpH5, refSeq)
  if (is.null(idx)) {
    idx <- seq.int(1, nrow(cmpH5))
  }
  w <- sort(intersect(idx, seq.int(refBlock[1], refBlock[2])))
  ranges <- alnIndex(cmpH5)[w, c("tStart", "tEnd", "nBackRead", "nOverlap")]
  if (is.na(refEnd)) refEnd <- getRefLength(cmpH5, refSeq)
  if (nrow(ranges) == 0) return(integer(refEnd - refStart + 1))
  .Call("PBR_coverage_within_range", ranges[, 1], ranges[,2], ranges[, 3], ranges[, 4],
        as.integer(refStart), as.integer(refEnd))
}

#############################################################################
##
## Template Position.
##
#############################################################################
templatePositionFromAln <- function(alnLst, strand, start, end) {
  .Call("PBR_get_template_position", alnLst, strand, start, end)
}

.reverse <- function(x, strand) {
  if (strand == 1) {
    if (! is.null(dim(x))) {
      x[nrow(x):1, ]
    } else {
      rev(x)
    }
  } else {
    x
  }
}

.reverseComplement <- function(x, strand) {
  if (strand == 0) {
    x
  } else {
    rev(c("T", "G", "C", "A", "-")[match(x, c("A", "C", "G", "T", "-"))])
  }
}

getTemplatePosition <- function(cmpH5, idx = 1:nrow(cmpH5), withAlignments = FALSE, asDataFrame = FALSE) {
  alns <- getAlignmentsRaw(cmpH5, idx)
  strands <- cmpH5$alignedStrand[idx]
  tpos <- templatePositionFromAln(alns, strands, cmpH5$tStart[idx], cmpH5$tEnd[idx])

  if (! withAlignments) {
    return(tpos)
  }
  mapply(function(a, strand, id, p) {
    ref  <- .reverse(.bMap[a+1, 2], strand)
    read <- .reverse(.bMap[a+1, 1], strand)
    z <- list(position = .reverse(p, strand), read = read, ref = ref,
              idx = rep(id, length(ref)), strand = rep(strand, length(ref)))
    if (asDataFrame) as.data.frame(z, stringsAsFactors = FALSE) else z
  }, alns, strands, idx, tpos, SIMPLIFY = FALSE)
}

getByTemplatePosition <- function(cmpH5, idx = 1:nrow(cmpH5), f = getIPD) {
  if (is.list(f)) {
    return(.getByTemplatePositionList(cmpH5, idx, f))
  }
  tpos <- getTemplatePosition(cmpH5, idx, withAlignments = TRUE)
  dta <- f(cmpH5, idx)
  strands <- getTemplateStrand(cmpH5, idx)
  x <- mapply(strands, dta, FUN = function(s, d) .reverse(d, s), SIMPLIFY = FALSE)
  x <- do.call(if (! is.null(dim(x[[1]]))) rbind else c, x)
  e0 <- tpos[[1]]
  c0 <- names(e0)
  d <- lapply(1:length(c0), function(i) do.call(c, lapply(tpos, function(e) e[[i]])))
  names(d) <- c0
  d$elt <- x
  as.data.frame(d, row.names = NULL, stringsAsFactors = FALSE)
}

.getByTemplatePositionList <- function (cmpH5, idx = 1:nrow(cmpH5),
                                        fxs = list(IPD = getIPD, PulseWidth = getPulseWidth, refName = getFullRefNames)) {
  tpos <- getTemplatePosition(cmpH5, idx, withAlignments = TRUE)
  strands <- getTemplateStrand(cmpH5, idx)
  dta <- lapply(fxs, function(f) {
    do.call(c, mapply(strands, tpos, f(cmpH5, idx), FUN = function(s, tp, d) {
      if (length(d) == 1) {
        ## to allow for the possibility of scalar vals per read.
        rep(d, length(tp[[1]]))
      } else {
        if (s == 0) d else d[length(d):1]
      }
    }, SIMPLIFY = FALSE))
  })
  e0 <- tpos[[1]]
  c0 <- names(e0)
  tpos <- lapply(1:length(c0), function(i) do.call(c, lapply(tpos, function(e) e[[i]])))
  names(tpos) <- c0
  as.data.frame.list(c(tpos, dta), row.names = NULL, stringsAsFactors = FALSE)
}


getAlignmentBlock <- function(cmpH5, ref, refStart, refEnd) {
  tpos <- getTemplatePosition(cmpH5, getReadsInRange(cmpH5, ref, refStart, refEnd), TRUE, TRUE)
  do.call(rbind, lapply(tpos, function(a) {
    w <- a$position %in% refStart:refEnd
    tapply(a$read[w], factor(a$position[w], refStart:refEnd), paste, collapse = "")
  }))
}


#############################################################################
##
## Consensus Code.
##
#############################################################################
getCallsForReads <- function(rawAlns, refStart, refEnd, strands, starts, ends) {
  stopifnot(length(rawAlns) == length(starts)  && length(rawAlns) == length(strands) &&
            length(rawAlns) == length(ends))
  .Call("PBR_compute_consensus", rawAlns, as.integer(refStart), as.integer(refEnd),
        as.integer(strands), as.integer(starts), as.integer(ends))
}

getConsensusForReads <- function(rawAlns, refStart, refEnd, strands, starts, ends, fast = TRUE) {
  if (fast) {
    stopifnot(length(rawAlns) == length(starts)  && length(rawAlns) == length(strands) &&
              length(rawAlns) == length(ends))
    .Call("PBR_compute_consensus2", rawAlns, as.integer(refStart), as.integer(refEnd),
          as.integer(strands), as.integer(starts), as.integer(ends))
  } else {
    computeConsensus(getCallsForReads(rawAlns, refStart, refEnd, strands, starts, ends))
  }
}

getConsensusForIdxs <- function(cmpH5, idx = 1:nrow(cmpH5)) {
  refNames <- getFullRefNames(cmpH5, idx)
  stopifnot(length(unique(refNames)) == 1)
  starts <- getTemplateStart(cmpH5, idx)
  ends <- getTemplateEnd(cmpH5, idx)
  strands <- getTemplateStrand(cmpH5, idx)
  getConsensusForReads(getAlignmentsRaw(cmpH5, idx), min(starts), max(ends), strands, starts, ends)
}

getConsensusInRange <- function(cmpH5, refSeq, refStart = 1, refEnd = getRefLength(cmpH5, refSeq), idx = 1:nrow(cmpH5),
                                blockSize = NA) {
  if (is.na(blockSize)) {
    idxs <- getReadsInRange(cmpH5, refSeq, refStart, refEnd, idx)
    getConsensusForReads(getAlignmentsRaw(cmpH5, idxs), refStart, refEnd,
                         getTemplateStrand(cmpH5, idxs), getTemplateStart(cmpH5, idxs),
                         getTemplateEnd(cmpH5, idxs))
  } else {
    intervals <- unique(c(seq.int(from=refStart, to=refEnd, by = blockSize),
                          refEnd+1))
    do.call(c, lapply(2:length(intervals), function(j) {
      idxs <- getReadsInRange(cmpH5, refSeq, refStart = intervals[j-1], refEnd = intervals[j]-1, idx)
      getConsensusForReads(getAlignmentsRaw(cmpH5, idxs), refStart=intervals[j-1],
                           refEnd = intervals[j]-1,
                           getTemplateStrand(cmpH5, idxs),
                           getTemplateStart(cmpH5, idxs),
                           getTemplateEnd(cmpH5, idxs))
    }))
  }
}

getCallsInRange <- function(cmpH5, refSeq, refStart = 1, refEnd = getRefLength(cmpH5, refSeq), idx = 1:nrow(cmpH5)) {
  idxs <- getReadsInRange(cmpH5, refSeq, refStart, refEnd, idx)
  getCallsForReads(getAlignmentsRaw(cmpH5, idxs), refStart, refEnd,
                   getTemplateStrand(cmpH5, idxs), getTemplateStart(cmpH5, idxs),
                   getTemplateEnd(cmpH5, idxs))
}

computeConsensus <- function(calls) {
  f <- base:::table # avoid dispatch and lookup.
  sapply(calls, function(a) {
    if (length(a) <= 0) return(NA)
    t <- f(a)
    m <- max(t) == t
    if (sum(m) == 1) {
      names(t)[m]
    } else { 
      sample(names(t)[m], size = 1)
    }
  })
}

.computeConsensusRle <- function(calls, reference) {
  runs <- rle(reference)
  rlabels <- rep(1:length(runs$length), runs$length)
  res <- tapply(calls, rlabels, function(a) {
    nbases <- round(sum(sapply(a, function(b) {
      if (is.null(b))
        0
      else
        mean(ifelse(b == '-', 0, ifelse(nchar(b) > 2, 2, 1)))
    })))
    
    calls <- do.call(c, a)
    calls <- calls[calls != '-']
    wbase <- names(sort(table(calls), decreasing = T))[1]
    bases <- rep(wbase, nbases)
    
    if (nbases < length(a)) {
      c(rep('-', length(a) - nbases), bases)
    } else if (nbases > length(a)) {
      if (length(a) == 1)
        paste(bases, collapse = '')
      else 
        c(bases[1:(length(a)-1)], paste(bases[length(a):length(bases)],
                                        collapse = ""))
    } else {
      bases
    }
  })
  as.character(do.call(c, res))
}

#############################################################################
##
## Access and Manipulation
##
#############################################################################
.getAlignmentsWithFeaturesAsDataFrame <- function(cmpH5, idx, features) {
   alns <- getAlignments(cmpH5, idx)
   ds <- lapply(features, function(fx) {
     do.call(c, fx(cmpH5, idx)) 
   })
   idx <- rep(idx, sapply(alns, nrow))
   alns <- do.call(rbind, alns)
   data.frame(alns, idx, ds, row.names = NULL, stringsAsFactors = FALSE)
}

getAlignmentsWithFeatures <- function(cmpH5, idx = seq.int(1, nrow(cmpH5)),
                                      fxs = list('ipd' = getIPD), collapse = FALSE) {
  if (collapse) {
    dta <- .getAlignmentsWithFeaturesAsDataFrame(cmpH5, idx, fxs)
  } else {
    alns <- getAlignments(cmpH5, idx)
    ds <- lapply(fxs, function(fx) {
      fx(cmpH5, idx)
    })
    ds <- lapply(1:length(idx), function(i) {
      do.call(cbind, lapply(ds, function(zds) zds[[i]]))
    })
    names(ds) <- names(alns)
    
    dta <- mapply(function(d, a) {
      data.frame(d, a, stringsAsFactors = FALSE)
    }, ds, alns, SIMPLIFY = FALSE)
  }
  return(dta)
}

#############################################################################
##
## Barcode Stuff
##
#############################################################################
getBarcodeLabels <- function(cmpH5, idx = 1:nrow(cmpH5)) {
  getH5Dataset(cmpH5, "AlnInfo/barcode")[idx]
}
getBarcodeScores <- function(cmpH5, idx = 1:nrow(cmpH5)) {
  getH5Dataset(cmpH5, "AlnInfo/barcodeScore")[idx]
}
