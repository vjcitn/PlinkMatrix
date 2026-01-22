# subsetByOverlaps(RangedPlink, gr) will subset to variants overlapping gr
# other indexing blocked

setClass("RangedPlink", slots=c(plinkCalls="DelayedMatrix",
    ranges="GRanges"))

#' concise presentation
#' @param object RangedPlink instance
#' @export
setMethod("show", "RangedPlink", function(object) {
 cat(sprintf("RangedPlink with %d variants on %d samples.\n", 
     ncol(slot(object, "plinkCalls")), nrow(slot(object, "plinkCalls"))))
 cat("  some sample ids:\n  ")
 cat(selectSome(rownames(slot(object, "plinkCalls"))), "\n")
 if (length(nn <- names(slot(object, "ranges")))>0) {
  cat("  some range ids:\n  ")
  cat(selectSome(nn), "\n")
  }
})



#' constructor for RangedPlink
#' @param delmat DelayedMatrix as returned by PlinkArray
#' @param ranges GenomicRanges checked for conformity with array
#' @examples
#' ex = example_PlinkArray()
#' data("example_GRanges", package="PlinkArray")
#' lit = RangedPlink(ex, example_GRanges)
#' lit
#' myGR = GenomicRanges::GRanges("chr18:1-500000")
#' IRanges::subsetByOverlaps(lit, myGR)
#' @export
RangedPlink = function(delmat, ranges) {
  stopifnot(length(ranges) == ncol(delmat))
  new("RangedPlink", plinkCalls=delmat, ranges=ranges)
}

#' implement a form of subsetting
#' @import IRanges
#' @import S4Vectors
#' @param x instance of RangedPlink
#' @param ranges instance of GenomicRanges
#' @param maxgap integer, see subsetByOverlaps, defaults to -1L
#' @param minoverlap integer, see subsetByOverlaps, defaults to 0L
#' @param type character, see subsetByOverlaps, defaults to 'any'
#' @param invert, logical, see subsetByOverlaps
#' @param \dots not used
#' @return RangedPlink instance
#' @export
setMethod("subsetByOverlaps", c("RangedPlink", "GRanges"), function (x, ranges, maxgap = -1L, minoverlap = 0L, type = c("any", 
    "start", "end", "within", "equal"), invert, ...) {
  if (!missing("invert")) stop("invert parameter not currently supported")
  ov = findOverlaps(ranges, slot(x, "ranges"), maxgap=maxgap, minoverlap=minoverlap, type=type)
  ansrngs = slot(x, "ranges")[subjectHits(ov)]
  ansarr = slot(x, "plinkCalls")[,subjectHits(ov)]
  RangedPlink(ansarr, ansrngs)
})
