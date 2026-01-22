#' S4 for bed reader
#' @import methods
#' @import DelayedArray
#' @export
setClass("PlinkPyArraySeed",
    contains="Array",
    slots=c(
        filepath="character",
        readerRef="ANY",  # will be the module with open_bed as method
        dim = "integer",
        dimnames = "list"
    )
)

#' seed
#' @param filepath character path to resource
#' @param readerRef reticulate-imported `bed_reader` module
#' @export
PlinkPyArraySeed <- function(filepath, readerRef) {
    filepath <- tools::file_path_as_absolute(filepath)
    stopifnot("open_bed" %in% names(readerRef))
    res = readerRef$open_bed(filepath)
    new("PlinkPyArraySeed", filepath=filepath, readerRef=readerRef,
       dim = as.integer(unlist(res$shape)), dimnames=list(res$iid, res$sid))
}


#' extractor
#' @import DelayedArray
#' @param x seed instance
#' @param index list of suitable values for extracting elements
#' @export
setMethod("extract_array", "PlinkPyArraySeed", function(x, index) {
      dims = slot(x, "dim")
      dimnames = slot(x, "dimnames")  # guaranteed
      stopifnot(length(index)==2)
      if (is.null(index[[1]])) i1 = 0:(dims[1]-1)  # should be efficient?
      if (is.null(index[[2]])) i2 = 0:(dims[2]-1)

      if (is.character(index[[1]])) {
          stopifnot(all(index[[1]] %in% dimnames[[1]]))
          i1 = match(index[[1]], dimnames[[1]])
          }
      else if (is.numeric(index[[1]])) {
          i1 = index[[1]]-1L  # python
          }

      if (is.character(index[[2]])) {
          stopifnot(all(index[[2]] %in% dimnames[[2]]))
          i2 = match(index[[2]], dimnames[[2]])
          }
      else if (is.numeric(index[[2]])) {
          i2 = index[[2]]-1L   # python
          }
      ans = slot(x, "readerRef")$open_bed(slot(x,"filepath"))$read(tuple(as.integer(i1), as.integer(i2)))
      ans
})

#' present seed concisely
#' @param object instance of PlinkPyArraySeed
setMethod("show", "PlinkPyArraySeed", function(object) {
 cat(sprintf("seed for Plink %d x %d\n", object@dim[1], object@dim[2]))
})

#' bed reader interface
#' @import DelayedArray
#' @import reticulate
#' @param path character path to resource
#' @examples
#' folder = tempdir()
#' pa = get_plink_example_path()
#' unzip(pa, exdir = folder)
#' tmp = dir(folder, full=TRUE)
#' ppath = grep("geuv445.bed", tmp, value=TRUE)
#' tst = PlinkPyArray(ppath)
#' tst
#' @export
PlinkPyArray = function(path) {
  reticulate::py_require("bed-reader")
  br = reticulate::import("bed_reader")
  DelayedArray(PlinkPyArraySeed(path, br))
}
