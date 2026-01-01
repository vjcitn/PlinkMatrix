#' Get sample metadata
#' @param x DelayedArray instance
#' @export
getSampleData <- function(x) {
  if (is(x, "DelayedArray")) {
    x <- seed(x)
  }
  x@fam
}

#' Get variant metadata
#' @param x DelayedArray instance
#' @export
getVariantData <- function(x) {
  if (is(x, "DelayedArray")) {
    x <- seed(x)
  }
  x@bim
}

## Get genotypes for specific samples and variants
#getGenotypes <- function(x, samples = NULL, variants = NULL) {
#  if (is.null(samples)) {
#    samples <- seq_len(nrow(x))
#  }
#  if (is.null(variants)) {
#    variants <- seq_len(ncol(x))
#  }
#  
#  x[samples, variants, drop = FALSE]
#}

