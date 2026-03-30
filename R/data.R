#' sample characteristics of 445 GEUVADIS samples
#' @importFrom utils data
#' @docType data
#' @format data.frame
#' @examples
#' data("g445samples", package="PlinkMatrix")
#' g445samples[seq_len(4), seq_len(4)]
#' @note Example data are those provided with tensorqtl, see \url{https://github.com/broadinstitute/tensorqtl/tree/0c4db65a0cdc47f3b824ae530b89d270ef5e0096/example/data}.
#' @usage data("g445samples", package="PlinkMatrix")
"g445samples"

#' sample GRanges coordinated with example_PlinkMatrix
#' @docType data
#' @format GRanges
#' @examples
#' data("example_GRanges", package="PlinkMatrix")
#' head(example_GRanges)
#' @usage data("example_GRanges", package="PlinkMatrix")
"example_GRanges"
