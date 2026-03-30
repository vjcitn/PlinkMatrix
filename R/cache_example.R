#' operate with BiocFileCache to retrieve zip file of plink example data
#' @import BiocFileCache
#' @param ca BiocFileCache instance
#' @return character string with path to example in cache
#' @examples
#' get_plink_example_path()
#' @export
get_plink_example_path <- function(ca = BiocFileCache::BiocFileCache()) {
  url <- "https://mghp.osn.xsede.org/bir190004-bucket01/BiocPlinkMatrix/geuv445_plink.zip"
  inf <- bfcquery(ca, "geuv445_plink.zip")
  nr <- nrow(inf)
  if (nr == 0) {
    bfcadd(ca, fpath = url, rname = "geuv445_plink.zip", download = TRUE)
  }
  inf <- bfcquery(ca, "geuv445_plink.zip")
  inf$rpath
}

#' produce PlinkMatrix from example data
#' @import SummarizedExperiment
#' @param folder a path where unzipped example data will be managed
#' @param as_RSE logical(1) if TRUE (default is FALSE) a RangedSummarizedExperiment is returned
#' @return a SummarizedExperiment or RangedSummarizedExperiment
#' with rowRanges calculated from SNP ids
#' @examples
#' example_PlinkMatrix()
#' @export
example_PlinkMatrix <- function(folder = tempdir(), as_RSE = FALSE) {
  data("g445samples", package = "PlinkMatrix")
  data("example_GRanges", package = "PlinkMatrix")
  pa <- get_plink_example_path()
  unzip(pa, exdir = folder)
  tmp <- dir(folder, full.names = TRUE)
  remsuff <- function(x) sub(".bed$", "", x)
  ppath <- remsuff(grep("geuv445.bed", tmp, value = TRUE))
  ans <- PlinkMatrix(ppath)
  if (as_RSE) {
    #   cn = colnames(ans) |> strsplit("_")
    colnames(ans) <- gsub("0_", "", colnames(ans))
    ans <- SummarizedExperiment(list(calls = ans), rowData = example_GRanges, colData = g445samples)
  }
  ans
}
