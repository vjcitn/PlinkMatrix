#' produce GRanges from variant notation for plink example from geuvadis
#' @import GenomicRanges
#' @import IRanges
#' @param x character vector of variant names
#' @param sepused single character, defaults to "_"
#' @return GRanges instance
#' @examples
#' plid2gr("chr18_80259028_AG_A_b38")
#' @export
plid2gr = function (x, sepused="_") 
{
    ss = strsplit(x, sepused)
    pos = as.integer(sapply(ss, "[", 2))
    ch = vapply(ss, function(x) x[1], character(1))
    ans = GRanges(ch, IRanges(pos, width = 1))
    names(ans) = x
    ans
}

