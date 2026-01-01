#' operate with BiocFileCache to retrieve zip file of plink example data
#' @import BiocFileCache
#' @param ca BiocFileCache instance
#' @examples
#' get_plink_example_path()
#' @export
get_plink_example_path = function(ca = BiocFileCache::BiocFileCache()) {
  url = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocPlinkArray/geuv445_plink.zip"
  inf = bfcquery(ca, "geuv445_plink.zip")
  nr = nrow(inf)
  if (nr == 0) {
    bfcadd(ca, fpath=url, rname="geuv445_plink.zip", download=TRUE)
  }
  inf = bfcquery(ca, "geuv445_plink.zip")
  inf$rpath
}
  
#' produce PlinkArray from example data
#' @param folder a path where unzipped example data will be managed
#' @examples
#' example_PlinkArray()
#' @export
example_PlinkArray = function(folder = tempdir()) {
  pa = get_plink_example_path()
  unzip(pa, exdir = folder)
  tmp = dir(folder, full=TRUE)
  remsuff = function(x) sub(".bed$", "", x)
  ppath = remsuff(grep("geuv445.bed", tmp, value=TRUE))
  PlinkArray(ppath)
}
