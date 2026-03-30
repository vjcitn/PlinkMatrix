#' abbreviated display
#' @keywords internal
selectSome <- function(obj, maxToShow = 5) {
  # from Biobase
  len <- length(obj)
  if (maxToShow < 3) {
    maxToShow <- 3
  }
  if (len > maxToShow) {
    maxToShow <- maxToShow - 1
    bot <- ceiling(maxToShow / 2)
    top <- len - (maxToShow - bot - 1)
    nms <- obj[c(seq_len(bot), top:len)]
    c(as.character(nms[seq_len(bot)]), "...", as.character(nms[-c(seq_len(bot))]))
  } else if (is.factor(obj)) {
    as.character(obj)
  } else {
    obj
  }
}
