#' browse and launch the learnr tutorials bundled with PlinkMatrix
#' @details Lists the \code{learnr} tutorials shipped with this package --
#' covering matrix foundations, the Singular Value Decomposition, Principal
#' Components Analysis, and ancestry discrimination using the package's
#' example genotype data. With \code{name} supplied, that tutorial is
#' launched directly. With \code{name} missing in an interactive session, a
#' menu is shown with a checkbox beside each tutorial already marked done;
#' the menu keeps returning after each tutorial finishes (offering to mark
#' it done) until you choose 0 to quit. Progress is kept in a small file in
#' your own R user data directory (see \code{\link[tools]{R_user_dir}}), so
#' it survives restarting R; it is self-reported and is not sent anywhere.
#' Each tutorial has its own "Exit tutorial" button that cleanly stops the
#' tutorial's web server and returns control to the R console, so there is
#' no need to interrupt with Ctrl-C. Requires the \code{learnr} package.
#' @param name character(1) name of a tutorial to launch directly (the
#' tutorial's directory name under \code{inst/tutorials}, as reported by
#' \code{learnr::available_tutorials(package = "PlinkMatrix")}). If missing,
#' a menu is presented in interactive sessions; in non-interactive sessions
#' the available tutorial names are printed and an error is raised.
#' @return invisible(NULL) if \code{name} is supplied; otherwise, invisibly,
#' a named logical vector recording which tutorials are marked done. As a
#' side effect, launches one or more \code{learnr} tutorials, typically in
#' the default web browser, via \code{learnr::run_tutorial}.
#' @importFrom utils menu askYesNo
#' @importFrom tools R_user_dir
#' @examples
#' if (interactive()) PlinkMatrix_tutorials()
#' @export
PlinkMatrix_tutorials <- function(name) {
  if (!requireNamespace("learnr", quietly = TRUE)) {
    stop(
      "Package 'learnr' is required to run PlinkMatrix tutorials. ",
      "Install it with install.packages('learnr')."
    )
  }
  avail <- as.data.frame(learnr::available_tutorials(package = "PlinkMatrix"))
  if (nrow(avail) == 0) {
    stop("No tutorials found for package PlinkMatrix.")
  }

  if (!missing(name) && !is.null(name)) {
    if (!name %in% avail$name) {
      stop(
        "Unknown tutorial name: '", name, "'. Available: ",
        paste(avail$name, collapse = ", ")
      )
    }
    learnr::run_tutorial(name, package = "PlinkMatrix")
    return(invisible(NULL))
  }

  if (!interactive()) {
    print(avail[, c("name", "title")])
    stop(
      "Session is not interactive; specify `name` explicitly. ",
      "See the printed table above for available tutorial names."
    )
  }

  repeat {
    done <- tutorial_progress_read(avail$name)
    choices <- c(
      sprintf("[%s] %s: %s", ifelse(done, "x", " "), avail$name, avail$title),
      "Clear my progress record"
    )
    title <- sprintf(
      "PlinkMatrix tutorials (%d of %d done). Pick one, or 0 to quit:",
      sum(done), length(done)
    )
    pick <- menu(choices, title = title)
    if (pick == 0L) {
      break
    }
    if (pick == length(choices)) {
      if (isTRUE(askYesNo("Clear all progress ticks?", default = FALSE))) {
        tutorial_progress_write(done & FALSE)
      }
      next
    }
    tname <- avail$name[pick]
    learnr::run_tutorial(tname, package = "PlinkMatrix")
    if (!done[[tname]] &&
      isTRUE(askYesNo(sprintf("Mark tutorial '%s' as done?", tname), default = FALSE))) {
      done[[tname]] <- TRUE
      tutorial_progress_write(done)
    }
  }
  invisible(tutorial_progress_read(avail$name))
}

tutorial_progress_file <- function() {
  getOption(
    "PlinkMatrix.progress_file",
    file.path(R_user_dir("PlinkMatrix", "data"), "tutorial-progress.rds")
  )
}

tutorial_progress_read <- function(names) {
  done <- rep(FALSE, length(names))
  names(done) <- names
  f <- tutorial_progress_file()
  if (file.exists(f)) {
    saved <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.character(saved)) {
      done[intersect(saved, names(done))] <- TRUE
    }
  }
  done
}

tutorial_progress_write <- function(done) {
  f <- tutorial_progress_file()
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(names(done)[done], f)
  invisible(done)
}
