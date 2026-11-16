# Warmup tutorials: run one locally, or pick from a menu that remembers progress.
# Not owned by any team.  After editing the roxygen comments run devtools::document().

tnames <- c(a = "warmup-a-beta", b = "warmup-b-simulation", c = "warmup-c-intervals",
            d = "warmup-d-loss-shortest", e = "warmup-e-densities", f = "warmup-f-ggplot",
            g = "warmup-g-tests-prs")

wtitles <- c(a = "Beta distributions as beliefs",
             b = "From Se, Sp and prevalence to the PPV",
             c = "Intervals: equal-tailed and shortest",
             d = "Mean, median, mode: different questions",
             e = "Densities, MAP and HPD regions",
             f = "Building a figure in ggplot2 layers",
             g = "Tests and the pull-request workflow")

#' Run a warmup tutorial locally
#'
#' @param which A single letter `"a"` to `"g"` selecting the warmup
#'   (see `screenbayes:::tnames`).
#' @return Called for its side effect of starting the tutorial. It returns when you stop it.
#' @examples
#' screenbayes:::tnames
#' @export
run_tut <- function(which = "a") {
  stopifnot(length(which) == 1, which %in% names(tnames))
  for (pkg in c("shiny", "learnr"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("please install the ", pkg, " package: install.packages(\"", pkg, "\")",
           call. = FALSE)
  learnr::run_tutorial(tnames[[which]], package = "screenbayes")
}

# ---- progress record: a small file in the user's own data directory -------------

progress_file <- function() {
  getOption("screenbayes.progress_file",
            file.path(tools::R_user_dir("screenbayes", "data"), "warmup-progress.rds"))
}

#' Which warmups have been marked done?
#' @return A named logical vector over `"a"` to `"g"`.
#' @keywords internal
read_progress <- function() {
  done <- stats::setNames(rep(FALSE, length(tnames)), names(tnames))
  f <- progress_file()
  if (file.exists(f)) {
    saved <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.character(saved)) done[intersect(saved, names(done))] <- TRUE
  }
  done
}

write_progress <- function(done) {
  f <- progress_file()
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(names(done)[done], f)
  invisible(done)
}

is_interactive <- function() interactive()

warmup_labels <- function(done) {
  sprintf("[%s] %s. %s", ifelse(done, "x", " "), toupper(names(tnames)), wtitles[names(tnames)])
}

#' Pick a warmup from a menu
#'
#' Shows the seven warmups with a tick beside each one you have marked done, runs the
#' one you choose with [run_tut()], and when you stop it asks whether to mark it done.
#' The menu comes back until you choose 0. Progress is kept in a small file in your own
#' R user data directory, so it survives restarting R. It is self-reported and is not
#' sent anywhere.
#'
#' @param menu_fun,ask_fun,run_fun Replaceable for testing; leave at their defaults.
#' @return Invisibly, a named logical vector of the warmups marked done.
#' @examples
#' if (interactive()) warmups()
#' @export
warmups <- function(menu_fun = utils::menu, ask_fun = utils::askYesNo, run_fun = run_tut) {
  if (!is_interactive())
    stop("warmups() is a menu and needs an interactive session; use run_tut(\"a\") etc.",
         call. = FALSE)
  repeat {
    done <- read_progress()
    choices <- c(warmup_labels(done), "Clear my progress record")
    title <- sprintf("Warmups (%d of %d done). Pick one, or 0 to quit:", sum(done), length(done))
    k <- menu_fun(choices, title = title)
    if (k == 0L) break
    if (k == length(choices)) {
      if (isTRUE(ask_fun("Clear all ticks?", default = FALSE))) write_progress(done & FALSE)
      next
    }
    w <- names(tnames)[k]
    run_fun(w)
    if (!done[[w]] && isTRUE(ask_fun(sprintf("Mark warmup %s as done?", toupper(w)), default = FALSE))) {
      done[[w]] <- TRUE
      write_progress(done)
    }
  }
  invisible(read_progress())
}
