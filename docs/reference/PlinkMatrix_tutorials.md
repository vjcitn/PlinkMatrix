# browse and launch the learnr tutorials bundled with PlinkMatrix

browse and launch the learnr tutorials bundled with PlinkMatrix

## Usage

``` r
PlinkMatrix_tutorials(name)
```

## Arguments

- name:

  character(1) name of a tutorial to launch directly (the tutorial's
  directory name under `inst/tutorials`, as reported by
  `learnr::available_tutorials(package = "PlinkMatrix")`). If missing, a
  menu is presented in interactive sessions; in non-interactive sessions
  the available tutorial names are printed and an error is raised.

## Value

invisible(NULL) if `name` is supplied; otherwise, invisibly, a named
logical vector recording which tutorials are marked done. As a side
effect, launches one or more `learnr` tutorials, typically in the
default web browser, via
[`learnr::run_tutorial`](https://pkgs.rstudio.com/learnr/reference/run_tutorial.html).

## Details

Lists the `learnr` tutorials shipped with this package – covering matrix
foundations, the Singular Value Decomposition, Principal Components
Analysis, and ancestry discrimination using the package's example
genotype data. With `name` supplied, that tutorial is launched directly.
With `name` missing in an interactive session, a menu is shown with a
checkbox beside each tutorial already marked done; the menu keeps
returning after each tutorial finishes (offering to mark it done) until
you choose 0 to quit. Progress is kept in a small file in your own R
user data directory (see
[`R_user_dir`](https://rdrr.io/r/tools/userdir.html)), so it survives
restarting R; it is self-reported and is not sent anywhere. Each
tutorial has its own "Exit tutorial" button that cleanly stops the
tutorial's web server and returns control to the R console, so there is
no need to interrupt with Ctrl-C. Requires the `learnr` package.

## Examples

``` r
if (interactive()) PlinkMatrix_tutorials()
```
