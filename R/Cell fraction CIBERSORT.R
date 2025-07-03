# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#' My function called sets
#'
#' @param x is a data frame contains at least 3 columns (sample, symbol and raw read count)
#' @param c1 is the name of the column contains sample names
#' @param c2 is the name of the column contains gene symbols
#' @param c3 is the name of the column contains raw read counts
#' @param ref is the data frame of signature matrix, see data/my_ref.rds
#' @param meth is the method, default "cibersort"
#' @paramact is the action var, default "get"
#' @return a data frame...
#'
#' @import dplyr
#' @import tidybulk
#' @import tibble
#' @import purrr
#' @import TCGAbiolinks
#'
#' @export

CIBERSORT <- function(x, c1, c2, c3, ref,
                      meth = "cibersort",
                      act = "get") {
  tibble::as_tibble(setDT(x)[, lapply(.SD, sum), by = c(c1, c2), .SDcols = c3]) %>%
    tidybulk::deconvolve_cellularity(
      .sample = !!sym(c1),
      .transcript = !!sym(c2),
      .abundance = !!sym(c3),
      reference = ref,
      method = meth,
      action = act
    )
}
