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
#' @param .data is a data frame contains a column named "sample" with TCGA barcodes
#' @param .foldername is the folder name to save and resource TCGA data
#' @param .project_name is the name of a TCGA project
#' @return a data frame...
#'
#' @import dplyr
#'
#' @export

clinical_combine <- function(.data, .foldername, .project_name) {
  .data |>
    dplyr::inner_join(read.csv(paste0(.foldername,"GDCdata/", .project_name, "/clinical.csv")), by = c("sample" = "submitter_id"))|>
    dplyr::mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) |>
    dplyr::mutate(na = is.na(total_living_days)) |>
    dplyr::mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(days_to_last_follow_up)), total_living_days)) |>
    dplyr::mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))}


