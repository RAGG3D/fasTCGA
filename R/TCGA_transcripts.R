# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
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
#' @param .foldername is the folder name to save and resource TCGA data
#' @param .project_name is the name of a TCGA project
#' @return a data frame...
#'
#' @import dplyr
#' @import tidyverse
#' @import tidybulk
#' @import tibble
#' @import purrr
#' @import TCGAbiolinks
#'
#' @export

#Download the STAR read count files from TCGA database
#If you get errors do:
#remotes::install_github("https://github.com/BioinformaticsFMRP/TCGAbiolinks",ref = "devel")
TCGA_STAR_download <- function(.foldername, .project_name){
  dir.create(.foldername)
  query <- TCGAbiolinks::GDCquery(project = .project_name,
                    data.category =  "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
  TCGAbiolinks::GDCdownload(query, .foldername)
  TCGAbiolinks::getResults(query) %>%
    write.csv(file = paste0(.foldername,"GDCdata/", .project_name, "/sample_sheet.csv"), col.names = F)
  GDCquery_clinic(project = .project_name, type = "clinical") %>%
    write.csv(file = paste0(.foldername,"GDCdata/", .project_name, "/clinical.csv"), col.names = F)

  }



#Aggregate RNA-seq read count files into a single data frame and do TMM normalization
TCGA_transcript <- function(.foldername, .project_name){
  .cancerfolder = paste0(.foldername,"GDCdata/", .project_name, "/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/")
  tibble::as_tibble(
    data.table::setDT(
      purrr::map_dfr(as.list(list.files(.cancerfolder)),
                     function(i){
                       j = list.files(paste0(.cancerfolder, "/", i))
                       data.table::fread(paste0(.cancerfolder,"/", i, "/", j)) %>%
                         na.omit() %>%
                         dplyr::mutate(id = i)
                     }) %>%
        dplyr::inner_join(read_csv(paste0(.foldername,"GDCdata/", .project_name, "/sample_sheet.csv"))) %>%
        tidybulk::rename(raw_count = unstranded,
                         symbol = gene_name,
                         sample = cases.submitter_id) %>%
        dplyr::select(sample, symbol, raw_count)
    )[, list(raw_count = sum(raw_count)), by = c("sample", "symbol")]) %>%  #Sum the raw_counts of duplicated rows
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol)
}
