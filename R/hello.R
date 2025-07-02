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

TCGA_transcript <- function(cancer){
  as_tibble(
    data.table::setDT(
      map_dfr(as.list(list.files(cancer)),
              function(i){
                data.table::fread(paste0(cancer, i)) %>%
                  mutate(sample = i)
              }) %>%
        tidybulk::rename(raw_count = `V2`) %>%
        mutate(ensembl_id = gsub("\\..*", "", V1)) %>%
        inner_join(toTable(org.Hs.egENSEMBL)) %>%
        inner_join(toTable(org.Hs.egSYMBOL)) %>%
        dplyr::select(sample, symbol, raw_count) %>%
        mutate(sample = gsub("counts.*", "counts.gz", sample)) %>%
        inner_join(
          read_csv(paste0(cancer, "/gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
        mutate(sample = `Case ID`) %>%
        dplyr::select(sample, symbol, raw_count)
    )[, list(raw_count = sum(raw_count)), by = c("sample", "symbol")]) %>%  #Sum the raw_counts of duplicated rows
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol)
}
