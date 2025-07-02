\name{fasTCGA}
\alias{fasTCGA}
\title{Tailored for TCGA data alignment and formation}
\usage{
TCGA_STAR_download(.folder, .project_name)
TCGA_transcript(.foldername, .project_name)
}
\description{
Aggregate individual patient read files into one unified data frame
}
\examples{
.project_name = "TCGA-ACC"
.foldername = "C:/Rstudio/test/"

TCGA_STAR_download(.folder, .project_name)
TCGA_transcript(.foldername, .project_name)
}
