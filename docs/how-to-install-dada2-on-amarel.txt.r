# how to install dada2 on Rstudio on ondemand
# R Wu Dec 2018

# For R version 3.5 and above
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")

# For R version below 3.5 (3.4.4 on Rstudio server as of 12/11/2018)
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite("dada2")
