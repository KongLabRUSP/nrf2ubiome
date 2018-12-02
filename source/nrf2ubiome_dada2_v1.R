# |----------------------------------------------------------------------------------|
# | Project: Nrf2 BL6 PEITC 16S microbiome data analysis                             |
# | Script: Data                                                                     |
# | Coordinator: Ran Yin, Renyi Wu                                                   |
# | Author: Davit Sargsyan                                                           |
# | Created: 11/28/2018                                                              |
# | Modified: 12/01/2018 (DS): installed and ran DADA2 v1.10.0                       |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_skin_uvb_dna_v2.txt")
date()

# # Source: https://benjjneb.github.io/dada2/index.html
# # DADA2 installation (12/01/2018)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.8")

# Follow the tutorial;
# https://benjjneb.github.io/dada2/tutorial.html

require(data.table)
require(dada2)
packageVersion("dada2")
# [1] ‘1.10.0’

path <- "fastq"
list.files(path = path,
           pattern = ".gz")

# Create sample description table----
# Source: /docs/legends_16s_11-28-2018.xlsx
samples <- data.table(Name = list.files(path = path,
                                        pattern = ".gz"))
samples$Name <- gsub(pattern = "_L001_R1_001.fastq.gz",
                     replacement = "",
                     x = samples$Name)
samples$Name <- gsub(pattern = "_L001_R2_001.fastq.gz",
                     replacement = "",
                     x = samples$Name)
samples <- unique(samples)
samples

samples$Week <- factor(substr(samples$Name, 1, 1))
samples$Cage <- factor(substr(samples$Name, 2, 2))
samples$Sex <- "NA"
samples$Sex[samples$Cage %in% c("A", "B")] <- "Male"
samples$Sex[samples$Cage =="C"] <- "Female"
samples$Sex <- factor(samples$Sex)
samples$MouseNum <- factor(as.numeric(substr(samples$Name, 3, 3)))
samples$Diet <- factor(substr(samples$Name, 5, 5),
                       levels = c("C", 
                                  "P"),
                       labels = c("AIN93M Control",
                                  "PEITC"))
save(samples,
     file = "data/samples.RData")

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

# Reads' quality
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Filter and trim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fwd = fnFs, 
                     filt = filtFs, 
                     rev = fnRs, 
                     filt.rev = filtRs,
                     truncLen = c(240,160),
                     maxN = 0, 
                     maxEE = c(2,2),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = FALSE) # On Windows set multithread=FALSE
head(out)
# 4A_S1_L001_R1_001.fastq.gz       602046    454068
# 4B_S2_L001_R1_001.fastq.gz       353067    274418
# 4C_S3_L001_R1_001.fastq.gz       638909    527488
# 5A1-C_S4_L001_R1_001.fastq.gz    377185    308643
# 5A1-P_S22_L001_R1_001.fastq.gz   906226    749294
# 5A2-C_S7_L001_R1_001.fastq.gz    430657    338116

# Learn the Error Rates
errF <- learnErrors(filtFs, 
                    multithread = FALSE)
plotErrors(errF, 
           nominalQ = TRUE)

errR <- learnErrors(filtRs, 
                    multithread = FALSE)
plotErrors(errR, 
           nominalQ = TRUE)

# Sample Inference
dadaFs <- dada(derepFs, 
               err = errF,
               multithread = FALSE)

dadaRs <- dada(derepRs, 
               err = errF,
               multithread = FALSE)


#sessionInfo()
# sink()