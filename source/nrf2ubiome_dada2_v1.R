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

# Increase mmemory size to 64 Gb----
invisible(utils::memory.limit(65536))

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
require(phyloseq)
require(ggplot2)

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
plotQualityProfile(fnFs[1:16])
plotQualityProfile(fnRs[1:16])

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

# Dereplication----
# NOTE: the derepFastq function does not perform garbage collection 
# so we will run each file separately and then merge them all into a list
# Forward----
derepFs <- list()
for (i in 1:length(filtFs)) {
  derepFs[[i]] <- derepFastq(filtFs[i],
                        verbose = TRUE)
  gc()
}
# Dereplicating sequence entries in Fastq file: fastq/filtered/4A_F_filt.fastq.gz
# Encountered 133082 unique sequences from 454068 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/4B_F_filt.fastq.gz
# Encountered 85322 unique sequences from 274418 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/4C_F_filt.fastq.gz
# Encountered 154819 unique sequences from 527488 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A1-C_F_filt.fastq.gz
# Encountered 113156 unique sequences from 308643 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A1-P_F_filt.fastq.gz
# Encountered 200364 unique sequences from 749294 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A2-C_F_filt.fastq.gz
# Encountered 113447 unique sequences from 338116 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A2-P_F_filt.fastq.gz
# Encountered 139753 unique sequences from 468031 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A3-C_F_filt.fastq.gz
# Encountered 245533 unique sequences from 996829 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A3-P_F_filt.fastq.gz
# Encountered 141311 unique sequences from 491831 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B1-C_F_filt.fastq.gz
# Encountered 156135 unique sequences from 505273 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B1-P_F_filt.fastq.gz
# Encountered 126238 unique sequences from 386585 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B2-C_F_filt.fastq.gz
# Encountered 100726 unique sequences from 281824 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B2-P_F_filt.fastq.gz
# Encountered 140012 unique sequences from 440874 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B3-C_F_filt.fastq.gz
# .Encountered 258331 unique sequences from 1018107 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B3-P_F_filt.fastq.gz
# Encountered 178440 unique sequences from 640215 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C1-C_F_filt.fastq.gz
# Encountered 157169 unique sequences from 463064 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C1-P_F_filt.fastq.gz
# Encountered 164074 unique sequences from 478623 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C2-C_F_filt.fastq.gz
# Encountered 229687 unique sequences from 714263 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C2-P_F_filt.fastq.gz
# Encountered 136507 unique sequences from 428389 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C3-C_F_filt.fastq.gz
# Encountered 204590 unique sequences from 642892 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C3-P_F_filt.fastq.gz
# Encountered 147402 unique sequences from 524426 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A1-C_F_filt.fastq.gz
# Encountered 176295 unique sequences from 677651 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A1-P_F_filt.fastq.gz
# Encountered 168997 unique sequences from 579665 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A2-C_F_filt.fastq.gz
# Encountered 196521 unique sequences from 564245 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A2-P_F_filt.fastq.gz
# Encountered 162393 unique sequences from 481984 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A3-C_F_filt.fastq.gz
# Encountered 161049 unique sequences from 648809 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A3-P_F_filt.fastq.gz
# Encountered 175602 unique sequences from 577087 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B1-C_F_filt.fastq.gz
# Encountered 253662 unique sequences from 692453 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B1-P_F_filt.fastq.gz
# Encountered 206536 unique sequences from 762689 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B2-C_F_filt.fastq.gz
# Encountered 157047 unique sequences from 428241 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B2-P_F_filt.fastq.gz
# Encountered 151187 unique sequences from 580128 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B3-C_F_filt.fastq.gz
# Encountered 118049 unique sequences from 395495 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B3-P_F_filt.fastq.gz
# Encountered 122109 unique sequences from 436206 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C1-C_F_filt.fastq.gz
# Encountered 237552 unique sequences from 639587 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C1-P_F_filt.fastq.gz
# Encountered 95174 unique sequences from 305144 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C2-C_F_filt.fastq.gz
# Encountered 138495 unique sequences from 365319 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C2-P_F_filt.fastq.gz
# Encountered 152687 unique sequences from 440962 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C3-C_F_filt.fastq.gz
# Encountered 151480 unique sequences from 469395 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C3-P_F_filt.fastq.gz
# Encountered 134245 unique sequences from 409244 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/Undetermined_F_filt.fastq.gz
# Encountered 540890 unique sequences from 987927 total sequences read.

derepFs <- do.call("c", derepFs)
# 40 elements, 16.2 Gb of data in memory
names(derepFs) <- sample.names
gc()
save(derepFs,
     file = "data/derepFs.RData")

# Reverse----
derepRs <- list()
for (i in 1:length(filtRs)) {
  derepRs[[i]] <- derepFastq(filtRs[i],
                             verbose = TRUE)
  gc()
}

# Dereplicating sequence entries in Fastq file: fastq/filtered/4A_R_filt.fastq.gz
# Encountered 187259 unique sequences from 454068 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/4B_R_filt.fastq.gz
# Encountered 86111 unique sequences from 274418 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/4C_R_filt.fastq.gz
# Encountered 190092 unique sequences from 527488 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A1-C_R_filt.fastq.gz
# Encountered 117403 unique sequences from 308643 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A1-P_R_filt.fastq.gz
# Encountered 238705 unique sequences from 749294 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A2-C_R_filt.fastq.gz
# Encountered 107327 unique sequences from 338116 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A2-P_R_filt.fastq.gz
# Encountered 145361 unique sequences from 468031 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A3-C_R_filt.fastq.gz
# Encountered 261260 unique sequences from 996829 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5A3-P_R_filt.fastq.gz
# Encountered 149806 unique sequences from 491831 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B1-C_R_filt.fastq.gz
# Encountered 164805 unique sequences from 505273 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B1-P_R_filt.fastq.gz
# Encountered 137163 unique sequences from 386585 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B2-C_R_filt.fastq.gz
# Encountered 105888 unique sequences from 281824 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B2-P_R_filt.fastq.gz
# Encountered 131631 unique sequences from 440874 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B3-C_R_filt.fastq.gz
# .Encountered 287131 unique sequences from 1018107 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5B3-P_R_filt.fastq.gz
# Encountered 192519 unique sequences from 640215 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C1-C_R_filt.fastq.gz
# Encountered 162026 unique sequences from 463064 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C1-P_R_filt.fastq.gz
# Encountered 173240 unique sequences from 478623 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C2-C_R_filt.fastq.gz
# Encountered 231903 unique sequences from 714263 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C2-P_R_filt.fastq.gz
# Encountered 140437 unique sequences from 428389 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C3-C_R_filt.fastq.gz
# Encountered 222656 unique sequences from 642892 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/5C3-P_R_filt.fastq.gz
# Encountered 153555 unique sequences from 524426 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A1-C_R_filt.fastq.gz
# Encountered 183898 unique sequences from 677651 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A1-P_R_filt.fastq.gz
# Encountered 174877 unique sequences from 579665 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A2-C_R_filt.fastq.gz
# Encountered 211843 unique sequences from 564245 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A2-P_R_filt.fastq.gz
# Encountered 151309 unique sequences from 481984 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A3-C_R_filt.fastq.gz
# Encountered 165212 unique sequences from 648809 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9A3-P_R_filt.fastq.gz
# Encountered 184669 unique sequences from 577087 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B1-C_R_filt.fastq.gz
# Encountered 223020 unique sequences from 692453 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B1-P_R_filt.fastq.gz
# Encountered 264229 unique sequences from 762689 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B2-C_R_filt.fastq.gz
# Encountered 168556 unique sequences from 428241 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B2-P_R_filt.fastq.gz
# Encountered 161363 unique sequences from 580128 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B3-C_R_filt.fastq.gz
# Encountered 122985 unique sequences from 395495 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9B3-P_R_filt.fastq.gz
# Encountered 130075 unique sequences from 436206 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C1-C_R_filt.fastq.gz
# Encountered 244020 unique sequences from 639587 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C1-P_R_filt.fastq.gz
# Encountered 105653 unique sequences from 305144 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C2-C_R_filt.fastq.gz
# Encountered 133821 unique sequences from 365319 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C2-P_R_filt.fastq.gz
# Encountered 140297 unique sequences from 440962 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C3-C_R_filt.fastq.gz
# Encountered 144144 unique sequences from 469395 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/9C3-P_R_filt.fastq.gz
# Encountered 193089 unique sequences from 409244 total sequences read.
# Dereplicating sequence entries in Fastq file: fastq/filtered/Undetermined_R_filt.fastq.gz
# Encountered 564523 unique sequences from 987927 total sequences read.

derepRs <- do.call("c", derepRs)
# 40 elements, 11.8 Gb of data in memory
names(derepRs) <- sample.names
gc()
save(derepRs,
     file = "data/derepRs.RData")

# rm(derepFs)
# rm(derepRs)
# gc()
# load("data/derepFs.RData")
# load("data/derepRs.RData")

# Sample Inference
dadaFs <- dada(derepFs, 
               err = errF,
               multithread = FALSE)
save(dadaFs,
     file = "data/dadaFs.RData")
head(dadaFs)

dadaRs <- dada(derepRs, 
               err = errF,
               multithread = FALSE)
save(dadaRs,
     file = "data/dadaRs.RData")

# Merge paired reads----
# load("data/dadaFs.RData")
# load("data/dadaRs.RData")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
save(mergers,
     file = "data/mergers.RData")

# load("data/mergers.RData")
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras----
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method = "consensus",
                                    multithread = FALSE, 
                                    verbose = TRUE)
dim(seqtab.nochim)
save(seqtab.nochim,
     file = "data/seqtab.nochim.RData")

# Assign taxonomy----
taxa <- assignTaxonomy(seqs = seqtab.nochim,
                       refFasta = "~/tax/silva_nr_v128_train_set.fa.gz",
                       multithread = TRUE)
dim(taxa)
View(head(taxa))
save(taxa,
     file = "data/taxa.RData")

# # SIDE ANALYSIS: save taxa obtained before removeing bimeras----
# taxa_with_bimera <- taxa
# save(taxa_with_bimera,
#      file = "data/taxa_with_bimera.RData")

# Removing sequence rownames for display only----
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

# Handoff to phyloseq----
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps

plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

# Barplots
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")

#sessionInfo()
# sink()