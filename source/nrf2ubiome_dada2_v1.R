# |----------------------------------------------------------------------------------|
# | Project: Nrf2 BL6 PEITC 16S microbiome data analysis                             |
# | Script: Data processing                                                          |
# | Coordinator: Rasika Hudlikar, Ran Yin                                            |
# | Author: Davit Sargsyan                                                           |
# | Created: 11/28/2018                                                              |
# | Modified: 12/01/2018 (DS): installed and ran DADA2 v1.10.0                       |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_nrf2ubiome_dada2_v1.txt")
date()

# Increase mmemory size to 64 Gb----
invisible(utils::memory.limit(65536))

# On Windows set multithread=FALSE----
mt <- TRUE

# # Source: https://benjjneb.github.io/dada2/index.html
# # DADA2 installation (12/01/2018)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.8")
# BiocManager::install("phyloseq", version = "3.8")

# # To install on server:
# # 1. Install External Dependencies
# source("http://bioconductor.org/biocLite.R")
# biocLite(suppressUpdates = FALSE)
# biocLite("ShortRead", suppressUpdates = FALSE)
# # You might need to install some more packages manually
# # 2. Install using devtools
# biocLite("devtools")
# library("devtools")
# devtools::install_github("benjjneb/dada2")

# Follow the tutorial;
# https://benjjneb.github.io/dada2/tutorial.html

require(data.table)
require(dada2)
require(phyloseq)
require(ggplot2)

# Get FastQ file names----
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

samples$Week <- factor(substr(samples$Name, 1, 1))
samples$Cage <- factor(substr(samples$Name, 2, 2))
samples$Sex <- "NA"
samples$Sex[samples$Cage %in% c("A", "B")] <- "Male"
samples$Sex[samples$Cage =="C"] <- "Female"
samples$Sex <- factor(samples$Sex)
samples$Diet <- factor(substr(samples$Name, 5, 5),
                       levels = c("C", 
                                  "P"),
                       labels = c("AIN93M Control",
                                  "PEITC"))
samples$MouseNum <- factor(as.numeric(substr(samples$Name, 3, 3)))
samples$ID <- paste(samples$Diet,
                    samples$Cage,
                    samples$MouseNum,
                    sep = "_")
samples$Diet_Week <- paste(samples$Diet,
                           samples$Week,
                           sep = "_")

samples <- data.frame(samples)
samples

save(samples,
     file = "data/samples.RData")

# Forward and reverse fastq filenames----
# Must have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

# Reads' quality----
p1 <- plotQualityProfile(fnFs[1:16])
tiff(filename = "tmp/qc16_forward.tiff",
     height = 7,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

p2 <- plotQualityProfile(fnRs[1:16])
tiff(filename = "tmp/qc16_reverse.tiff",
     height = 7,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Filter and trim----
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

p3 <- plotQualityProfile(fnFs[1]) +
  scale_x_continuous(breaks = seq(0, 300, 10)) +
  geom_vline(xintercept = c(6, 285),
             linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
tiff(filename = "tmp/trim_forward.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p3)
graphics.off()
  
p4 <- plotQualityProfile(fnRs[1]) +
  scale_x_continuous(breaks = seq(0, 300, 10)) +
  geom_vline(xintercept = c(6, 195),
             linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
tiff(filename = "tmp/trim_reverse.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p4)
graphics.off()

out <- filterAndTrim(fwd = fnFs, 
                     filt = filtFs, 
                     rev = fnRs, 
                     filt.rev = filtRs,
                     trimLeft = c(6, 6),
                     truncLen = c(280, 190),
                     rm.phix = TRUE,
                     compress = TRUE,
                     verbose = TRUE,
                     multithread = mt) 
head(out)
# 4A_S1_L001_R1_001.fastq.gz       602046    454068
# 4B_S2_L001_R1_001.fastq.gz       353067    274418
# 4C_S3_L001_R1_001.fastq.gz       638909    527488
# 5A1-C_S4_L001_R1_001.fastq.gz    377185    308643
# 5A1-P_S22_L001_R1_001.fastq.gz   906226    749294
# 5A2-C_S7_L001_R1_001.fastq.gz    430657    338116

# Learn the Error Rates
errF <- learnErrors(filtFs, 
                    multithread = mt)
plotErrors(errF, 
           nominalQ = TRUE)

errR <- learnErrors(filtRs, 
                    multithread = mt)
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
                       refFasta = "fastq/tax/silva_nr_v128_train_set.fa.gz",
                       multithread = TRUE)
dim(taxa)
View(head(taxa))
save(taxa,
     file = "data/taxa.RData")

# # SIDE ANALYSIS: save taxa obtained before removeing bimeras----
# taxa_with_bimera <- assignTaxonomy(seqs = seqtab,
#                                    refFasta = "fastq/tax/silva_nr_v128_train_set.fa.gz",
#                                    multithread = TRUE)
# save(taxa_with_bimera,
#      file = "data/taxa_with_bimera.RData")
# dim(taxa_with_bimera)
# View(head(taxa_with_bimera))

# load("data/taxa_with_bimera.RData")
# taxa <- taxa_with_bimera

# Add species----
load("data/taxa.RData")
dim(taxa)
View(head(taxa))

load("data/seqtab.nochim.RData")
dim(seqtab.nochim)
View(head(t(seqtab.nochim)))
length(colnames(seqtab.nochim))

# Keep only sequences found in data----
taxa.tmp <- taxa[rownames(taxa) %in% colnames(seqtab.nochim), ]
View(head(taxa.tmp))

# Add species----
taxa.plus <- addSpecies(taxtab = taxa.tmp,
                        refFasta = "fastq/tax/silva_species_assignment_v128.fa.gz",
                        verbose = TRUE)
dim(taxa.plus)
View(head(taxa.plus))
View(taxa.plus[!is.na(taxa.plus[, 7]), ])
write.csv(taxa.plus[!is.na(taxa.plus[, 7]), ],
          file = "tmp/species.csv")

table(!is.na(taxa.plus[, 7]))
# FALSE  TRUE 
# 10754     5 
save(taxa.plus,
     file = "data/taxa.plus.RData")
gc()

# Removing sequence rownames for display only----
taxa.print <- taxa.plus 
rownames(taxa.print) <- NULL
head(taxa.print)

# Handoff to phyloseq----
load("data/samples.RData")
samples
load("data/seqtab.nochim.RData")
load("data/taxa.plus.RData")

# Keep only mapped species data----
taxa.spc <- taxa.plus[!is.na(taxa.plus)]
ps <- phyloseq(otu_table(seqtab.nochim, 
                         taxa_are_rows = FALSE), 
               sample_data(samples), 
               tax_table(taxa.plus))
ps
save(ps,
     file = "data/ps.RData")

# Keep Week5 samples only
# ps <- prune_samples(sample_names(ps) %in% rownames(samples)[samples$Week %in% c(5, 9)], ps)
ps <- prune_samples(sample_names(ps) %in% rownames(samples)[samples$Week != "U"], ps)
ps

p1 <- plot_richness(ps,
                    x = "Diet_Week", 
                    measures = "Shannon",
                    color = "Week") +
  geom_line(aes(group = ID),
            color = "black") +
  geom_point(shape = 21,
             size = 3,
             color = "black")
p1

tiff(filename = "tmp/richness.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps,
                                   function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop,
                          method = "NMDS",
                          distance = "bray")

p1 <- plot_ordination(ps.prop,
                      ord.nmds.bray,
                      color = "Diet_Week",
                      title = "Bray NMDS") +
  geom_point(size = 3)
p1
tiff(filename = "tmp/nmds.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Barplots
View(head(taxa_sums(ps)))
names(taxa_sums(ps))[1:10]
top20 <- names(sort(taxa_sums(ps),
                    decreasing = TRUE))[1:10]
ps.top20 <- transform_sample_counts(ps,
                                    function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20,
                       ps.top20)
ps.top20@tax_table@.Data

plot_bar(ps.top20, 
         x = "Diet", 
         fill = "Order") + 
  facet_wrap( ~ Week, 
              scales = "free_x")

plot_bar(ps.top20, 
         x = "Diet", 
         fill = "Family") + 
  facet_wrap( ~ Week, 
              scales = "free_x") +
  geom_tile(color = "white")

plot_bar(ps.top20, 
         x = "Diet", 
         fill = "Genus") + 
  facet_wrap( ~ Week, 
              scales = "free_x")

ps.species <- ps
top20 <- names(sort(taxa_sums(ps),
                    decreasing = TRUE))[1:10]
ps.top20 <- transform_sample_counts(ps,
                                    function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20,
                       ps.top20)
ps.top20@tax_table@.Data
plot_bar(ps.top20, 
         x = "Diet", 
         fill = "Species") + 
  facet_wrap( ~ Week, 
              scales = "free_x")

# OTU table----
tn <- taxa_names(ps)
length(tn)
tn[1:10]

?otu_table
dotu <- otu_table(object = ps)
dim(dotu)
rownames(dotu) <- NULL
head(dotu)

t1 <- taxa_sums(ps)
View(t1)
#sessionInfo()
# sink()