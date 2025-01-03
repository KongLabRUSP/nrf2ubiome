# |----------------------------------------------------------------------------------|
# | Project: Nrf2 BL6 PEITC 16S microbiome data analysis                             |
# | Script: Data visualization                                                       |
# | Coordinator: Rasika Hudlikar, Ran Yin                                            |
# | Author: Davit Sargsyan                                                           |
# | Created: 02/07/2019                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_nrf2ubiome_data_simulation_v1.txt")
date()

# Increase mmemory size to 64 Gb----
invisible(utils::memory.limit(65536))

require(data.table)
require(phyloseq)
require(ggplot2)
require(plotly)

# Load data----
# 1. Counts
load("data/ps.RData")

# 2. Taxonomy
load("data/taxa.plus.RData")
taxa <- data.table(seq16s = rownames(taxa.plus),
                   taxa.plus)

# 3. Samples
# ps@sam_data
load("data/samples.RData")
samples$Sample <- substr(x = samples$Name,
                         start = 1,
                         stop = 5)
samples$Sample[samples$Sample %in% c("4A_S1",
                                     "4B_S2",
                                     "4C_S3")] <- 
  substr(x = samples$Sample[samples$Sample %in% c("4A_S1",
                                                  "4B_S2",
                                                  "4C_S3")],
         start = 1,
         stop = 2)
samples
samples$Diet <- as.character(samples$Diet)
samples$Diet[is.na(samples$Diet)] <- "None"

# Part I: keep controls and PEITCs at weeks 5 and 9 only
# w59 <- prune_samples(samples = sample_names(ps)[!(sample_names(ps) %in%
#                                                     c("4A",
#                                                       "4B",
#                                                       "4C",
#                                                       "Undetermined"))], 
#                     x = ps)
w59 <- prune_samples(samples = sample_names(ps)[sample_names(ps) !="Undetermined"], 
                     x = ps)
w59
plot_richness(w59,
              x = "Diet_Week", 
              measures = "Shannon",
              color = "Week") +
  geom_line(aes(group = ID),
            color = "black") +
  geom_point(shape = 21,
             size = 3,
             color = "black")

# OTU table
otu <- t(w59@otu_table@.Data)
otu <- data.table(seq16s = rownames(otu),
                  otu)

# Merge taxonomy and counts tables----
dt1 <- merge(taxa,
             otu,
             by = "seq16s")
dt1$seq16s <- NULL
dt1
table(dt1$Kingdom)
# Archaea  Bacteria Eukaryota 
# 4        10197    472 

# Remove archea and eucaryota----
dt1 <- droplevels(dt1[Kingdom == "Bacteria", ])

# TAXONOMY LEVELS----
# **K**ing **P**hillip **C**an n**O**t **F**ind **G**een **S**ocks
# Kingdom           
# Phylum               
# Class              
# Order              
# Family
# Genus
# Species
# NOTE: work only on Class to Genus levels!

counts <- dt1[, 
              8:ncol(dt1),
              with = FALSE]
counts

slegend <- ps@sam_data

# 1. Aggregate counts by Class----
dt.c <- lapply(counts,
               function(a) {
                 out <- aggregate(x = a,
                                  by = list(Class = dt1$Class),
                                  FUN = "sum")
               })
dt.c

dt.c <- Reduce(f = function(...){merge(..., 
                                       by = "Class",
                                       all = TRUE)},
               x = dt.c)
colnames(dt.c)[-1] <- colnames(counts)
dt.c

# Relative abundance----
dtr.c <- data.table(apply(dt.c[, -1],
                          2,
                          function(a) {
                            a/sum(a)
                          }))
dtr.c$Class <- dt.c$Class
dtr.c

dtr.c <- melt.data.table(data = dtr.c,
                         id.vars = ncol(dtr.c),
                         measure.vars = 1:(ncol(dtr.c) - 1),
                         variable.name = "Sample",
                         value.name = "RA")
dtr.c

# Merge counts and sample info----
dtr.c <- merge(samples,
               dtr.c,
               by = "Sample")
#dtr.c$Week <- factor(dtr.c$Week)
dtr.c$Sample <- factor(dtr.c$Sample,
                       levels = unique(dtr.c$Sample))
dtr.c

p1 <- ggplot(dtr.c[order(RA,
                         decreasing = TRUE), ],
       aes(x = Sample,
           y = RA,
           fill = Class,
           group = Diet)) +
  facet_wrap(~ Week + Diet,
             scales = "free_x",
             nrow = 1) +
  geom_bar(stat = "identity") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggplotly(p1)

# 2. Aggregate counts by Order----
dt.o <- lapply(counts,
               function(a) {
                 out <- aggregate(x = a,
                                  by = list(Class = dt1$Class,
                                            Order = dt1$Order),
                                  FUN = "sum")
               })

dt.o <- Reduce(f = function(...){merge(..., 
                                       by = c("Class",
                                              "Order"),
                                       all = TRUE)},
               x = dt.o)
colnames(dt.o)[-c(1:2)] <- colnames(counts)

# Relative abundance----
dtr.o <- data.table(apply(dt.o[, -c(1:2)],
                          2,
                          function(a) {
                            a/sum(a)
                          }))
dtr.o$Class <- dt.o$Class
dtr.o$Order <- dt.o$Order

dtr.o <- melt.data.table(data = dtr.o,
                         id.vars = (ncol(dtr.o) - 1):ncol(dtr.o),
                         measure.vars = 1:(ncol(dtr.o) - 2),
                         variable.name = "Sample",
                         value.name = "RA")

# Merge counts and sample info----
dtr.o <- merge(samples,
               dtr.o,
               by = "Sample")
dtr.o$Sample <- factor(dtr.o$Sample,
                       levels = unique(dtr.o$Sample))

p1 <- ggplot(dtr.o[order(RA,
                         decreasing = TRUE), ],
             aes(x = Sample,
                 y = RA,
                 color = Class,
                 fill = Order,
                 group = Diet)) +
  facet_wrap(~ Week + Diet,
             scales = "free_x",
             nrow = 1) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggplotly(p1)

# 3. Aggregate counts by Family----
dt.f <- lapply(counts,
               function(a) {
                 out <- aggregate(x = a,
                                  by = list(Class = dt1$Class,
                                            Family = dt1$Family),
                                  FUN = "sum")
               })

dt.f <- Reduce(f = function(...){merge(..., 
                                       by = c("Class",
                                              "Family"),
                                       all = TRUE)},
               x = dt.f)
colnames(dt.f)[-c(1:2)] <- colnames(counts)

# Relative abundance----
dtr.f <- data.table(apply(dt.f[, -c(1:2)],
                          2,
                          function(a) {
                            a/sum(a)
                          }))
dtr.f$Class <- dt.f$Class
dtr.f$Family <- dt.f$Family

dtr.f <- melt.data.table(data = dtr.f,
                         id.vars = (ncol(dtr.f) - 1):ncol(dtr.f),
                         measure.vars = 1:(ncol(dtr.f) - 2),
                         variable.name = "Sample",
                         value.name = "RA")

# Merge counts and sample info----
dtr.f <- merge(samples,
               dtr.f,
               by = "Sample")
dtr.f$Sample <- factor(dtr.f$Sample,
                       levels = unique(dtr.f$Sample))

p1 <- ggplot(dtr.f[order(RA,
                         decreasing = TRUE), ],
             aes(x = Sample,
                 y = RA,
                 color = Class,
                 fill = Family,
                 group = Diet)) +
  facet_wrap(~ Week + Diet,
             scales = "free_x",
             nrow = 1) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggplotly(p1)

# 4. Aggregate counts by Genus----
dt.g <- lapply(counts,
               function(a) {
                 out <- aggregate(x = a,
                                  by = list(Class = dt1$Class,
                                            Genus = dt1$Genus),
                                  FUN = "sum")
               })

dt.g <- Reduce(f = function(...){merge(..., 
                                       by = c("Class",
                                              "Genus"),
                                       all = TRUE)},
               x = dt.g)
colnames(dt.g)[-c(1:2)] <- colnames(counts)

# Relative abundance----
dtr.g <- data.table(apply(dt.g[, -c(1:2)],
                          2,
                          function(a) {
                            a/sum(a)
                          }))
dtr.g$Class <- dt.g$Class
dtr.g$Genus <- dt.g$Genus

dtr.g <- melt.data.table(data = dtr.g,
                         id.vars = (ncol(dtr.g) - 1):ncol(dtr.g),
                         measure.vars = 1:(ncol(dtr.g) - 2),
                         variable.name = "Sample",
                         value.name = "RA")

# Merge counts and sample info----
dtr.g <- merge(samples,
               dtr.g,
               by = "Sample")
dtr.g$Sample <- factor(dtr.g$Sample,
                       levels = unique(dtr.g$Sample))

p1 <- ggplot(dtr.g[order(RA,
                         decreasing = TRUE), ],
             aes(x = Sample,
                 y = RA,
                 color = Class,
                 fill = Genus,
                 group = Diet)) +
  facet_wrap(~ Week + Diet,
             scales = "free_x",
             nrow = 1) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggplotly(p1)

# # Sample legend
# ps@sam_data
# 
# p1 <- plot_richness(ps,
#                     x = "Diet_Week", 
#                     measures = "Shannon",
#                     color = "Week") +
#   geom_line(aes(group = ID),
#             color = "black") +
#   geom_point(shape = 21,
#              size = 3,
#              color = "black")
# p1
# 
# tiff(filename = "tmp/richness.tiff",
#      height = 5,
#      width = 6,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# print(p1)
# graphics.off()
# 
# # Transform data to proportions as appropriate for Bray-Curtis distances
# ps.prop <- transform_sample_counts(ps,
#                                    function(otu) otu/sum(otu))
# ord.nmds.bray <- ordinate(ps.prop,
#                           method = "NMDS",
#                           distance = "bray")
# 
# p1 <- plot_ordination(ps.prop,
#                       ord.nmds.bray,
#                       color = "Diet_Week",
#                       title = "Bray NMDS") +
#   geom_point(size = 3)
# p1
# tiff(filename = "tmp/nmds.tiff",
#      height = 5,
#      width = 7,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# print(p1)
# graphics.off()
# 
# # Barplots
# View(head(taxa_sums(ps)))
# names(taxa_sums(ps))[1:10]
# top20 <- names(sort(taxa_sums(ps),
#                     decreasing = TRUE))[1:10]
# ps.top20 <- transform_sample_counts(ps,
#                                     function(OTU) OTU/sum(OTU))
# ps.top20 <- prune_taxa(top20,
#                        ps.top20)
# ps.top20@tax_table@.Data
# 
# plot_bar(ps.top20, 
#          x = "Diet", 
#          fill = "Order") + 
#   facet_wrap( ~ Week, 
#               scales = "free_x")
# 
# plot_bar(ps.top20, 
#          x = "Diet", 
#          fill = "Family") + 
#   facet_wrap( ~ Week, 
#               scales = "free_x") +
#   geom_tile(color = "white")
# 
# plot_bar(ps.top20, 
#          x = "Diet", 
#          fill = "Genus") + 
#   facet_wrap( ~ Week, 
#               scales = "free_x")
# 
# ps.species <- ps
# top20 <- names(sort(taxa_sums(ps),
#                     decreasing = TRUE))[1:10]
# ps.top20 <- transform_sample_counts(ps,
#                                     function(OTU) OTU/sum(OTU))
# ps.top20 <- prune_taxa(top20,
#                        ps.top20)
# ps.top20@tax_table@.Data
# plot_bar(ps.top20, 
#          x = "Diet", 
#          fill = "Species") + 
#   facet_wrap( ~ Week, 
#               scales = "free_x")
# 
# # OTU table----
# tn <- taxa_names(ps)
# length(tn)
# tn[1:10]
# 
# ?otu_table
# dotu <- otu_table(object = ps, taxa_are_rows = )
# dim(dotu)
# rownames(dotu) <- NULL
# head(dotu)
# 
# t1 <- taxa_sums(ps)
# View(t1)
#sessionInfo()
# sink()