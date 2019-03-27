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

# Data----
require(data.table)
require(phyloseq)
require(ggplot2)
require(plotly)
require(DT)
require(shiny)
source("source/functions_v2.R")

# Load data----
# Counts
load("data/ps.RData")

# Taxonomy
load("data/taxa.plus.RData")
taxa <- data.table(seq16s = rownames(taxa.plus),
                   taxa.plus)

# Samples
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
DT::datatable(samples)

# Part I: keep controls and PEITCs at weeks 5 and 9 only
w59 <- prune_samples(samples = sample_names(ps)[!(sample_names(ps) %in%
                                                    c("4A",
                                                      "4B",
                                                      "4C",
                                                      "Undetermined"))],
                    x = ps)

p1 <- plot_richness(w59,
                    x = "Diet_Week", 
                    measures = "Shannon",
                    color = "Week") +
  geom_point(shape = 21,
             size = 3,
             color = "black") +
  geom_line(aes(group = ID),
            color = "black")
p1

# Save the plot as a TIFF file
tiff(filename = "tmp/nrf2ubiome_shannon.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# OTU table
otu <- t(w59@otu_table@.Data)
otu <- data.table(seq16s = rownames(otu),
                  otu)

# Merge taxonomy and counts tables----
dt1 <- merge(taxa,
             otu,
             by = "seq16s")
dt1$seq16s <- NULL

# Remove archea and eucaryota----
dt1 <- droplevels(dt1[Kingdom == "Bacteria", ])
t1 <- colSums(dt1[, 8:ncol(dt1)])
t1 <- data.table(Sample = names(t1),
                 Total = t1)
t1 <- merge(unique(samples[, c("Sample",
                               "Diet",
                               "Week")]),
            t1,
            by = "Sample")

p1 <- ggplot(t1,
             aes(x = Sample,
                 y = Total,
                 fill = Diet,
                 colour = Week)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggplotly(p1)

# Class----
tax <- "Class"
dt2 <- counts_ra_by_tax_rank(dt1, tax)

# Plot samples
p1 <- ggplot(dt2[order(RA,
                       decreasing = TRUE), ],
             aes(x = Sample,
                 y = RA,
                 fill = Tax,
                 group = Diet)) +
  facet_wrap(~ Week + Diet,
             scales = "free_x",
             nrow = 1) +
  geom_bar(stat = "identity") +
  scale_fill_discrete(tax) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggplotly(p1)

# Plot means
p2 <- ggplot_mean_ra(dt2, tax)
# Save the plot as a TIFF file
tiff(filename = "tmp/nrf2ubiome_class.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

mu.class <- data.table(mean_ra(dt2, tax, facet_sex = FALSE))
mu.class <- dcast.data.table(mu.class,
                             Tax ~ Diet + Week,
                             fill = "x")
mu.class <- mu.class[order(mu.class$`AIN93M Control_5`,
                           decreasing = TRUE)]
mu.class$mu <- 100*(rowSums(mu.class[, -1])/4)
mu.class
DT::datatable(mu.class)

#sessionInfo()
# sink()