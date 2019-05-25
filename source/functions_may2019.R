# Functions
# Created (DS): 02/09/2019
# Modified (DS): 05/24/2019

# 1. Aggregate counts by a taxonomic rank
# Input:
## dt1 = OTU table, with first 6 columns being taonomic ranks
## aggr_by = name of taxonomic level to aggregate by
## counts_start = number of column containing first sample; 
##                assuming all following columns are samples (counts)
# Output:
## Table of counts
# Example:
## counts <- counts_by_tax_rank(dt1 = out.table,
##                              arrg_by = "Phylum")
#
counts_by_tax_rank <- function(dt1, aggr_by, counts_start = 7) {
  # Counts
  counts <- dt1[, 
                counts_start:ncol(dt1),
                with = FALSE]
  
  dt.t <- lapply(counts,
                 function(a) {
                   out <- aggregate(x = a,
                                    by = list(dt1[[which(colnames(dt1) == aggr_by)]]),
                                    FUN = "sum")
                 })
  
  dt.t <- Reduce(function(x, y) {
    merge(x, 
          y, 
          all = TRUE,
          by = "Group.1")}, 
    dt.t)
  colnames(dt.t) <- c(aggr_by,
                      colnames(counts))
  dt.t <- data.table(dt.t)
  return(dt.t)
}

# 2. Compute relative abundance from counts
# Input:
## counts = output from abouve function (counts_by_tax_rank)
## pct = should the function output percents?
## digit = should the output be rounded, and if yes (non-NA), to what decimal place?
# Output:
## Table of counts
# Example:
## counts <- counts_by_tax_rank(dt1 = out.table,
##                              arrg_by = "Phylum")
#
ra_by_tax_rank <- function(counts, pct = TRUE, digit = 1) {
  ra <- data.table(counts[, 1],
                   apply(X = counts[, -1],
                         MARGIN = 2,
                         FUN = function(a) {
                           a <- a/sum(a)
                           if (pct) a <- 100*a
                           if (!is.na(round)) a <- round(a, digit)
                         }))
  return(ra)
}



#   # Melt
#   ldt.t <- melt.data.table(data = dt.t,
#                            id.vars = 1,
#                            measure.vars = 2:ncol(dt.t),
#                            variable.name = "Sample",
#                            value.name = "Counts")
#   
#   # Relative abundance----
#   dtr.t <- data.table(apply(dt.t[, -1],
#                             2,
#                             function(a) {
#                               a/sum(a)
#                             }))
#   dtr.t$Tax <- dt.t$Tax
#   
#   ldtr.t <- melt.data.table(data = dtr.t,
#                             id.vars = ncol(dtr.t),
#                             measure.vars = 1:(ncol(dtr.t) - 1),
#                             variable.name = "Sample",
#                             value.name = "RA")
#   
#   # Merge counts and relative abundance info----
#   dt2 <- merge(ldt.t,
#                ldtr.t,
#                by = c("Tax",
#                       "Sample"))
#   # Merge counts and sample info----
#   dt2 <- merge(samples,
#                dt2,
#                by = "Sample")
#   
#   # Merge Class info----
#   if (arrg_by != "Class") {
#     tmp <- droplevels(unique(dt1[, c(arrg_by, 
#                                      "Class"),
#                                  with = FALSE]))
#     colnames(tmp)[1] <- "Tax"
#     dt2 <- merge(tmp,
#                  dt2,
#                  by = "Tax")
#   }
#   
#   dt2$Diet <- factor(dt2$Diet)
#   return(dt2)
# }
# 
# # 2. Mean abundance
# mean_ra <- function(dt2, arrg_by, semi_log_x = FALSE, facet_sex = FALSE) {
#   # Mean abundances----
#   mu <- aggregate(dt2$RA,
#                   by = list(Tax = dt2$Tax,
#                             Week = dt2$Week,
#                             Diet = dt2$Diet),
#                   FUN = "mean")
#   
#   if (facet_sex) {
#     mu <- aggregate(dt2$RA,
#                     by = list(Tax = dt2$Tax,
#                               Week = dt2$Week,
#                               Diet = dt2$Diet,
#                               Sex = dt2$Sex),
#                     FUN = "mean")
#   }
#   
#   lvls <- aggregate(dt2$RA,
#                     by = list(Tax = dt2$Tax),
#                     FUN = "mean")
#   lvls <- lvls$Tax[order(lvls$x)]
#   mu$Tax <- factor(mu$Tax,
#                    levels = lvls)
#   return(mu)
# }
# 
# # 3. Plot mean abundances
# ggplot_mean_ra <- function(dt2, arrg_by, semi_log_x = FALSE, facet_sex = FALSE) {
#   # Mean abundances----
#   mu <- mean_ra(dt2 = dt2,
#                 arrg_by = arrg_by, 
#                 semi_log_x = semi_log_x, 
#                 facet_sex = facet_sex)
#   
#   p1 <- ggplot(mu,
#                aes(x = x,
#                    y = Tax,
#                    fill = Diet,
#                    color = Week)) +
#     geom_point(shape = 21,
#                size = 3,
#                alpha = 0.5) +
#     geom_vline(xintercept = 0.01,
#                linetype = "dashed")
#   scale_y_discrete(arrg_by)
#   if (semi_log_x) {
#     p1 <- p1 + scale_x_log10("Relative Abundance")
#   } else {
#     p1 <- p1 + scale_x_continuous("Relative Abundance")
#   }
#   
#   if (facet_sex) {
#     p1 <- p1 + facet_wrap(~ Sex, nrow = 1)
#   }
#   return(p1)
# }