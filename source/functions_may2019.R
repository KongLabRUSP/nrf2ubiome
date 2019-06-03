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
## counts_by_tax_rank(dt1 = out.table,
##                    arrg_by = "Phylum")

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
## counts = output from above function (counts_by_tax_rank)
## pct = should the function output percents?
## digit = should the output be rounded, and if yes (non-NA), to what decimal place?
# Output:
## Table of relative abundances at a given taxonomic level (defined by 'counts')
# Example:
## ra_by_tax_rank(counts = counts,
##                pct = TRUE,
##                digit = 2)

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

# 3. Melt relative abundance table and merge with sample info
# Input:
## ra = output from above function (ra_by_tax_rank)
## samples = meta data about samples
## sample_name = name of the column containing sample names
# Output:
## Long table of RA with sample info
# Example:
## ra_melt(ra = ra,
##         samples = samples,
##         sample_name = "SAMPLE_NAME")

ra_melt <- function(ra, samples, sample_name) {
  samples <- data.frame(samples)
  samples[, sample_name] <- as.character(samples[, sample_name])
  
  lra <- melt.data.table(data = ra,
                         id.vars = 1:2,
                         measure.vars = 3:ncol(ra),
                         variable.name = sample_name,
                         value.name = "RA")
  lra[[sample_name]] <- as.character(lra[[sample_name]])
  
  # Merge counts and relative abundance info----
  out <- merge(lra,
               samples,
               by = "SAMPLE_NAME")
  return(out)
}
