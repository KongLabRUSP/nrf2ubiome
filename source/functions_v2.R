# Functions
# DS, 02/09/2019

# 1. Aggregate data by a taxonomy rank and compute relative abundances
counts_ra_by_tax_rank <- function(dt1, arrg_by) {
  # Counts
  counts <- dt1[, 
                8:ncol(dt1),
                with = FALSE]
  
  dt.t <- lapply(counts,
                 function(a) {
                   out <- aggregate(x = a,
                                    by = list(Tax = dt1[[which(colnames(dt1) == arrg_by)]]),
                                    FUN = "sum")
                 })
  
  dt.t <- Reduce(f = function(...){merge(..., 
                                         by = "Tax",
                                         all = TRUE)},
                 x = dt.t)
  colnames(dt.t)[-1] <- colnames(counts)
  dt.t <- data.table(dt.t)
  
  # Melt
  ldt.t <- melt.data.table(data = dt.t,
                           id.vars = 1,
                           measure.vars = 2:ncol(dt.t),
                           variable.name = "Sample",
                           value.name = "Counts")
  
  # Relative abundance----
  dtr.t <- data.table(apply(dt.t[, -1],
                            2,
                            function(a) {
                              a/sum(a)
                            }))
  dtr.t$Tax <- dt.t$Tax
  
  ldtr.t <- melt.data.table(data = dtr.t,
                            id.vars = ncol(dtr.t),
                            measure.vars = 1:(ncol(dtr.t) - 1),
                            variable.name = "Sample",
                            value.name = "RA")
  
  # Merge counts and relative abundance info----
  dt2 <- merge(ldt.t,
               ldtr.t,
               by = c("Tax",
                      "Sample"))
  # Merge counts and sample info----
  dt2 <- merge(samples,
               dt2,
               by = "Sample")
  
  # Merge Class info----
  if (arrg_by != "Class") {
    tmp <- droplevels(unique(dt1[, c(arrg_by, 
                                     "Class"),
                                 with = FALSE]))
    colnames(tmp)[1] <- "Tax"
    dt2 <- merge(tmp,
                 dt2,
                 by = "Tax")
  }
  
  dt2$Diet <- factor(dt2$Diet)
  return(dt2)
}

# 2. Plot mean abundances
ggplot_mean_ra <- function(dt2, arrg_by, semi_log_x = FALSE, facet_sex = FALSE) {
  # Mean abundances----
  mu <- aggregate(dt2$RA,
                  by = list(Tax = dt2$Tax,
                            Week = dt2$Week,
                            Diet = dt2$Diet),
                  FUN = "mean")
  
  if (!is.null(facet_sex)) {
    mu <- aggregate(dt2$RA,
                    by = list(Tax = dt2$Tax,
                              Week = dt2$Week,
                              Diet = dt2$Diet,
                              Sex = dt2$Sex),
                    FUN = "mean")
  }
  
  lvls <- aggregate(dt2$RA,
                    by = list(Tax = dt2$Tax),
                    FUN = "mean")
  lvls <- lvls$Tax[order(lvls$x)]
  mu$Tax <- factor(mu$Tax,
                   levels = lvls)
  
  p1 <- ggplot(mu,
               aes(x = x,
                   y = Tax,
                   fill = Diet,
                   color = Week)) +
    geom_point(shape = 21,
               size = 3,
               alpha = 0.5) +
    geom_vline(xintercept = 0.01,
               linetype = "dashed")
    scale_y_discrete(arrg_by)
  if (semi_log_x) {
    p1 <- p1 + scale_x_log10("Relative Abundance")
  } else {
    p1 <- p1 + scale_x_continuous("Relative Abundance")
  }
    
  if (facet_sex) {
    p1 <- p1 + facet_wrap(~ Sex, nrow = 1)
  }
  return(p1)
}
