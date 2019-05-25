
```{r MeanClass, warning=FALSE,echo=FALSE,message=FALSE,fig.width=10,fig.height=10}
# Plot means
ggplotly(ggplot_mean_ra(dt2, tax, facet_sex = TRUE))

# Plot means, semi-log
ggplotly(ggplot_mean_ra(dt2, tax, TRUE, TRUE))
```

### 2. Order
```{r Order, warning=FALSE,echo=FALSE,message=FALSE,fig.width=10,fig.height=6}
tax <- "Order"
dt2 <- counts_ra_by_tax_rank(dt1, tax)

# Plot samples
p1 <- ggplot(dt2[order(RA,
                       decreasing = TRUE), ],
             aes(x = Sample,
                 y = RA,
                 fill = Tax,
                 colour = Class,
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
```

```{r MeanOrder, warning=FALSE,echo=FALSE,message=FALSE,fig.width=10,fig.height=12}
# Plot means
ggplotly(ggplot_mean_ra(dt2, tax, facet_sex = TRUE))

# Plot means, semi-log
ggplotly(ggplot_mean_ra(dt2, tax, TRUE, TRUE))
```

### 3. Family
```{r Family, warning=FALSE,echo=FALSE,message=FALSE,fig.width=10,fig.height=6}
tax <- "Family"
dt2 <- counts_ra_by_tax_rank(dt1, tax)

# Plot samples
p1 <- ggplot(dt2[order(RA,
                       decreasing = TRUE), ],
             aes(x = Sample,
                 y = RA,
                 fill = Tax,
                 colour = Class,
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


```

```{r MeanFamily, warning=FALSE,echo=FALSE,message=FALSE,fig.width=10,fig.height=15}
# Plot means
ggplotly(ggplot_mean_ra(dt2, tax, facet_sex = TRUE))

# Plot means, semi-log
ggplotly(ggplot_mean_ra(dt2, tax, TRUE, TRUE))
```

### 4. Genus
```{r Genus, warning=FALSE,echo=FALSE,message=FALSE,fig.width=10,fig.height=6}
tax <- "Genus"
dt2 <- counts_ra_by_tax_rank(dt1, tax)

# Plot samples
p1 <- ggplot(dt2[order(RA,
                       decreasing = TRUE), ],
             aes(x = Sample,
                 y = RA,
                 fill = Tax,
                 colour = Class,
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
```

```{r MeanGenera, warning=FALSE,echo=FALSE,message=FALSE,fig.width=10,fig.height=20}
# Plot means
ggplotly(ggplot_mean_ra(dt2, tax, facet_sex = TRUE))

# Plot means, semi-log
ggplotly(ggplot_mean_ra(dt2, tax, TRUE, TRUE))
```
