library(readr)
library(tidyverse)
library(viridis)
library(ggthemes)
library(gridExtra)
library(magrittr)

phenotype <- read_delim("~/R/MAGIC_Analysis/GWAS/phenotype.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
summary(phenotype)

score <- tibble(
  n = 491,
  Min = -1.935,
  '1st Qu.' = -0.317,
  Median = 0.330,
  Mean = 0.299,
  '3rd Qu.' = 0.847,
  Max = 3.184,
  "NA's" = 35
)

tbl <- tableGrob(score, rows=NULL)
grid.table(score, rows=NULL )
# save from R graphics device

View(phenotype)

ggplot(phenotype, aes(x = score)) +
  geom_histogram(fill = rev(viridis(20)), na.rm = T, bins = 20) +
  xlab('Root Score') +
  ggtitle("Distribution of 491 Root Scores", subtitle = "Negative = tolerant\nPositive = sensitive\n") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=13,face="bold"),
        axis.title.x = element_text(size=13,face="bold"))

ggsave("Distribution_Root_Scores.tiff", plot = last_plot(), device = "tiff",
       dpi = 'retina')


