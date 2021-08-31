library(tidyverse)
library(viridis)

mosaic <- read_tsv("http://mtweb.cs.ucl.ac.uk/mus/www/POOLING/ARABIDOPSIS/FOUNDER/GENOTYPES/mosaic10.lowcov.txt", col_names = T)
 
mosaic <- mosaic %>% 
  separate(col = magic, into = c('x', 'magic_line')) %>% 
  select(-x)

mosaic$magic_line <- as.integer(mosaic$magic_line)

# Just want MAGIC lines 1-527
mosaic <- filter(mosaic, magic_line < 528)

# this shows only 391 of the 527 MAGIC lines are present in the mosaic
n_distinct(mosaic$magic_line)

# phenotype file containing a score for 491 out of the 526 MAGIC lines we received
phenotype <- read_csv("~/R/MAGIC_Analysis/GWAS/Sept_2019/Sept19_updated_figs/phenotype.csv") # update this path

phenotype <- phenotype %>% 
  separate(col = SUBJECT.NAME, into = c('x', 'magic_line')) %>% 
  select(-x)

phenotype$magic_line <- as.integer(phenotype$magic_line)

# this tells you which MAGIC lines in my phenotype file are not found in the mosaic file.  136 total!
left_out <- anti_join(phenotype, mosaic, by = "magic_line")
527-136

# 109 MAGIC lines from my phenotype data aren't being used
sum(!is.na(left_out$score))


####################################################################################################################

# t shows there are MAGIC lines that have up to 19 founders. I thought max was 16.   
t <- mosaic %>% 
  group_by(magic_line) %>% 
  summarise(n_distinct(acc)) %>% 
  rename(Founders = 'n_distinct(acc)')

# t2 shows which MAGIC lines have 19 founders
t2 <- filter(t, Founders == 19)
View(t2)

# t3 shows mosaic data for just MAGIC.13
t3 <- filter(mosaic, magic_line == 13)
View(t3)
# double checking that MAGIC.13 has 19 founders 
n_distinct(t3$acc)

# histogram showing many MAGIC lines have more than 16 founders
ggplot(t, aes(x = t$Founders)) +
  geom_histogram(fill = rev(viridis(13)), na.rm = T, bins = 13) +
  xlab('\nnMAGIC Founders contributing to a MAGIC line ') +
  ylab('Count\n') +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=13,face="bold"),
        axis.title.y = element_text(size=13,face="bold")) +
  scale_x_continuous(breaks=seq(6,19,1)) +
  scale_y_continuous(breaks=seq(0,70,10))

ggsave("magic_genotypes_per_line_problem.tiff", plot = last_plot(), device = "tiff",
       dpi = 600)

