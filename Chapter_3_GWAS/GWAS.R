# [twarczak@polaris ~]$ cd "/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS_2/VARIANT.TABLES/"
# [twarczak@polaris VARIANT.TABLES]$ PATH=$PATH:"/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS_2/SRC/"
# [twarczak@polaris VARIANT.TABLES]$ genome_scan -f phenotypemodified.txt -p score -n 1000 -w  "/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS_2/VARIANT.TABLES/Apr10/" -t 0.2

library(tidyverse)
library(qqman)
library(ggrepel)
library(ggpubr)

MAGIC_Parentals<-c("Bur", "Can", "Col", "Ct", "Edi", "Hi", "Kn", "Ler", "Mt", "No", "Oy", "Po", "Rsch", "Sf", "Tsu", "Wil", "Ws", "Wu", "Zu")

score_annotated <- read_delim("~/R/MAGIC_Analysis/GWAS/April_10/score.annotated.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

## Careful, rs_ID is specific for the SNPs in this subset.  i.e. top 10000 logP values
score_annotated <- score_annotated %>%
  subset(select = -c(trait, logP))

SNP_top20_2 <- left_join(SNP_top20, score_annotated, by = c('chr', 'bp')) %>%
  subset(select = -c(is_highlight, is_annotate, tot, BPcum))

score_annotated <- score_annotated %>%
  mutate(SNP = paste0('rs', row_number()), P = 10^-(score_annotated$logP))

## Order variables: rs_ID, chr, bp, P-value, -logP-value
score_annotated <- score_annotated[, c(4, 1:2, 5, 3)]

## qqman package wants class of chr and bp to be integer
score_annotated$chr <- as.integer(score_annotated$chr)
score_annotated$bp <- as.integer(score_annotated$bp)

glimpse(score_annotated)
score_annotated_Apr11 <- filter(score_annotated, -log10(P)>0.3)

manhattan(score_annotated_Apr11, chr  = 'chr', bp = 'bp')

# QC for GWAS  qq(score_annotated$P)

# View most significant SNPs on chromosome 4
score_annotated_Apr11 %>%
  filter(chr == 4) %>%
  arrange(desc(logP)) %>%
  View()

# df for top 20 SNPs of interest on chromosome 4
snpsOfInterest <- filter(score_annotated, logP > 9.8) %>%
  select(SNP)

# vector of SNP IDs for top 20                   
snpsOfInterest <- snpsOfInterest$SNP

# Build more informative Manhattan plot
BetterGWAS <- score_annotated_Apr11 %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(score_annotated_Apr11, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot) %>%
  
  # Add highlight and annotation information
  mutate(is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  mutate(is_annotate=ifelse(-log10(P)>9.8, "yes", "no") ) 

# Prepare X axis
axisdf <- BetterGWAS %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
ggplot(BetterGWAS, aes(x=BPcum, y=-log10(P))) +
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 5 )) +
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  labs(x = "Chromosome") +
  # Add highlighted points
  geom_point(data=subset(BetterGWAS, is_highlight=="yes"), color="orange", size=1.5) +
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=subset(BetterGWAS, is_annotate=="yes"), aes(label=SNP), size=2) +
  # Add horizontal lines for -log(P) = 5 & 7.5
  geom_hline(yintercept= -log10(1e-5), linetype="dashed", color = "red")  +
  geom_hline(yintercept= -log10(5*10^-8), linetype= 1, color = "blue")  +
  # Customize theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank() 
  )

#==========================================================================================================
SNP_top20 <- BetterGWAS %>%
  filter(is_annotate == 'yes')