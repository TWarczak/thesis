# [twarczak@polaris ~]$ cd "/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS_2/VARIANT.TABLES/"
# [twarczak@polaris VARIANT.TABLES]$ PATH=$PATH:"/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS_2/SRC/"
# [twarczak@polaris VARIANT.TABLES]$ genome_scan -f phenotypemodified.txt -p score -n 100 -w  "/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS_2/VARIANT.TABLES/Apr10/" -t 0.2

library(ggplot2)
library(dplyr)
library(qqman)
library(readr)
library(ggrepel)score_annotated <- read_delim("~/R/MAGIC_Analysis/GWAS/score.annotated2.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)


## Adding a rs_ID for each SNP/marker, eliminate 'trait' variable
## Careful, rs_ID is specific for the SNPs in this subset.  i.e. top 10000 logP values
score_annotated <- score_annotated %>%
  mutate(SNP = paste0('rs', row_number()), P = 10^-(score_annotated$logP)) %>%
  subset(select = -trait)

## Order variables: rs_ID, chr, bp, P-value, -logP-value
score_annotated <- score_annotated[, c(23, 1:2, 24, 3, 4:22)]

## qqman package wants class of chr and bp to be integer
score_annotated$chr <- as.integer(score_annotated$chr)
score_annotated$bp <- as.integer(score_annotated$bp)

str(score_annotated)

#all 5 chromosomes, labels for only logP > 10.3
ggplot(score_annotated, aes(x=bp, y=logP, color= factor(chr))) + 
  geom_point() + 
  facet_wrap(~chr, nrow=1) +
  theme(legend.position = "none") +
  geom_text(aes(label=ifelse(logP>10.3,as.character(rsID),'')),hjust=0,vjust=0, position=position_jitter(width=0,height=0))

# Chr4 w/ top labels
score_annotated %>%
  filter(chr == 4) %>%
  arrange(desc(logP)) %>%
  View()

score_annotated_chr4 <- score_annotated %>%
  filter(chr == 4)

ggplot(score_annotated_chr4, aes(x=bp, y=logP, color= chr)) + 
  geom_point() + 
  theme(legend.position = "none") +
  geom_text_repel(aes(label=ifelse(logP>10.2,as.character(rsID),'')))

# Chr3 w/ top labels
score_annotated %>%
  filter(chr == 3) %>%
  arrange(desc(logP)) %>%
  View()

score_annotated_chr3 <- score_annotated %>%
  filter(chr == 3)

ggplot(score_annotated_chr3, aes(x=bp, y=logP, color= chr)) + 
  geom_point() + 
  theme(legend.position = "none") +
  geom_text_repel(aes(label=ifelse(logP>7,as.character(rsID),'')))



manhattan(score_annotated, bp= "bp", chr = 'chr', p= "logP", snp = "bp",  logp = F,  main = "MAGIC AsIII GWAS", ylim = c(0, 15), cex = 0.6, 
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F)
