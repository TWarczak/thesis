# [twarczak@polaris ~]$ cd "/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS_2/VARIANT.TABLES/"
  # Set up your directory 
# [twarczak@polaris VARIANT.TABLES]$ wget -r --no-directories "http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/variants.tables/chr1.alleles.txt"
# [twarczak@polaris VARIANT.TABLES]$ wget -r --no-directories "http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/variants.tables/chr2.alleles.txt"
# [twarczak@polaris VARIANT.TABLES]$ wget -r --no-directories "http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/variants.tables/chr3.alleles.txt"
# [twarczak@polaris VARIANT.TABLES]$ wget -r --no-directories "http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/variants.tables/chr4.alleles.txt"
# [twarczak@polaris VARIANT.TABLES]$ wget -r --no-directories "http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/variants.tables/chr5.alleles.txt"
  # download the 5 chromosome variant data
# [twarczak@polaris MAGIC_GWAS_2]$ reconstruction  -a "/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS_2/VARIANT.TABLES/" -m "/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS_2/transposons.gff" -c "/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS" -s "arabidopsis"
  # contruct the mosaic (SNP genomes of all MAGIC lines)
# [twarczak@polaris VARIANT.TABLES]$ PATH=$PATH:"/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS_2/SRC/"
# [twarczak@polaris VARIANT.TABLES]$ genome_scan -f phenotypemodified.txt -p score -n 1000 -w  "/afs/northstar.dartmouth.edu/users/t/twarczak/MAGIC_GWAS_2/VARIANT.TABLES/Sept24_2019/" -t 0
  # use phenotype data to perform association tests
library(tidyverse)
library(qqman)
library(ggrepel)
library(ggpubr)
library(viridis)

MAGIC_Parentals<-c("Bur", "Can", "Col", "Ct", "Edi", "Hi", "Kn", "Ler", "Mt", "No", "Oy", "Po", "Rsch", "Sf", "Tsu", "Wil", "Ws", "Wu", "Zu")

score_logP <- read_delim("~/R/MAGIC_Analysis/GWAS/Sept_2019/Sept19_updated_figs/score.logP.txt", 
                              " ", escape_double = FALSE, trim_ws = TRUE)

## Careful, rs_ID is specific for the SNPs in this subset.  
score_logP <- score_logP %>%
  mutate(SNP = paste0('rs', row_number()), P = 10^-(score_logP$logP))

## Order variables: rs_ID, chr, bp, P-value, -logP-value
score_logP <- score_logP[, c(4, 1:2, 5, 3)]

## qqman package wants class of chr and bp to be integer
score_logP$chr <- as.integer(score_logP$chr)
score_logP$bp <- as.integer(score_logP$bp)

# Skip this filter to get entire manhattan plot
#score_logP_brooklyn_plot <- filter(score_logP, -log10(P)>0.3)

# View most significant SNPs on chromosome 4
#score_logP %>%
  #filter(chr == 4) %>%
  #arrange(desc(logP)) %>%
  #View()

# df for top 20 SNPs of interest on chromosome 4
snpsOfInterest <- filter(score_logP, chr== 4, logP > -log10(0.05/3300000))
locus_size <- max(snpsOfInterest$bp) - min(snpsOfInterest$bp)

27284-26492
10430319-10423508

snpsOfInterest <- filter(score_logP, chr== 4, logP > 10.5) #%>%
                  select(SNP)
                  
# vector of SNP IDs for top 20                   
snpsOfInterest <- snpsOfInterest$SNP
snpsOfInterest

# Build more informative Manhattan plot
BetterGWAS <- score_logP %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(score_logP, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot) %>%
  # Add highlight and annotation information
  mutate(is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  mutate(is_annotate=ifelse(-log10(P)>9.8, "yes", "no")) 

# Prepare X axis
axisdf <- BetterGWAS %>% group_by(chr) %>% summarize(center=(max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
ggplot(BetterGWAS, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = viridis(5), 5) +
  # scale_color_manual(values = rep(c("grey", "skyblue"), 5 )) +
  # custom X axis:
  scale_x_continuous(label = axisdf$chr, 
                     breaks= axisdf$center ) +
  # remove space between plot area and x axis, set y-axis limits
  scale_y_continuous(expand = c(0, 0), limits = c(0, 12.5) ) +     
  labs(x = "Chromosome",
       title = expression(paste("QTL Mapping of Arsenic Tolerance/Sensitivity in ", italic('Arabidopsis thaliana'), ' MAGIC Population'))) +
  # Add highlighted points
  # geom_point(data=subset(BetterGWAS, is_highlight=="yes"), color="orange", size=1.5) +
  # Add label using ggrepel to avoid overlapping
  # geom_label_repel(data=subset(BetterGWAS, is_highlight=="yes"), aes(label=SNP), size=4) +
  # Add horizontal lines for -log(P) = 5 & 7.5
  geom_hline(yintercept= -log10(1e-5), linetype="longdash", color = "#33FF33", size = 0.7)  +
  geom_hline(yintercept= -log10(5*10^-8), linetype= 'longdash', color = "blue", size = 0.65)  +
  geom_hline(yintercept= -log10(0.05/3326146), linetype= 1, color = "red", size = 0.65)  +
  # Customize theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y=element_text(size=14,face="bold"),
    axis.title.x=element_text(size=14,face="bold"),
    axis.text.y=element_text(size = 13),
    axis.text.x=element_text(size = 13),
    plot.title = element_text(size=18,face = 'bold', margin = margin(t= -5)))


ggsave("Manhattan_MAGIC_3.3M.tiff", plot = last_plot(), device = "tiff",
       dpi = 'retina')
-log10(5*10^-8) # Conventional alpha for GWAS
-log10(1e-5) # Conventional lower confidence interval for GWAS
-log10(0.05/3300000) # alpha value adjusted with the Bonferroni correction.  a = 0.05 goes to a = (0.05/k) = (0.05/3,300,000)

# ZOOM in on chr 4 QTL
x <- ggplot(BetterGWAS, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=2.3) +
  scale_color_manual(values = viridis(5), 5) +
  # scale_color_manual(values = rep(c("grey", "skyblue"), 5 )) +
  # custom X axis:
  scale_x_continuous(label= axisdf$chr,
                     breaks= axisdf$center) +
  # remove space between plot area and x axis, set y-axis limits
  scale_y_continuous(expand = c(0, 0)) +     
  # Add horizontal lines for -log(P) = 5 & 7.5
  geom_hline(yintercept= -log10(0.05/3326146), linetype= 1, color = "red", size = 0.65)  +
  # Customize theme:
  labs(x = "Chromosome") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y=element_text(size=14,face="bold"),
    axis.title.x=element_text(size=14,face="bold"),
    axis.text.x=element_text(size = 13),
    axis.text.y=element_text(size = 13))
  
x2 <- x + 
  coord_cartesian(xlim = c(78585303, 92570326), ylim = c(7.5,13)) +
  geom_text_repel(data          = subset(BetterGWAS, is_highlight=="yes"), aes(label=SNP),
                  nudge_y       = 13 - subset(BetterGWAS, is_highlight=="yes")$logP,
                  size          = 4,
                  box.padding   = 1.5,
                  point.padding = 0.5,
                  force         = 10,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x")

ggsave("chr4_QTL_Zoom.tiff", plot = last_plot(), device = "tiff",
       dpi = 'retina')

############

snpsOfInterest_QTL3 <- filter(score_logP, chr == 3, logP > 8) %>%
  select(SNP)

snpsOfInterest_QTL3 <- snpsOfInterest_QTL3$SNP
snpsOfInterest_QTL3

# Just to highlight variants in chr3 QTL
BetterGWAS_QTL3 <- score_logP %>% 
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(score_logP, ., by=c("chr"="chr")) %>%
  arrange(chr, bp) %>%
  mutate(BPcum=bp+tot) %>%
  mutate(is_highlight=ifelse(SNP %in% snpsOfInterest_QTL3, "yes", "no"))

axisdf_QTL3 <- BetterGWAS_QTL3 %>% group_by(chr) %>% summarize(center=(max(BPcum) + min(BPcum) ) / 2 )

x3 <- ggplot(BetterGWAS_QTL3, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=3) +
  scale_color_manual(values = viridis(5), 5) +
  # scale_color_manual(values = rep(c("grey", "skyblue"), 5 )) +
  # custom X axis:
  scale_x_continuous(label= axisdf$chr,
                     breaks= axisdf$center) +
  # remove space between plot area and x axis, set y-axis limits
  scale_y_continuous(expand = c(0, 0)) +     
  # Add horizontal lines for -log(P) = 5 & 7.5
  geom_hline(yintercept= -log10(0.05/3326146), linetype= 1, color = "red", size = 0.65)  +
  # Customize theme:
  labs(x = "Chromosome") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y=element_text(size=14,face="bold"),
    axis.title.x=element_text(size=14,face="bold"),
    axis.text.x=element_text(size = 13),
    axis.text.y=element_text(size = 13))

# QTL3 Zoom
QTLch3 <- x3 + 
  coord_cartesian(xlim = c(47980240, 60980240), ylim = c(5.0,9.5)) +
  geom_text_repel(data          = subset(BetterGWAS_QTL3, is_highlight=="yes"), aes(label=SNP),
                  nudge_y       = 9.5 - subset(BetterGWAS_QTL3, is_highlight=="yes")$logP,
                  size          = 4,
                  box.padding   = 1.5,
                  point.padding = 0.5,
                  force         = 10,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x")

ggsave("chr3_QTL_Zoom.tiff", plot = last_plot(), device = "tiff",
       dpi = 'retina')
#############################################################################
# Table of SNPs of interest

top_variants_QTL4 <- score_logP %>% 
  filter(SNP == "rs2408584"| SNP == "rs2412612" | SNP == "rs2412782" | SNP == "rs2412888" | SNP == "rs2420041"
         | SNP == "rs2420043"| SNP == "rs2423698")

top_variants_QTL3 <- score_logP %>% 
  filter(SNP == "rs1438314"| SNP == "rs1451683" | SNP == "rs1451684" | SNP == "rs1451685" | SNP == "rs1451715"
         | SNP == "rs1451716"| SNP == "rs1451722")
# don't know why this doesn't work anymore
# rename(top_variants_QTL3, 'Variant' = SNP, '-logP' = logP)

top_variants_QTL4 <- top_variants_QTL4 %>% 
  select(SNP, chr, bp, logP)
  
top_variants_QTL4 <- top_variants_QTL4 %>%   
  arrange(desc(logP))

library(gridExtra)

# Set theme to allow for plotmath expressions
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(top_variants_QTL4, rows=NULL)
# Plot chart and table into one object
# grid.arrange(x2, tbl, ncol=2)

# or just plot table
grid.table(top_variants_QTL4, rows=NULL)
# save from R graphics device

  #=============================================================================================================================================

MAGIC_Parentals<-c("Bur", "Can", "Col", "Ct", "Edi", "Hi", "Kn", "Ler", "Mt", "No", "Oy", "Po", "Rsch", "Sf", "Tsu", "Wil", "Ws", "Wu", "Zu")
phenotype<-c(8, 1, 12, 6, 19, 13, 18, 17, 16, 15, 2, 5, 4, 14, 9, 11, 7, 10, 3) #rank from most sensitive (1) to most tolerant (19)

NIP1_Expression_As10_14d<-c(2, 19, 4, 16, 1, 6, 5, 9, 17, 3, 18, 13, 8, 7, 12, 10, 11, 15, 14) # ranked NIP1 expression from high (1) to low (19)
NIP1_Expression_Control_7d <- c(16, 19, 11, 17, 1, 9, 4, 12, 7, 3, 14, 8, 13, 2, 18, 10, 5, 15, 6)

MAGIC_matrix <- matrix(c(phenotype, NIP1_Expression_As10_14d), nrow=19, byrow = FALSE, dimnames = list(MAGIC_Parentals))
colnames(MAGIC_matrix)<- c("Phenotype_Rank", "NIP1_Exp_As10")
rownames(MAGIC_matrix)<-MAGIC_Parentals

MAGIC_table<-data.frame(MAGIC_matrix)

MAGIC_table2 <- cbind(MAGIC_table, NIP1_Expression_Control_7d)
colnames(MAGIC_table2)<- c("Phenotype_Rank", "NIP1_Exp_As10", "NIP1_Exp_Con_7d")


ggscatter(MAGIC_table2, x = "Phenotype_Rank", y = "NIP1_Exp_As10", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coef.coord = c(10, 19), cor.method = "spearman",
          xlab = "Phenoptype Rank \n(Sensitive-Tolerant)", ylab = "NIP1_Expression", label = MAGIC_Parentals)

ggscatter(MAGIC_table2, x = "Phenotype_Rank", y = "NIP1_Exp_Con_7d", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coef.coord = c(10, 19), cor.method = "spearman",
          xlab = "Phenoptype Rank \n(Sensitive-Tolerant)", ylab = "NIP1_Expression", label = MAGIC_Parentals)

#=====================================================================================================================================================

# Measure expression of NIP1 after 7d, then expose plants to arsenic, measure NIP1 expression every day for 5 days
MAGIC_10As_Assay <- read_csv("~/R/qPCR/March2019/MAGIC_10As_Assay.csv")

# Only want 3 variables and don't need Col 15As data
Tidy_MAGIC <- MAGIC_10As_Assay %>%
  select(Target, Sample, Cq) %>%  
  filter(!grepl('15', Sample))    

# Col data from Feb was in different Excel format; matching to the new qPCR data
Tidy_MAGIC$Sample[649:654]="Col-Con"
Tidy_MAGIC$Sample[655:660]="Col-1dAs"
Tidy_MAGIC$Sample[661:666]="Col-2dAs"
Tidy_MAGIC$Sample[667:672]="Col-3dAs"
Tidy_MAGIC$Sample[673:678]="Col-4dAs"
Tidy_MAGIC$Sample[679:684]="Col-5dAs"

# Separating the Sample variable into Parental_ID and Condition variables
Tidy_MAGIC <- separate(Tidy_MAGIC, col = Sample, into = c("Parental_ID", "Condition"), sep = "-")

Tidy_MAGIC$Target <- as.factor(Tidy_MAGIC$Target)
Tidy_MAGIC$Condition <- as.factor(Tidy_MAGIC$Condition)
Tidy_MAGIC$Parental_ID <- as.factor(Tidy_MAGIC$Parental_ID)  # for plot at end

levels(Tidy_MAGIC$Condition)[levels(Tidy_MAGIC$Condition)=="Con"] <- "0dAs"
levels(Tidy_MAGIC$Target)[levels(Tidy_MAGIC$Target)=="UBQ1"] <- "UBQ"   # Used UBQ1 in the Col data

# Assign levels to Condition factor
Tidy_MAGIC$Condition <- factor(Tidy_MAGIC$Condition, levels = c("0dAs", "1dAs", "2dAs", "3dAs", "4dAs", "5dAs"))

# Averaging replicates 
Tidy_MAGIC <- arrange(Tidy_MAGIC, Parental_ID, Condition, Target, Cq) %>%
  group_by(Parental_ID, Condition, Target) %>%
  summarise(avg_Cq = mean(Cq))
# Want to compare UBQ and NIP values
Tidy_MAGIC <- spread(Tidy_MAGIC, Target, avg_Cq)

# Data Cleaning: removing obvious outlier reps, averaging two remaining reps
Tidy_MAGIC$NIP1[65] <- mean(19.46345, 19.13845)  # Oy 4dAs NIP1
Tidy_MAGIC$UBQ[66] <- mean(16.04281, 16.11887)   # Oy 5dAs UBQ
Tidy_MAGIC$NIP1[6] <- mean(20.11393, 19.69400)   # Bur 5dAs NIP1
Tidy_MAGIC$UBQ[4] <- mean(16.41937, 16.36203)    # Bur 3dAs UBQ
Tidy_MAGIC$NIP1[4] <- mean(20.27851, 21.53514)   # Bur 3dAs NIP1
Tidy_MAGIC$UBQ[41] <- mean(16.44878, 16.71500)   # Kn 4dAs UBQ
Tidy_MAGIC$NIP1[41] <- mean(20.77655, 20.89758)  # Kn 4dAs NIP1
Tidy_MAGIC$UBQ[88] <- mean(16.33316, 15.96877)   # Tsu 3dAs UBQ
Tidy_MAGIC$NIP1[88] <- mean(19.47061, 19.54333)  # Tsu 3dAs NIP1
Tidy_MAGIC$NIP1[10] <- mean(19.08664, 18.78850)  # Can 3dAs NIP1
Tidy_MAGIC$NIP1[23] <- mean(19.65905, 19.74175)  # Ct 4dAs NIP1
Tidy_MAGIC$UBQ[24] <- mean(15.69372, 16.09460)   # Ct 5dAs UBQ
Tidy_MAGIC$UBQ[97] <- mean(15.41159, 15.39738)   # Ws 0dAs UBQ
Tidy_MAGIC$NIP1[106] <- mean(19.65544, 19.28824) # Wu 3dAs NIP1

# Calculate fold change of NIP1 expression
Tidy_MAGIC <- mutate(Tidy_MAGIC,NIP1_Exp = 2^(UBQ - NIP1)) 
# Exp. relative to Col 0dAs
Tidy_MAGIC <- mutate(Tidy_MAGIC, NIP1_Exp_Norm = NIP1_Exp / Tidy_MAGIC$NIP1_Exp[13])  

# Order levels of Parental_ID by sensitivity to arsenic 
Tidy_MAGIC$Parental_ID <- factor(Tidy_MAGIC$Parental_ID, levels = c("Can","Oy","Zu","Rsch","Po","Ct","Ws","Bur","Tsu","Wu","Wil","Col","Hi","Sf","No","Mt","Ler","Kn","Edi"))

ggplot(Tidy_MAGIC, aes(x = Condition, y = NIP1_Exp_Norm, color = Parental_ID)) +
  geom_line(aes(group = Parental_ID)) +
  geom_point() +
  facet_wrap(. ~ Parental_ID) +
  theme(legend.position = "none") +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  labs(x = "Days Exposed to arsenic") 


