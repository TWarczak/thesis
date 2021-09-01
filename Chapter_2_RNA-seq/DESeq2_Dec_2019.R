library(tidyverse)
BiocManager::install('SummarizedExperiment')
BiocManager::install("DESeq2")
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(patchwork)
library(BiocManager)
library(biomaRt)
library(grid)
library(venn)


FACS_Todd_counts<- read_delim("C:/Users/Todd/OneDrive for Business/Todd_RNAseq/July_2018_Data/sharing_Data_RNA_Guerinot/sharing_Data_RNA_Guerinot2/RNAseq_WER_Cortex_SCR_genes_count_July2018.txt", delim = "\t")
glimpse(FACS_Todd_counts)

celltype <- c("WER", "WER", "WER", "WER", "WER", "WER", "Cortex", "Cortex", "Cortex", "Cortex", "Cortex", "Cortex",
              "SCR", "SCR", "SCR", "SCR", "SCR", "SCR")
condition <- c("0_As", "0_As","0_As","50_As","50_As","50_As","0_As","0_As","0_As","50_As","50_As","50_As","0_As","0_As","0_As",
               "50_As","50_As","50_As")

FACS_metadata <- data.frame(celltype, condition)

tracking_id <-FACS_Todd_counts$tracking_id

FACS_Todd_counts <- FACS_Todd_counts[, -1]

FACS_Todd_counts <- as.matrix.data.frame(round(FACS_Todd_counts))

rownames(FACS_Todd_counts) <- tracking_id
rownames(FACS_metadata) <- colnames(FACS_Todd_counts)
# Important for DESeq2 that rownames of metadata == colnames of countdata

class(FACS_Todd_counts) #make sure its a count matrix

# Create a DESeq2 object
# DESeq object to test for the effect of celltype on the effect of arsenic
dds_FACS <- DESeqDataSetFromMatrix(countData = FACS_Todd_counts,
                                             colData = FACS_metadata,
                                             design = ~ condition + celltype + condition:celltype)

# Determine the size factors to use for normalization
dds_FACS <- estimateSizeFactors(dds_FACS)

# Extract the normalized counts
#norm_counts_dds_FACS <- round(counts(dds_FACS, normalized =T))
norm_counts_dds_FACS <- counts(dds_FACS, normalized =T)

myTable <- tableGrob(filter(FACS_norm_counts, X1 == 'AT4G14880'))

grid.draw(myTable)
########################################################################
# Hierarchical Clustering with correlation heatmaps
vsd_FACS <- vst(dds_FACS, blind = TRUE)

vsd_mat <- assay(vsd_FACS)
vsd_cor <- cor(vsd_mat)
View(vsd_cor)
pheatmap <- pheatmap(vsd_cor, annotation = dplyr::select(FACS_metadata, condition))
# ggsave("hier_clus_cor_heatmap.tiff", plot = pheatmap, device = "tiff",
#       dpi = 'retina')

# PCA Analysis
# Transform the normalized counts   
# Create the PCA plot for PC1 and PC2 and color by condition
pcaData <- plotPCA(vsd_FACS, intgroup=c("condition", "celltype"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=celltype, shape=condition)) +
  geom_point(size=4) +
  scale_shape_manual(values=c(15, 17)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# ggsave("PCA.tiff", plot = last_plot(), device = "tiff",
#       dpi = 'retina')


#########################################################################
# Need to make groups to use in comparisons. i.e. the 3 reps of WER-As become the group: WER0_As
dds_FACS$group <- factor(paste0(dds_FACS$celltype, dds_FACS$condition))
design(dds_FACS) <- ~ group
?DESeq
# Run Differential Expression analysis
# This uses a negative binomial model (nbGLM or Gamma-Poisson distribution), 
# K_ij ~ NB(mu_ij, α_i), K = raw count for i (gene) j (sample); mu_ij = s_i x q_ij, s_i = sample specific size factor, q_ij = normalized count parameter; α_i = gene-specific dispersion parameter
# It accounts for the additional variation in the data added by the small number of biol. reps. 
# K_ij for each gene is modeled using the fitted mean mu_ir and αi, as input to the neg binomial model to fit the raw count data
# log2q_ij = Σ_rX_jrß_ir ; Σ_rX_jr = sum of the samples of conditional interest, ß_ir = log2 foldchange estimates.   ΣrXjrßir = log2 foldchange between conditions
# For each gene, the model uses the log2 normalized counts (log2qij), to determine the log2 foldchange between conditions (ΣrXjrßir) 
# DESeq2 will perform the Wald test for pairwise comparisons to test for differences in expression between two sample groups for the condition of interest
# My sample groups for condition are 0_As and 50_As
# Results of Wald test can be extracted with the results() function
dds_FACS_DESeq <- DESeq(dds_FACS)
resultsNames(dds_FACS_DESeq)

# Check to see if dispersions decrease with increasing mean and raw dispersions seem to cluster around the max likelihood line.
plotDispEsts(dds_FACS_DESeq)
# Looks good      save from plot

# Make contrasts between conditions of celltype groups.  i.e WER50_As_vs_WER0_As 
# Select significant genese with padj < 0.05
results_WER <- results(dds_FACS_DESeq, contrast=c("group","WER50_As","WER0_As"), alpha = 0.05) 
results_SCR <- results(dds_FACS_DESeq, contrast=c("group","SCR50_As","SCR0_As"), alpha = 0.05) 
results_Cortex <- results(dds_FACS_DESeq, contrast=c("group","Cortex50_As","Cortex0_As"), alpha = 0.05) 

mcols(results_WER, use.names=T)
# baseMean: mean of normalized counts for all samples
# log2FoldChange: log2 fold change
# lfcSE: standard error of the log2foldchange 
# stat: Wald statistic (log2fc/lfcSE), used to generate a two-tailed pvalue
# pvalue: Wald test p-value
# padj: p-values corrected for multiple testing using Benjamini and Hochberg method
summary(results_WER)
summary(results_Cortex)
summary(results_SCR)
# To generate more accurate log2 foldchange estimates, use Shrunken log2 foldchanges (LFC)
# for rationale go to:
# https://hbctraining.github.io/DGE_workshop_salmon/lessons/05_DGE_DESeq2_analysis2.html
# produces the same DESeqResults object, but with the log2FoldChange and lfcSE columns replaced with shrunken (more accurate) results.  baseMean, stat, pvalue, & padj unchanged

results_WER_shr <- lfcShrink(dds_FACS_DESeq, contrast=c("group","WER50_As","WER0_As"), res=results_WER)
results_SCR_shr <- lfcShrink(dds_FACS_DESeq, contrast=c("group","SCR50_As","SCR0_As"), res=results_SCR)
results_Cortex_shr <- lfcShrink(dds_FACS_DESeq, contrast=c("group","Cortex50_As","Cortex0_As"), res=results_Cortex)
summary(results_WER_shr)

plotMA(results_WER_shr, alpha = 0.05, ylim = c(-10,15), xlab="",ylab="", frame.plot = F)   
  mytitle = "MA plot, Epidermal Genes" 
  mysubtitle = "0 AsIII -> 50 AsIII" 
  mtext(side=3, line=-1, at=0.07, font= 2, adj=0, cex=1.3, mytitle) 
  mtext(side=3, line=-2, at=0.07, adj=0, cex=0.9, mysubtitle) 
  mtext(side=1, line=2, "Mean of Normalized Counts", font=2,cex=1.1) 
  mtext(side=2, line=3, "Log2 Fold Change", font=2, cex=1.1) 
# save from plot window
plotMA(results_SCR_shr, alpha = 0.05, ylim = c(-10,15), xlab="",ylab="", frame.plot = F) 
  mytitle = "MA plot, Endodermal Genes"
  mysubtitle = "0 AsIII -> 50 AsIII"
  mtext(side=3, line=-1, at=0.07, font= 2, adj=0, cex=1.3, mytitle)
  mtext(side=3, line=-2, at=0.07, adj=0, cex=0.9, mysubtitle)
  mtext(side=1, line=2, "Mean of Normalized Counts", font=2,cex=1.1)
  mtext(side=2, line=3, "Log2 Fold Change", font=2, cex=1.1)
# save from plot window
plotMA(results_Cortex_shr, alpha = 0.05, ylim = c(-10,15), xlab="",ylab="", frame.plot = F) 
  mytitle = "MA plot, Cortex Genes"
  mysubtitle = "0 AsIII -> 50 AsIII"
  mtext(side=3, line=-1, at=0.07, font= 2, adj=0, cex=1.3, mytitle)
  mtext(side=3, line=-2, at=0.07, adj=0, cex=0.9, mysubtitle)
  mtext(side=1, line=2, "Mean of Normalized Counts", font=2,cex=1.1)
  mtext(side=2, line=3, "Log2 Fold Change", font=2, cex=1.1)
# save from plot window

# Obtain logical vector regarding whether padj values are less than 0.05 
results_WER2 <- data.frame(results_WER_shr) %>%
  rownames_to_column(var = "GeneID") %>%
  mutate(threshold = padj < 0.05)
results_SCR2 <- data.frame(results_SCR_shr) %>%
  rownames_to_column(var = "GeneID") %>%
  mutate(threshold = padj < 0.05)
results_Cortex2 <- data.frame(results_Cortex_shr) %>%
  rownames_to_column(var = "GeneID") %>%
  mutate(threshold = padj < 0.05)

# Volcano plots
p_W <- ggplot(results_WER2) +
        geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold), alpha = 0.5, size = 1.5) +
        theme_minimal() +
        labs(x="", y ="-log10 adjusted p-value", title = "Volcano plot, Epidermal Genes", subtitle =  "0 AsIII -> 50 AsIII") +
        theme(legend.position = "none",
          plot.title = element_text(size = rel(1.2), face = 'bold'),
          axis.title = element_text(size = rel(1), face = 'bold'))
p_S <- ggplot(results_SCR2) +
        geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold), alpha = 0.5, size = 1.5) +
        theme_minimal() +
        labs(x="", y ="", title = "Volcano plot, Endodermal Genes", subtitle =  "0 AsIII -> 50 AsIII") +
        theme(legend.position = "none",
          plot.title = element_text(size = rel(1.2), face = 'bold'),
          axis.title = element_text(size = rel(1), face = 'bold'))
p_C <- ggplot(results_Cortex2) +
        geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold), alpha = 0.5, size = 1.5) +
        theme_minimal() +
        labs(x="log2 fold change", y ="", title = "Volcano plot, Cortex Genes", subtitle =  "0 AsIII -> 50 AsIII") +
        theme(legend.position = "none",
          plot.title = element_text(size = rel(1.2), face = 'bold'),
          axis.title = element_text(size = rel(1), face = 'bold'))

p_W + p_C + p_S
#ggsave("volcano_all.tiff", plot = last_plot(), device = "tiff",
#       dpi = 'retina')

###############################################################################

# create tibbles of the shrunken results
res_WER_tb <- results_WER_shr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
res_SCR_tb <- results_SCR_shr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
res_Cortex_tb <- results_Cortex_shr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# filter out non sig & abs(log2FC) < 1 genes
sig_WER <- res_WER_tb %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)
sig_SCR <- res_SCR_tb %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)
sig_Cortex <- res_Cortex_tb %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

#################################################################
FPKM_complete <- read_csv("~/R/DESeq2_Todd_RNAseq/FPKM_complete.csv")
glimpse(FPKM_complete)
View(FPKM_complete)

FPKM_avg <- FPKM_complete %>% 
  transmute(gene_id = gene_id,
            WER_minus_As = (FPKM_complete$Wer_minus_AS_0 + FPKM_complete$Wer_minus_AS_1 + FPKM_complete$Wer_minus_AS_2)/3,
            WER_plus_As = (FPKM_complete$Wer_plus_AS_0 + FPKM_complete$Wer_plus_AS_1 + FPKM_complete$Wer_plus_AS_2)/3,
            Cortex_minus_As = (FPKM_complete$Cortex_minus_AS_0 + FPKM_complete$Cortex_minus_AS_1 + FPKM_complete$Cortex_minus_AS_2)/3, 
            Cortex_plus_As = (FPKM_complete$Cortex_plus_AS_0 + FPKM_complete$Cortex_plus_AS_1 + FPKM_complete$Cortex_plus_AS_2)/3,  
            SCR_minus_As = (FPKM_complete$SCR_minus_AS_0 + FPKM_complete$SCR_minus_AS_1 + FPKM_complete$SCR_minus_AS_2)/3, 
            SCR_plus_As = (FPKM_complete$SCR_plus_AS_0 + FPKM_complete$SCR_plus_AS_1 + FPKM_complete$SCR_plus_AS_2)/3)

#######################
WER_up_FPKM <- sig_WER %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange > 1.3) %>% 
  inner_join(FPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, WER_minus_As, WER_plus_As) %>% 
  filter(WER_plus_As > 8) %>% 
  dplyr::rename(WER_0As_FPKM = WER_minus_As, WER_50As_FPKM = WER_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   WER_0As_FPKM = round(WER_0As_FPKM, digits = 2), 
                   WER_50As_FPKM = round(WER_50As_FPKM, digits = 2))

WER_down_FPKM <- sig_WER %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange < -1.3) %>% 
  inner_join(FPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, WER_minus_As, WER_plus_As) %>% 
  filter(WER_minus_As > 8) %>% 
  dplyr::rename(WER_0As_FPKM = WER_minus_As, WER_50As_FPKM = WER_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   WER_0As_FPKM = round(WER_0As_FPKM, digits = 2), 
                   WER_50As_FPKM = round(WER_50As_FPKM, digits = 2))

SCR_up_FPKM <- sig_SCR %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange > 1.3) %>% 
  inner_join(FPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, SCR_minus_As, SCR_plus_As) %>% 
  filter(SCR_plus_As > 8) %>% 
  dplyr::rename(SCR_0As_FPKM = SCR_minus_As, SCR_50As_FPKM = SCR_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   SCR_0As_FPKM = round(SCR_0As_FPKM, digits = 2), 
                   SCR_50As_FPKM = round(SCR_50As_FPKM, digits = 2))

SCR_down_FPKM <- sig_SCR %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange < -1.3) %>% 
  inner_join(FPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, SCR_minus_As, SCR_plus_As) %>% 
  filter(SCR_minus_As > 8) %>% 
  dplyr::rename(SCR_0As_FPKM = SCR_minus_As, SCR_50As_FPKM = SCR_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   SCR_0As_FPKM = round(SCR_0As_FPKM, digits = 2), 
                   SCR_50As_FPKM = round(SCR_50As_FPKM, digits = 2))

Cortex_up_FPKM <- sig_Cortex %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange > 1.3) %>% 
  inner_join(FPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, Cortex_minus_As, Cortex_plus_As) %>% 
  filter(Cortex_plus_As > 8) %>% 
  dplyr::rename(Cortex_0As_FPKM = Cortex_minus_As, Cortex_50As_FPKM = Cortex_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   Cortex_0As_FPKM = round(Cortex_0As_FPKM, digits = 2), 
                   Cortex_50As_FPKM = round(Cortex_50As_FPKM, digits = 2))

Cortex_down_FPKM <- sig_Cortex %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange < -1.3) %>% 
  inner_join(FPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, Cortex_minus_As, Cortex_plus_As) %>% 
  filter(Cortex_minus_As > 8) %>% 
  dplyr::rename(Cortex_0As_FPKM = Cortex_minus_As, Cortex_50As_FPKM = Cortex_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   Cortex_0As_FPKM = round(Cortex_0As_FPKM, digits = 2), 
                   Cortex_50As_FPKM = round(Cortex_50As_FPKM, digits = 2))


# venn = 
#   list(WER_up_FPKM$gene_id, Cortex_up_FPKM$gene_id, SCR_up_FPKM$gene_id, WER_down_FPKM$gene_id, Cortex_down_FPKM$gene_id, SCR_down_FPKM$gene_id)
#     
# tiff('ven5.tiff', width = 1200, height = 1200)
# 
# venn.result =
#   venn(venn, ilabels = T, 
#        zcolor = "style", size = 25, cexil = 1.3, cexsn = 1.3,
#        snames = c('Epiderm_Up', 'Cortex_Up', 'Endoderm_Up', 'Epiderm_Down', 'Cortex_Down', 'Endoderm_Down'));
# 
# dev.off()

###################################################

# write_delim(WER_up_FPKM$gene_id, '~\WER_up_geneID.csv', 
#             delim = " ", na = "NA", append = FALSE,col_names = F, quote_escape = "double")
# 
# path = '~/R/DESeq2_Todd_RNAseq'
# write.csv(WER_up_FPKM,'WER_up_geneID.csv',sep = '')

##########################
# FPKM_density plot, faceted by up/down regulation & cell-type
myColors <- c('#4E84C4', 'tomato')
my_title <-"FPKM Values of Differentially Expressed Genes"

WER_FPKM_hist <- bind_rows(WER_up_FPKM, WER_down_FPKM) %>%
  dplyr::rename('0As'= WER_0As_FPKM, '50As'= WER_50As_FPKM) %>% 
  pivot_longer(
    cols = contains("0As"),
    names_to = "Condition",
    values_to = "FPKM",
    values_drop_na = TRUE)

WER_FPKM_hist <- WER_FPKM_hist %>%     
  mutate(As_regulated = ifelse(log2FoldChange > 0, "Up-regulated", "Down-regulated"))

WER_den_FPKM <- ggplot(WER_FPKM_hist, aes(x= FPKM)) + 
  geom_density(aes( fill = Condition),alpha=0.5) +  
  labs(y='Density', title =my_title, subtitle =expression(bolditalic("Epidermal Cells"))) +
  scale_x_continuous(breaks=c(0,1,2,5,10, 25,50, 100,250, 700, 2000, 5000), trans="log1p", expand=c(0,0)) +
  scale_y_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), limits=c(0,0.5), expand=c(0,0)) +
  scale_fill_manual(values = myColors) +
  facet_grid(rows = vars(As_regulated)) +
  theme(strip.text = element_text(face="bold", size=12,lineheight=5.0),
        strip.background = element_rect(fill="gold", colour="black",size=1),
        legend.position = "none",
        axis.title.y = element_text(face = 'bold', size = 14),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank())

##############################  
  
Cortex_FPKM_hist <- bind_rows(Cortex_up_FPKM, Cortex_down_FPKM) %>%
  dplyr::rename('0As'= Cortex_0As_FPKM, '50As'= Cortex_50As_FPKM) %>% 
  pivot_longer(
    cols = contains("0As"),
    names_to = "Condition",
    values_to = "FPKM",
    values_drop_na = TRUE)

Cortex_FPKM_hist <- Cortex_FPKM_hist %>%     
  mutate(As_regulated = ifelse(log2FoldChange > 0, "Up-regregulated", "Down-regulated"))

Cortex_den_FPKM <- ggplot(Cortex_FPKM_hist, aes(x= FPKM)) + 
  geom_density(aes( fill = Condition),alpha=0.5) +
  labs(subtitle =expression(bolditalic("Cortex Cells"))) +
  scale_x_continuous(breaks=c(0,1,2,5,10, 25,50, 100,250, 700, 2000, 5000), trans="log1p", expand=c(0,0)) +
  scale_y_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), limits=c(0,0.5), expand=c(0,0)) +
  scale_fill_manual(values = myColors) +
  facet_grid(rows = vars(As_regulated)) +
  theme(strip.text = element_text(face="bold", size=12,lineheight=5.0),
        strip.background = element_rect(fill="gold", colour="black",size=1),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(face='bold', size=14),
        panel.grid.minor = element_blank())

#########################################################
SCR_FPKM_hist <- bind_rows(SCR_up_FPKM, SCR_down_FPKM) %>%
  dplyr::rename('0As'= SCR_0As_FPKM, '50As'= SCR_50As_FPKM) %>% 
  pivot_longer(
    cols = contains("0As"),
    names_to = "Condition",
    values_to = "FPKM",
    values_drop_na = TRUE)

SCR_FPKM_hist <- SCR_FPKM_hist %>%     
  mutate(As_regulated = ifelse(log2FoldChange > 0, "Up-regregulated", "Down-regulated"))

SCR_den_FPKM <- ggplot(SCR_FPKM_hist, aes(x= FPKM)) + 
  geom_density(aes( fill = Condition),alpha=0.6, color='black', weight=5) +
  labs(subtitle =expression(bolditalic("Endodermal Cells"))) +
  scale_x_continuous(breaks=c(0,1,2,5,10, 25,50, 100,250, 700, 2000, 5000), trans="log1p", expand=c(0,0)) +
  scale_y_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), limits=c(0,0.5), expand=c(0,0)) +
  scale_fill_manual(values = myColors) +
  facet_grid(rows = vars(As_regulated)) +
  theme(strip.text = element_text(face="bold", size=12,lineheight=5.0),
        strip.background = element_rect(fill='gold', colour="black",size=1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank())
WER_den_FPKM+Cortex_den_FPKM+SCR_den_FPKM
# ggsave('FPKM_density.tiff', plot=last_plot())
#################################################################################

