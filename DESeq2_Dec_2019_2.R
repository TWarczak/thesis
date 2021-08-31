library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(patchwork)
library(BiocManager)
library(biomaRt)
library(grid)
library(venn)
library(ggpolypath)
library(viridis)
#BiocManager::install("org.At.tair.db")
library(dendextend)  # use this to determine which dendogram method to use for clustering rows and how many clusters you want
library(org.At.tair.db)
library(DOSE)
library(clusterProfiler)
library(gridExtra)


#FACS_Todd_counts<- read_delim("C:/Users/Todd/OneDrive for Business/Todd_RNAseq/July_2018_Data/sharing_Data_RNA_Guerinot/sharing_Data_RNA_Guerinot2/RNAseq_WER_Cortex_SCR_genes_count_July2018.txt", delim = "\t")
FACS_Todd_counts<- read_delim("~/R/DESeq2_Todd_RNAseq/RNAseq_WER_Cortex_SCR_genes_count_July2018.txt", delim = "\t")
# if you want to save raw counts as its own file, execute the code below, but with new path
#FACS_Todd_counts_tbl <- as_tibble(FACS_Todd_counts, rownames = NA) # default rownames are removed with as_tibble, need rownames = NA
#write.csv(FACS_Todd_counts_tbl, "FACS_RNAseq_Todd_counts_tbl.csv", col.names = TRUE, row.names = TRUE)  # can't use write_csv, it removes rownames

#glimpse(FACS_Todd_counts)

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

#class(FACS_Todd_counts) #make sure its a count matrix

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
#norm_counts_tbl <- as_tibble(norm_counts_dds_FACS, rownames = NA) # default rownames are removed with as_tibble, need rownames = NA
#write.csv(norm_counts_tbl, "FACS_RNAseq_norm_counts_tbl.csv", col.names = TRUE, row.names = TRUE)  # can't use write_csv, it removes rownames

#class(norm_counts_dds_FACS)
library(EDASeq)
#par(mfrow = c(1, 2))
#plotRLE(FACS_Todd_counts, outline=FALSE, ylim=c(-4, 4), col=as.numeric(colData$group), main = 'Raw Counts')
#plotRLE(norm_counts_dds_FACS, outline=FALSE, ylim=c(-4, 4), col = as.numeric(colData$group), main = 'Normalized Counts')

########################################################################
# Hierarchical Clustering with correlation heatmaps
vsd_FACS <- vst(dds_FACS, blind = TRUE)

vsd_mat <- assay(vsd_FACS)
vsd_cor <- cor(vsd_mat)
#pheatmap <- pheatmap(vsd_cor, annotation = dplyr::select(FACS_metadata, condition))
# ggsave("hier_clus_cor_heatmap.tiff", plot = pheatmap, device = "tiff", dpi = 'retina')

#######################################################################
# PCA Analysis

# Transform the normalized counts   
# Create the PCA plot for PC1 and PC2 and color by condition

# pcaData <- plotPCA(vsd_FACS, intgroup=c("condition", "celltype"), returnData=TRUE)

#percentVar <- round(100 * attr(pcaData, "percentVar"))

# ggplot(pcaData, aes(PC1, PC2, color=celltype, shape=condition)) +
#   geom_point(size=4) +
#   scale_shape_manual(values=c(15, 17)) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#   coord_fixed()

# ggsave("PCA.tiff", plot = last_plot(), device = "tiff", dpi = 'retina')

#########################################################################
# Need to make groups to use in comparisons. i.e. the 3 reps of WER-As become the group: WER0_As
dds_FACS$group <- factor(paste0(dds_FACS$celltype, dds_FACS$condition))
design(dds_FACS) <- ~ group

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
# plotDispEsts(dds_FACS_DESeq)
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

# summary(results_WER)
# summary(results_Cortex)
# summary(results_SCR)

# To generate more accurate log2 foldchange estimates, use Shrunken log2 foldchanges (LFC)
# for rationale go to:
# https://hbctraining.github.io/DGE_workshop_salmon/lessons/05_DGE_DESeq2_analysis2.html
# produces the same DESeqResults object, but with the log2FoldChange and lfcSE columns replaced with shrunken (more accurate) results.  baseMean, stat, pvalue, & padj unchanged

results_WER_shr <- lfcShrink(dds_FACS_DESeq, contrast=c("group","WER50_As","WER0_As"), res=results_WER)
results_SCR_shr <- lfcShrink(dds_FACS_DESeq, contrast=c("group","SCR50_As","SCR0_As"), res=results_SCR)
results_Cortex_shr <- lfcShrink(dds_FACS_DESeq, contrast=c("group","Cortex50_As","Cortex0_As"), res=results_Cortex)

# summary(results_WER_shr)

# plotMA(results_WER_shr, alpha = 0.05, ylim = c(-10,15), xlab="",ylab="", frame.plot = F)   
# mytitle = "MA plot, Epidermal Genes" 
# mysubtitle = "0 AsIII -> 50 AsIII" 
# mtext(side=3, line=-1, at=0.07, font= 2, adj=0, cex=1.3, mytitle) 
# mtext(side=3, line=-2, at=0.07, adj=0, cex=0.9, mysubtitle) 
# mtext(side=1, line=2, "Mean of Normalized Counts", font=2,cex=1.1) 
# mtext(side=2, line=3, "Log2 Fold Change", font=2, cex=1.1) 
# # save from plot window
# plotMA(results_SCR_shr, alpha = 0.05, ylim = c(-10,15), xlab="",ylab="", frame.plot = F) 
# mytitle = "MA plot, Endodermal Genes"
# mysubtitle = "0 AsIII -> 50 AsIII"
# mtext(side=3, line=-1, at=0.07, font= 2, adj=0, cex=1.3, mytitle)
# mtext(side=3, line=-2, at=0.07, adj=0, cex=0.9, mysubtitle)
# mtext(side=1, line=2, "Mean of Normalized Counts", font=2,cex=1.1)
# mtext(side=2, line=3, "Log2 Fold Change", font=2, cex=1.1)
# # save from plot window
# plotMA(results_Cortex_shr, alpha = 0.05, ylim = c(-10,15), xlab="",ylab="", frame.plot = F) 
# mytitle = "MA plot, Cortex Genes"
# mysubtitle = "0 AsIII -> 50 AsIII"
# mtext(side=3, line=-1, at=0.07, font= 2, adj=0, cex=1.3, mytitle)
# mtext(side=3, line=-2, at=0.07, adj=0, cex=0.9, mysubtitle)
# mtext(side=1, line=2, "Mean of Normalized Counts", font=2,cex=1.1)
# mtext(side=2, line=3, "Log2 Fold Change", font=2, cex=1.1)
# # save from plot window

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
# p_W <- ggplot(results_WER2) +
#   geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold), alpha = 0.5, size = 1.5) +
#   theme_minimal() +
#   labs(x="", y ="-log10 adjusted p-value", title = "Epidermal Genes", subtitle =  "0 AsIII -> 50 AsIII") +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.2), face = 'bold'),
#         axis.title = element_text(size = rel(1), face = 'bold'))
# p_S <- ggplot(results_SCR2) +
#   geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold), alpha = 0.5, size = 1.5) +
#   theme_minimal() +
#   labs(x="", y ="", title = "Endodermal Genes", subtitle =  "0 AsIII -> 50 AsIII") +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.2), face = 'bold'),
#         axis.title = element_text(size = rel(1), face = 'bold'))
# p_C <- ggplot(results_Cortex2) +
#   geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold), alpha = 0.5, size = 1.5) +
#   theme_minimal() +
#   labs(x="log2 fold change", y ="", title = "Cortex Genes", subtitle =  "0 AsIII -> 50 AsIII") +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.2), face = 'bold'),
#         axis.title = element_text(size = rel(1), face = 'bold'))
# 
# p_W + p_C + p_S
#ggsave("volcano_all.tiff", plot = last_plot(), device = "tiff", dpi = 'retina')

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
RPKM_complete <- read_csv("~/R/DESeq2_Todd_RNAseq/FPKM_complete.csv")
write_csv(RPKM_complete, "rpkm_complete.csv", col_names = TRUE)

RPKM_avg <- RPKM_complete %>% 
  transmute(gene_id = gene_id,
            WER_minus_As = (RPKM_complete$Wer_minus_AS_0 + RPKM_complete$Wer_minus_AS_1 + RPKM_complete$Wer_minus_AS_2)/3,
            WER_plus_As = (RPKM_complete$Wer_plus_AS_0 + RPKM_complete$Wer_plus_AS_1 + RPKM_complete$Wer_plus_AS_2)/3,
            Cortex_minus_As = (RPKM_complete$Cortex_minus_AS_0 + RPKM_complete$Cortex_minus_AS_1 + RPKM_complete$Cortex_minus_AS_2)/3, 
            Cortex_plus_As = (RPKM_complete$Cortex_plus_AS_0 + RPKM_complete$Cortex_plus_AS_1 + RPKM_complete$Cortex_plus_AS_2)/3,  
            SCR_minus_As = (RPKM_complete$SCR_minus_AS_0 + RPKM_complete$SCR_minus_AS_1 + RPKM_complete$SCR_minus_AS_2)/3, 
            SCR_plus_As = (RPKM_complete$SCR_plus_AS_0 + RPKM_complete$SCR_plus_AS_1 + RPKM_complete$SCR_plus_AS_2)/3)
# round and put in order for heatmap or table
RPKM_avg_round <- RPKM_avg %>% 
  transmute(gene_id = gene_id,
            WER_minus_As = round(WER_minus_As, digits = 0),
            Cortex_minus_As = round(Cortex_minus_As, digits = 0),
            SCR_minus_As = round(SCR_minus_As, digits = 0),
            WER_plus_As = round(WER_plus_As, digits = 0),
            Cortex_plus_As = round(Cortex_plus_As, digits = 0),
            SCR_plus_As = round(SCR_plus_As, digits = 0)) 
View(RPKM_avg_round)
write_csv(RPKM_avg_round, "rpkm_complete_avg_round.csv", col_names = TRUE)

#######################
WER_up_RPKM <- sig_WER %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange > 1) %>% 
  inner_join(RPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, WER_minus_As, WER_plus_As) %>% 
  filter(WER_plus_As > 6) %>% 
  dplyr::rename(WER_0As_RPKM = WER_minus_As, WER_50As_RPKM = WER_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   WER_0As_RPKM = round(WER_0As_RPKM, digits = 2), 
                   WER_50As_RPKM = round(WER_50As_RPKM, digits = 2))

WER_down_RPKM <- sig_WER %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange < -1) %>% 
  inner_join(RPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, WER_minus_As, WER_plus_As) %>% 
  filter(WER_minus_As > 6) %>% 
  dplyr::rename(WER_0As_RPKM = WER_minus_As, WER_50As_RPKM = WER_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   WER_0As_RPKM = round(WER_0As_RPKM, digits = 2), 
                   WER_50As_RPKM = round(WER_50As_RPKM, digits = 2))

SCR_up_RPKM <- sig_SCR %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange > 1) %>% 
  inner_join(RPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, SCR_minus_As, SCR_plus_As) %>% 
  filter(SCR_plus_As > 6) %>% 
  dplyr::rename(SCR_0As_RPKM = SCR_minus_As, SCR_50As_RPKM = SCR_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   SCR_0As_RPKM = round(SCR_0As_RPKM, digits = 2), 
                   SCR_50As_RPKM = round(SCR_50As_RPKM, digits = 2))

SCR_down_RPKM <- sig_SCR %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange < -1) %>% 
  inner_join(RPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, SCR_minus_As, SCR_plus_As) %>% 
  filter(SCR_minus_As > 6) %>% 
  dplyr::rename(SCR_0As_RPKM = SCR_minus_As, SCR_50As_RPKM = SCR_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   SCR_0As_RPKM = round(SCR_0As_RPKM, digits = 2), 
                   SCR_50As_RPKM = round(SCR_50As_RPKM, digits = 2))

Cortex_up_RPKM <- sig_Cortex %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange > 1) %>% 
  inner_join(RPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, Cortex_minus_As, Cortex_plus_As) %>% 
  filter(Cortex_plus_As > 6) %>% 
  dplyr::rename(Cortex_0As_RPKM = Cortex_minus_As, Cortex_50As_RPKM = Cortex_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   Cortex_0As_RPKM = round(Cortex_0As_RPKM, digits = 2), 
                   Cortex_50As_RPKM = round(Cortex_50As_RPKM, digits = 2))

Cortex_down_RPKM <- sig_Cortex %>% 
  dplyr::rename(gene_id = gene) %>%   
  filter(log2FoldChange < -1) %>% 
  inner_join(RPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, log2FoldChange, Cortex_minus_As, Cortex_plus_As) %>% 
  filter(Cortex_minus_As > 6) %>% 
  dplyr::rename(Cortex_0As_RPKM = Cortex_minus_As, Cortex_50As_RPKM = Cortex_plus_As ) %>% 
  dplyr::transmute(gene_id = gene_id,
                   log2FoldChange = round(log2FoldChange, digits = 2), 
                   Cortex_0As_RPKM = round(Cortex_0As_RPKM, digits = 2), 
                   Cortex_50As_RPKM = round(Cortex_50As_RPKM, digits = 2))


# venn = 
#    list(WER_up_RPKM$gene_id, Cortex_up_RPKM$gene_id, SCR_up_RPKM$gene_id, WER_down_RPKM$gene_id, Cortex_down_RPKM$gene_id, SCR_down_RPKM$gene_id)
#      
# tiff('ven_march12.tiff', width = 3700, height = 3700,res = 800)
# 
# venn.result =
#    venn(venn, ilabels = T, 
#         zcolor = "style", size = 85, cexil = 5, cexsn = 2.3, box=F, ilcs=0.55, sncs = 0.6, opacity = 0.22,
#         snames = c('Epiderm_Up', 'Cortex_Up', 'Endoderm_Up', 'Epiderm_Down', 'Cortex_Down', 'Endoderm_Down'));
#  
#  dev.off()
# 
# View(venn.result)
# filt_venn.result <- venn.result %>% 
#   rownames_to_column('group') %>% 
#   filter(counts>0)
# 

End_down_841 <- SCR_down_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(End_down_841, 'End_down_841.csv', col_names = T)
# 
Cor_down_314 <- Cortex_down_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Cor_down_314, 'Cor_down_314.csv', col_names = T)
# 
Cor_down_End_down_464 <- Cortex_down_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  semi_join(SCR_down_RPKM, by ='gene_id') %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Cor_down_End_down_464, 'Cor_down_End_down_464.csv', col_names = T)
# 
Epi_down_885 <- WER_down_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Epi_down_885, 'Epi_down_885.csv', col_names = T)
# 
Epi_down_End_down_277 <- WER_down_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  semi_join(SCR_down_RPKM, by ='gene_id') %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Epi_down_End_down_277, 'Epi_down_End_down_277.csv', col_names = T)
# 
Epi_down_Cor_down_472 <- WER_down_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  semi_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Epi_down_Cor_down_472, 'Epi_down_Cor_down_472.csv', col_names = T)
# 
Epi_down_Cor_down_End_down_1676 <- WER_down_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  semi_join(Cortex_down_RPKM, by ='gene_id') %>%
  semi_join(SCR_down_RPKM, by ='gene_id') %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Epi_down_Cor_down_End_down_1676, 'Epi_down_Cor_down_End_down_1676.csv', col_names = T)

End_up_306 <- SCR_up_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(WER_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(End_up_306, 'End_up_306.csv', col_names = T)
# 
End_up_Epi_down_10 <- SCR_up_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(WER_up_RPKM, by = 'gene_id') %>%
  semi_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(End_up_Epi_down_10, 'End_up_Epi_down_10.csv', col_names = T)
# 
Cor_up_247 <- Cortex_up_RPKM %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Cor_up_247, 'Cor_up_247.csv', col_names = T)
# 
Cor_up_Epi_down_4 <- Cortex_up_RPKM %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  semi_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Cor_up_Epi_down_4, 'Cor_up_Epi_down_4.csv', col_names = T)
# 
Cor_up_End_down_4 <- Cortex_up_RPKM %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  semi_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Cor_up_End_down_4, 'Cor_up_End_down_4.csv', col_names = T)
# 
Cor_up_End_up_475 <- Cortex_up_RPKM %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  semi_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Cor_up_End_up_475, 'Cor_up_End_up_475.csv', col_names = T)
# 
Cor_up_End_up_Epi_down_2 <- Cortex_up_RPKM %>%
  anti_join(WER_up_RPKM, by ='gene_id') %>%
  semi_join(SCR_up_RPKM, by = 'gene_id') %>%
  semi_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Cor_up_End_up_Epi_down_2, 'Cor_up_End_up_Epi_down_2.csv', col_names = T)
# 
Epi_up_669 <- WER_up_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Epi_up_669, 'Epi_up_669.csv', col_names = T)
# 
Epi_up_End_down_59 <- WER_up_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  semi_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Epi_up_End_down_59, 'Epi_up_End_down_59.csv', col_names = T)
# 
Epi_up_Cor_down_End_down_7 <- WER_up_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  semi_join(Cortex_down_RPKM, by ='gene_id') %>%
  semi_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Epi_up_Cor_down_End_down_7, 'Epi_up_Cor_down_End_down_7.csv', col_names = T)
# 
Epi_up_End_up_100 <- WER_up_RPKM %>%
  anti_join(Cortex_up_RPKM, by ='gene_id') %>%
  semi_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Epi_up_End_up_100, 'Epi_up_End_up_100.csv', col_names = T)
# 
Epi_up_Cor_up_520 <- WER_up_RPKM %>%
  semi_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Epi_up_Cor_up_520, 'Epi_up_Cor_up_520.csv', col_names = T)
# 
Epi_up_Cor_up_End_down_5 <- WER_up_RPKM %>%
  semi_join(Cortex_up_RPKM, by ='gene_id') %>%
  anti_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  semi_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Epi_up_Cor_up_End_down_5, 'Epi_up_Cor_up_End_down_5.csv', col_names = T)
# 
Epi_up_Cor_up_End_up_1859 <- WER_up_RPKM %>%
  semi_join(Cortex_up_RPKM, by ='gene_id') %>%
  semi_join(SCR_up_RPKM, by = 'gene_id') %>%
  anti_join(WER_down_RPKM, by ='gene_id') %>%
  anti_join(Cortex_down_RPKM, by ='gene_id') %>%
  anti_join(SCR_down_RPKM, by ='gene_id') %>%
  dplyr::select(gene_id) %>% 
  left_join(RPKM_avg_round, by = 'gene_id')
write_csv(Epi_up_Cor_up_End_up_1859, 'Epi_up_Cor_up_End_up_1859.csv', col_names = T)

###################################################

# write_delim(WER_up_RPKM$gene_id, '~\WER_up_geneID.csv', 
#             delim = " ", na = "NA", append = FALSE,col_names = F, quote_escape = "double")
# 
# path = '~/R/DESeq2_Todd_RNAseq'
# write.csv(WER_up_RPKM,'WER_up_geneID.csv',sep = '')

##########################
# RPKM_density plot, faceted by up/down regulation & cell-type

myColors <- c('#4E84C4', 'tomato')
my_title <-"RPKM Values of Differentially Expressed Genes"

WER_RPKM_hist <- bind_rows(WER_up_RPKM, WER_down_RPKM) %>%
  dplyr::rename('0As'= WER_0As_RPKM, '50As'= WER_50As_RPKM) %>%
  pivot_longer(
    cols = contains("0As"),
    names_to = "Condition",
    values_to = "RPKM",
    values_drop_na = TRUE)

WER_RPKM_hist <- WER_RPKM_hist %>%
  mutate(As_regulated = ifelse(log2FoldChange > 0, "Up-regulated", "Down-regulated"))

# reorder levels so up-regulated genes are top plot
WER_RPKM_hist$As_regulated = factor(WER_RPKM_hist$As_regulated, levels=c("Up-regulated","Down-regulated"))

WER_den_RPKM <- ggplot(WER_RPKM_hist, aes(x= RPKM)) +
  geom_density(aes( fill = Condition),alpha=0.5) +
  labs(y='Density',
       #title =my_title,
       subtitle =expression(bolditalic("Epidermal Cells"))) +
  scale_x_continuous(breaks=c(0,1,3,6,10, 25,50, 100,250, 700, 2000, 5000),
                     trans="log1p", expand=c(0,0),
                     guide = guide_axis(angle=45)) +
  scale_y_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), limits=c(0,0.46), expand=c(0,0)) +
  scale_fill_manual(values = myColors) +
  facet_grid(rows = vars(As_regulated)) +
  theme(strip.text = element_text(face="bold", size=12,lineheight=5.0),
        strip.background = element_rect(fill="gold", colour="black",size=1),
        legend.position = "none",
        axis.title.y = element_text(face = 'bold', size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        panel.grid.minor = element_blank())
#
# ##############################
#
Cortex_RPKM_hist <- bind_rows(Cortex_up_RPKM, Cortex_down_RPKM) %>%
  dplyr::rename('0As'= Cortex_0As_RPKM, '50As'= Cortex_50As_RPKM) %>%
  pivot_longer(
    cols = contains("0As"),
    names_to = "Condition",
    values_to = "RPKM",
    values_drop_na = TRUE)

Cortex_RPKM_hist <- Cortex_RPKM_hist %>%
  mutate(As_regulated = ifelse(log2FoldChange > 0, "Up-regulated", "Down-regulated"))

# reorder levels so up-regulated genes are top plot
Cortex_RPKM_hist$As_regulated = factor(Cortex_RPKM_hist$As_regulated, levels=c("Up-regulated","Down-regulated"))

Cortex_den_RPKM <- ggplot(Cortex_RPKM_hist, aes(x= RPKM)) +
  geom_density(aes( fill = Condition),alpha=0.5) +
  labs(subtitle =expression(bolditalic("Cortex Cells"))) +
  scale_x_continuous(breaks=c(0,1,3,6,10, 25,50, 100,250, 700, 2000, 5000),
                     trans="log1p", expand=c(0,0),
                     guide = guide_axis(angle=45)) +
  scale_y_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), limits=c(0,0.46), expand=c(0,0)) +
  scale_fill_manual(values = myColors) +
  facet_grid(rows = vars(As_regulated)) +
  theme(strip.text = element_text(face="bold", size=12,lineheight=5.0),
        strip.background = element_rect(fill="gold", colour="black",size=1),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(face='bold', size=14),
        axis.text.x = element_text(size=13),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank())

# #########################################################
SCR_RPKM_hist <- bind_rows(SCR_up_RPKM, SCR_down_RPKM) %>%
  dplyr::rename('0As'= SCR_0As_RPKM, '50As'= SCR_50As_RPKM) %>%
  pivot_longer(
    cols = contains("0As"),
    names_to = "Condition",
    values_to = "RPKM",
    values_drop_na = TRUE)

SCR_RPKM_hist <- SCR_RPKM_hist %>%
  mutate(As_regulated = ifelse(log2FoldChange > 0, "Up-regulated", "Down-regulated"))

# reorder levels so up-regulated genes are top plot
SCR_RPKM_hist$As_regulated = factor(SCR_RPKM_hist$As_regulated, levels=c("Up-regulated","Down-regulated"))

SCR_den_RPKM <- ggplot(SCR_RPKM_hist, aes(x= RPKM)) +
  geom_density(aes( fill = Condition),alpha=0.6, color='black', weight=5) +
  labs(subtitle =expression(bolditalic("Endodermal Cells"))) +
  scale_x_continuous(breaks=c(0,1,3,6,10, 25,50, 100,250, 700, 2000, 5000),
                     trans="log1p", expand=c(0,0),
                     guide = guide_axis(angle=45)) +
  scale_y_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), limits=c(0,0.46), expand=c(0,0)) +
  scale_fill_manual(values = myColors) +
  facet_grid(rows = vars((As_regulated))) +
  theme(strip.text = element_text(face="bold", size=12,lineheight=5.0),
        strip.background = element_rect(fill='gold', colour="black",size=1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=13),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank())
WER_den_RPKM+Cortex_den_RPKM+SCR_den_RPKM
ggsave('RPKM_density.tiff', plot=last_plot(), dpi = 650)
# 2 supplementary tables for density plots. 1 table for up DEGs, 1 table for down DEGs
# 1st variable if RPKM values, binned 0-6, 6-10, 10-25, 25-50, 50-100, 100-250, 250-700, 700-2000, 2000-5000, 5000+
# Variables 2-4 are Epi, Cor, Endo for control condition, Variables 5-7 are Epi, Cor, Endo for AsIII condition

install.packages('rbin')
library(rbin)
bins <- rbin_manual(WER_up_RPKM, gene_id, WER_0As_RPKM, cut_points=c(6, 10, 25, 50, 100, 250, 700, 2000, 5000))
bins
?rbin
View(mbank)
?rbin_manual

b <- c(-Inf, 10000, 31000, Inf)
names <- c("Low", "Medium", "High")
students$Income.cat <- cut(students$Income, breaks = b, labels = names)

bins <- c(-Inf, 6, 10, 25, 50, 100, 250, 700, 2000, 5000, Inf)
RPKM_Range <- c("0-6", "6-10", "10-25", "25-50", "50-100", "100-250","250-700", "700-2000", "2000-5000", "5000+" )

WER_up_RPKM$RPKM_Range_control <- cut(WER_up_RPKM$WER_0As_RPKM, breaks = bins, labels = RPKM_Range)
Epid_up_Control <- dplyr::count(WER_up_RPKM, RPKM_Range_control)

WER_up_RPKM$RPKM_Range_AsIII <- cut(WER_up_RPKM$WER_50As_RPKM, breaks = bins, labels = RPKM_Range)
Epid_up_AsIII <- dplyr::count(WER_up_RPKM, RPKM_Range_AsIII)

Cortex_up_RPKM$RPKM_Range_control <- cut(Cortex_up_RPKM$Cortex_0As_RPKM, breaks = bins, labels = RPKM_Range)
Cort_up_Control <- dplyr::count(Cortex_up_RPKM, RPKM_Range_control)

Cortex_up_RPKM$RPKM_Range_AsIII <- cut(Cortex_up_RPKM$Cortex_50As_RPKM, breaks = bins, labels = RPKM_Range)
Cort_up_AsIII <- dplyr::count(Cortex_up_RPKM, RPKM_Range_AsIII)

SCR_up_RPKM$RPKM_Range_control <- cut(SCR_up_RPKM$SCR_0As_RPKM, breaks = bins, labels = RPKM_Range)
Endo_up_Control <- dplyr::count(SCR_up_RPKM, RPKM_Range_control)

SCR_up_RPKM$RPKM_Range_AsIII <- cut(SCR_up_RPKM$SCR_50As_RPKM, breaks = bins, labels = RPKM_Range)
Endo_up_AsIII <- dplyr::count(SCR_up_RPKM, RPKM_Range_AsIII)
cumsum(Endo_up_AsIII$n)

Sup_Table3.1 <- Epid_up_Control %>% 
  left_join(Cort_up_Control, by= 'RPKM_Range_control') %>% 
  left_join(Endo_up_Control, by= 'RPKM_Range_control') %>% 
  dplyr::rename('RPKM_Range' = RPKM_Range_control, 'Epid_Control' = n.x, 'Cort_Control' = n.y, 'Endo_Control' = n)

Sup_Table3.2 <- Endo_up_AsIII %>% 
  full_join(Cort_up_AsIII, by= 'RPKM_Range_AsIII') %>% 
  full_join(Epid_up_AsIII, by= 'RPKM_Range_AsIII') %>% 
  dplyr::rename('RPKM_Range' = RPKM_Range_AsIII,'Endo_AsIII' = n.x, 'Cort_AsIII' = n.y, 'Epid_AsIII' = n)

Sup_Table3 <- Sup_Table3.1 %>% 
  full_join(Sup_Table3.2, by = 'RPKM_Range') %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(RPKM_Range, Epid_Control, Cort_Control, Endo_Control, Epid_AsIII, Cort_AsIII, Endo_AsIII)

write_csv(Sup_Table3, "Sup_Table3.csv", col_names = T)
# 
WER_down_RPKM$RPKM_Range_control <- cut(WER_down_RPKM$WER_0As_RPKM, breaks = bins, labels = RPKM_Range)
Epid_down_Control <- dplyr::count(WER_down_RPKM, RPKM_Range_control)

WER_down_RPKM$RPKM_Range_AsIII <- cut(WER_down_RPKM$WER_50As_RPKM, breaks = bins, labels = RPKM_Range)
Epid_down_AsIII <- dplyr::count(WER_down_RPKM, RPKM_Range_AsIII)

Cortex_down_RPKM$RPKM_Range_control <- cut(Cortex_down_RPKM$Cortex_0As_RPKM, breaks = bins, labels = RPKM_Range)
Cort_down_Control <- dplyr::count(Cortex_down_RPKM, RPKM_Range_control)

Cortex_down_RPKM$RPKM_Range_AsIII <- cut(Cortex_down_RPKM$Cortex_50As_RPKM, breaks = bins, labels = RPKM_Range)
Cort_down_AsIII <- dplyr::count(Cortex_down_RPKM, RPKM_Range_AsIII)

SCR_down_RPKM$RPKM_Range_control <- cut(SCR_down_RPKM$SCR_0As_RPKM, breaks = bins, labels = RPKM_Range)
Endo_down_Control <- dplyr::count(SCR_down_RPKM, RPKM_Range_control)

SCR_down_RPKM$RPKM_Range_AsIII <- cut(SCR_down_RPKM$SCR_50As_RPKM, breaks = bins, labels = RPKM_Range)
Endo_down_AsIII <- dplyr::count(SCR_down_RPKM, RPKM_Range_AsIII)

Sup_Table4.1 <- Epid_down_Control %>% 
  left_join(Cort_down_Control, by= 'RPKM_Range_control') %>% 
  left_join(Endo_down_Control, by= 'RPKM_Range_control') %>% 
  dplyr::rename('RPKM_Range' = RPKM_Range_control, 'Epid_Control' = n.x, 'Cort_Control' = n.y, 'Endo_Control' = n)

Sup_Table4.2 <- Epid_down_AsIII %>% 
  left_join(Cort_down_AsIII, by= 'RPKM_Range_AsIII') %>% 
  left_join(Endo_down_AsIII, by= 'RPKM_Range_AsIII') %>% 
  dplyr::rename('RPKM_Range' = RPKM_Range_AsIII,'Epid_AsIII' = n.x, 'Cort_AsIII' = n.y, 'Endo_AsIII' = n)

Sup_Table4 <- Sup_Table4.1 %>% 
  full_join(Sup_Table4.2, by = 'RPKM_Range') %>% 
  replace(is.na(.), 0)

write_csv(Sup_Table4, "Sup_Table4.csv", col_names = T)

# Took these tables to excel to customize 
# Final product in manuscript_figures



###### GO term enrichment #######################################################################################

 # might need to click from Packages tab 

columns(org.At.tair.db)
#keys <- head( keys(org.At.tair.db) )

#### WER
sig_WER <- sig_WER %>% 
  dplyr::rename(gene_id = gene)

id_list <- sig_WER$gene_id

symbol_list <- mapIds(org.At.tair.db, keys= id_list, column=c("SYMBOL"), keytype="TAIR", multiVals="first")
ENTREZID <- mapIds(org.At.tair.db, keys= id_list, column=c("ENTREZID"), keytype="TAIR", multiVals="first")
ENTREZID <- as.data.frame(ENTREZID)
description_list <- mapIds(org.At.tair.db, keys= id_list, column=c("GENENAME"), keytype="TAIR", multiVals="first")
descriptions <- as.data.frame(description_list)
sig_WER_info <- cbind(sig_WER, symbol_list, descriptions,ENTREZID)

sig_WER_all_info <- sig_WER_info %>% 
  inner_join(RPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj, 
                symbol_list, description_list, ENTREZID, WER_minus_As, WER_plus_As) %>% 
  dplyr::transmute(gene_id = gene_id,
                   baseMean = round(baseMean, digits = 2),
                   log2FoldChange = round(log2FoldChange, digits = 2),
                   lfcSE = round(lfcSE, digits = 3),
                   stat = round(stat, digits = 3),
                   pvalue = pvalue,
                   padj = padj,
                   symbol_list = symbol_list,
                   description_list = description_list,
                   ENTREZID = ENTREZID,
                   WER_0As_RPKM = round(WER_minus_As, digits = 2), 
                   WER_50As_RPKM = round(WER_plus_As, digits = 2))

as.data.frame(sig_WER_all_info) %>%
  write_csv('sig_WER_all_info.csv',col_names = T) # Saving all genes differentially expressed in WER

#### Cortex 
sig_Cortex <- sig_Cortex %>% 
  dplyr::rename(gene_id = gene)

id_list <- sig_Cortex$gene_id

symbol_list <- mapIds(org.At.tair.db, keys= id_list, column=c("SYMBOL"), keytype="TAIR", multiVals="first")
ENTREZID <- mapIds(org.At.tair.db, keys= id_list, column=c("ENTREZID"), keytype="TAIR", multiVals="first")
ENTREZID <- as.data.frame(ENTREZID)
description_list <- mapIds(org.At.tair.db, keys= id_list, column=c("GENENAME"), keytype="TAIR", multiVals="first")
descriptions <- as.data.frame(description_list)
sig_Cortex_info <- cbind(sig_Cortex, symbol_list, descriptions,ENTREZID)

sig_Cortex_all_info <- sig_Cortex_info %>% 
  inner_join(RPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj, 
                symbol_list, description_list, ENTREZID, Cortex_minus_As, Cortex_plus_As) %>% 
  dplyr::transmute(gene_id = gene_id,
                   baseMean = round(baseMean, digits = 2),
                   log2FoldChange = round(log2FoldChange, digits = 2),
                   lfcSE = round(lfcSE, digits = 3),
                   stat = round(stat, digits = 3),
                   pvalue = pvalue,
                   padj = padj,
                   symbol_list = symbol_list,
                   description_list = description_list,
                   ENTREZID = ENTREZID,
                   Cortex_0As_RPKM = round(Cortex_minus_As, digits = 2), 
                   Cortex_50As_RPKM = round(Cortex_plus_As, digits = 2))

as.data.frame(sig_Cortex_all_info) %>%
  write_csv('sig_Cortex_all_info.csv',col_names = T) # Saving all genes differentially expressed in Cortex

#########SCR
sig_SCR <- sig_SCR %>% 
  dplyr::rename(gene_id = gene)

id_list <- sig_SCR$gene_id

symbol_list <- mapIds(org.At.tair.db, keys= id_list, column=c("SYMBOL"), keytype="TAIR", multiVals="first")
ENTREZID <- mapIds(org.At.tair.db, keys= id_list, column=c("ENTREZID"), keytype="TAIR", multiVals="first")
ENTREZID <- as.data.frame(ENTREZID)
description_list <- mapIds(org.At.tair.db, keys= id_list, column=c("GENENAME"), keytype="TAIR", multiVals="first")
descriptions <- as.data.frame(description_list)
sig_SCR_info <- cbind(sig_SCR, symbol_list, descriptions,ENTREZID)

sig_SCR_all_info <- sig_SCR_info %>% 
  inner_join(RPKM_avg, by = 'gene_id') %>% 
  dplyr::select(gene_id, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj, 
                symbol_list, description_list, ENTREZID, SCR_minus_As, SCR_plus_As) %>% 
  dplyr::transmute(gene_id = gene_id,
                   baseMean = round(baseMean, digits = 2),
                   log2FoldChange = round(log2FoldChange, digits = 2),
                   lfcSE = round(lfcSE, digits = 3),
                   stat = round(stat, digits = 3),
                   pvalue = pvalue,
                   padj = padj,
                   symbol_list = symbol_list,
                   description_list = description_list,
                   ENTREZID = ENTREZID,
                   SCR_0As_RPKM = round(SCR_minus_As, digits = 2), 
                   SCR_50As_RPKM = round(SCR_plus_As, digits = 2))

# as.data.frame(sig_SCR_all_info) %>%
#   write_csv('sig_SCR_all_info.csv',col_names = T) # Saving all genes differentially expressed in SCR

###############################################
sig_WER_all_info_filter <- sig_WER_all_info %>% 
  filter(abs(log2FoldChange) > 1, (WER_0As_RPKM > 6 | WER_50As_RPKM > 6))
sig_Cortex_all_info_filter <- sig_Cortex_all_info %>% 
  filter(abs(log2FoldChange) > 1, (Cortex_0As_RPKM > 6 | Cortex_50As_RPKM > 6))
sig_SCR_all_info_filter <- sig_SCR_all_info %>% 
  filter(abs(log2FoldChange) > 1, (SCR_0As_RPKM > 6 | SCR_50As_RPKM > 6))

# How many DEGs with abs(log2FoldChange) > 1, RPKM > 6 have any GO term annotation?

# Violin plots of BP, MF, CC terms and number of unique gene counts




############################################## GO analysis
GO_term_id <- read_tsv("https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt", col_names = F, skip=4)
View(GO_term_id)


sig_WER_down <- sig_WER_all_info %>% 
  filter(log2FoldChange < -1.3, WER_0As_RPKM > 8)
sig_Cortex_down <- sig_Cortex_all_info %>% 
  filter(log2FoldChange < -1.3, Cortex_0As_RPKM > 8)
sig_SCR_down <- sig_SCR_all_info %>% 
  filter(log2FoldChange < -1.3, SCR_0As_RPKM > 8)

sig_WER_up <- sig_WER_all_info %>% 
  filter(log2FoldChange > 1.3, WER_50As_RPKM > 8)
sig_Cortex_up <- sig_Cortex_all_info %>% 
  filter(log2FoldChange > 1.3, Cortex_50As_RPKM > 8)
sig_SCR_up <- sig_SCR_all_info %>% 
  filter(log2FoldChange > 1.3, SCR_50As_RPKM > 8)

# summary table of sig up/down genes per cell type
GO_sum_table <-tibble('Cell-Type'=c('Epidermal','Cortex','Endodermal'),'Up-Reg'=c(2564,2459,2145),'Down-Reg'=c(2465,2042,2454))
grid.table(GO_sum_table, rows = NULL)
#save from image

# Make enrichResult objects for up/down, cell-type, & subontology categories (Biological process = BP, Molecular Function = MF, Cellular Components = CC)
# GO_WER_down_bp, GO_WER_down_mf, GO_WER_down_cc, GO_WER_up_bp, GO_WER_up_mf, GO_WER_up_cc, ...18 total objects
#### BUG #### In order to change ont ('MF', 'BP', 'CC'), you need to run the example org.Hs.eg.db below, then go back to our org.At.tair.db data
GO_Cortex_up_bp <- clusterProfiler::enrichGO(
                   gene          = sig_Cortex_up$ENTREZID,
                   OrgDb         = 'org.At.tair.db',
                   ont           = "BP", # BP, MF, or CC 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.01, 
                   readable      = T,
                   pool = F)
View(GO_WER_down_bp)

library(org.Hs.eg.db)
data(geneList, package = "DOSE")
de <- names(geneList)[1:100]

yy <- enrichGO(de, 'org.Hs.eg.db', ont="CC", pvalueCutoff=0.01)
yy <- as_tibble(yy)
View(yy)   

# class(GO_WER_up_bp)
# head(GO_WER_down)
### the merge_results function doesn't seem to be working in the clusterProfiler package, so we need to compile the results manually
# 1) change each enrichResult object to a tibble 
# 2) create new columns for groups (cell_type = c(Epidermal,Cortex,Endodermal), subontology = c(Biological_Process, Molecular_Function, Cellular_Components), regulation = c(UP, DOWN))
# 3) paste tables together 
# 4) ggplot with facets

 GO_Cortex_down_bp <- as_tibble(GO_Cortex_down_bp) %>%
   mutate(cell_type = 'Cortex', Regulation = 'Down', subontology = 'Biological_Process')
 GO_Cortex_down_mf <- as_tibble(GO_Cortex_down_mf) %>%
   mutate(cell_type = 'Cortex', Regulation = 'Down', subontology = 'Molecular_Function')
 GO_Cortex_down_cc <- as_tibble(GO_Cortex_down_cc) %>%
   mutate(cell_type = 'Cortex', Regulation = 'Down', subontology = 'Cellular_Component')
 GO_Cortex_up_bp <- as_tibble(GO_Cortex_up_bp) %>%
   mutate(cell_type = 'Cortex', Regulation = 'Up', subontology = 'Biological_Process')
 GO_Cortex_up_mf <- as_tibble(GO_Cortex_up_mf) %>%
   mutate(cell_type = 'Cortex', Regulation = 'Up', subontology = 'Molecular_Function')
 GO_Cortex_up_cc <- as_tibble(GO_Cortex_up_cc) %>%
   mutate(cell_type = 'Cortex', Regulation = 'Up', subontology = 'Cellular_Component')

 GO_WER_down_bp <- as_tibble(GO_WER_down_bp) %>%
   mutate(cell_type = 'Epidermal', Regulation = 'Down', subontology = 'Biological_Process')
 GO_WER_down_mf <- as_tibble(GO_WER_down_mf) %>%
   mutate(cell_type = 'Epidermal', Regulation = 'Down', subontology = 'Molecular_Function')
 GO_WER_down_cc <- as_tibble(GO_WER_down_cc) %>%
   mutate(cell_type = 'Epidermal', Regulation = 'Down', subontology = 'Cellular_Component')
 GO_WER_up_bp <- as_tibble(GO_WER_up_bp) %>%
   mutate(cell_type = 'Epidermal', Regulation = 'Up', subontology = 'Biological_Process')
 GO_WER_up_mf <- as_tibble(GO_WER_up_mf) %>%
   mutate(cell_type = 'Epidermal', Regulation = 'Up', subontology = 'Molecular_Function')
 GO_WER_up_cc <- as_tibble(GO_WER_up_cc) %>%
   mutate(cell_type = 'Epidermal', Regulation = 'Up', subontology = 'Cellular_Component')

 GO_SCR_down_bp <- as_tibble(GO_SCR_down_bp) %>%
   mutate(cell_type = 'Endodermal', Regulation = 'Down', subontology = 'Biological_Process')
 GO_SCR_down_mf <- as_tibble(GO_SCR_down_mf) %>%
   mutate(cell_type = 'Endodermal', Regulation = 'Down', subontology = 'Molecular_Function')
 GO_SCR_down_cc <- as_tibble(GO_SCR_down_cc) %>%
   mutate(cell_type = 'Endodermal', Regulation = 'Down', subontology = 'Cellular_Component')
 GO_SCR_up_bp <- as_tibble(GO_SCR_up_bp) %>%
   mutate(cell_type = 'Endodermal', Regulation = 'Up', subontology = 'Biological_Process')
 GO_SCR_up_mf <- as_tibble(GO_SCR_up_mf) %>%
   mutate(cell_type = 'Endodermal', Regulation = 'Up', subontology = 'Molecular_Function')
 GO_SCR_up_cc <- as_tibble(GO_SCR_up_cc) %>%
   mutate(cell_type = 'Endodermal', Regulation = 'Up', subontology = 'Cellular_Component')

full_GO <- bind_rows(GO_Cortex_down_bp, GO_Cortex_down_cc, GO_Cortex_down_mf,
                     GO_Cortex_up_bp,   GO_Cortex_up_cc,   GO_Cortex_up_mf, 
                     GO_WER_down_bp,    GO_WER_down_cc,    GO_WER_down_mf, 
                     GO_WER_up_bp,      GO_WER_up_cc,      GO_WER_up_mf, 
                     GO_SCR_down_bp,    GO_SCR_down_cc,    GO_SCR_down_mf, 
                     GO_SCR_up_bp,      GO_SCR_up_cc,      GO_SCR_up_mf, 
                     .id = "dataset")
View(full_GO)
full_GO$cell_type <- as_factor(full_GO$cell_type, levels())
full_GO$cell_type = factor(full_GO$cell_type, levels=c("Epidermal","Cortex", "Endodermal"))
full_GO$subontology = factor(full_GO$subontology, levels=c("Biological_Process","Cellular_Component", "Molecular_Function"))
full_GO$Regulation = factor(full_GO$Regulation, levels=c("Up","Down"))


full_GO <- full_GO %>% 
  arrange(subontology, Regulation, cell_type, desc(Count)) %>% 
  select(subontology, Regulation, cell_type, Count, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID)

write_csv(full_GO,'full_GO.csv',col_names = T)

# need this function to make a reverse log scale x-axis
library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

GO_cc_up <- ggplot(filter(full_GO, subontology == 'Cellular_Component', Regulation == 'Up'),#Description %in% cc_up_categories), 
                   aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
  geom_point(aes(size = Count, color = cell_type),position = position_jitter(height = 0.25)) +
  theme_dark(base_size = 12) +
  ylab(NULL) +
  scale_x_continuous(trans=reverselog_trans(10) ) +
  guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        panel.grid.minor = element_blank(),
        legend.position = "left") +
  labs(title = "Cellular Component",
       subtitle = "Up-regulated Genes")

GO_cc_down <- ggplot(filter(full_GO, subontology == 'Cellular_Component', Regulation == 'Down'),#Description %in% cc_down_categories), 
                   aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
  geom_point(aes(size = Count, color = cell_type), position = position_jitter(height = 0.25)) +
  theme_dark(base_size = 12) +
  ylab(NULL) +
  scale_x_continuous(trans=reverselog_trans(10)) +
  scale_y_discrete(position = 'right') +
  guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(title = "Cellular Component",
       subtitle = "Down-regulated Genes") 

GO_cc <- GO_cc_up+GO_cc_down
ggsave('GO_cc_unfiltered.tiff', plot=last_plot(), device = 'tiff', dpi = 400)
cc_up_categories <- c('ribosome','nuclear lumen','vacuolar membrane','nucleolus','ribosomal subunit',
                      'mitochondrial part','peroxisome','cytosolic small ribosomal subunit',
                      'cytosolic large ribosomal subunit','proteasome complex') # 10 terms

cc_down_categories <- c('vacuolar membrane','Golgi apparatus part','endosome','trans-Golgi network',
                        'intrinsic component of membrane','plant-type cell wall','anchored component of membrane',
                        'apoplast','mitochondrial membrane','integral component of membrane') # 10 terms
#ggsave('GO_cc.tiff', plot=GO_cc, width = 10.94, height = 3.4,  device = 'tiff', dpi = 'retina')

############ molecular function
GO_mf_up <- ggplot(filter(full_GO, subontology == 'Molecular_Function', Regulation == 'Up'), #Description %in% mf_up_categories), 
                   aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
  geom_point(aes(size = Count, color = cell_type),position = position_jitter(height = 0.25)) +
  theme_dark(base_size = 12) +
  ylab(NULL) +
  scale_x_continuous(trans=reverselog_trans(10) ) +
  scale_y_discrete(position = 'left') +
  guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        panel.grid.minor = element_blank(),
        legend.position = "left") +
  labs(title = "Molecular Function",
       subtitle = "Up-regulated Genes")

GO_mf_down <- ggplot(filter(full_GO, subontology == 'Molecular_Function', Regulation == 'Down'), #Description %in% mf_down_categories), 
                     aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
  geom_point(aes(size = Count, color = cell_type), position = position_jitter(height = 0.25)) +
  theme_dark(base_size = 12) +
  ylab(NULL) +
  scale_x_continuous(trans=reverselog_trans(10)) +
  scale_y_discrete(position = 'right') +
  guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(title = "Molecular Function",
       subtitle = "Down-regulated Genes") 

GO_mf <- GO_mf_up+GO_mf_down
ggsave('GO_mf_unfiltered.tiff', plot=last_plot(), device = 'tiff', dpi = 400)
# ggsave('GO_mf_up.tiff', plot=GO_mf_up, width = 5.5, height = 2.5,  device = 'tiff', dpi = 'retina')
# ggsave('GO_mf_down.tiff', plot=GO_mf_down, width = 7, height = 3.38,  device = 'tiff', dpi = 'retina')
# these two plot width/heights allow the grids to match

mf_up_categories <- c('structural molecule activity','structural constituent of ribosome','ATP binding',
                      'ubiquitin-protein transferase activity','glutathione transferase activity',
                      'protein heterodimerization activity') # 6 terms

mf_down_categories <- c('UDP-glycosyltransferase activity','copper ion binding','isomerase activity','calcium ion binding',
                        'antioxidant activity','peroxidase activity','water transmembrane transporter activity',
                        'intramolecular oxidoreductase activity','structural constituent of cytoskeleton',
                        'polygalacturonate 4-alpha-galacturonosyltransferase activity',
                        'ATPase activity, couple to transmembrane movement of ions, rotational mechanism',
                        'protein disulfide isomerase activity') # 12 terms

########### biological process
GO_bp_up <- ggplot(filter(full_GO, subontology == 'Biological_Process', Regulation == 'Up' ), #Description %in% bp_up_categories), 
                   aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
  geom_point(aes(size = Count, color = cell_type),position = position_jitter(height = 0.25)) +
  theme_dark(base_size = 12) +
  ylab(NULL) +
  scale_x_continuous(trans=reverselog_trans(10) ) +
  guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        panel.grid.minor = element_blank(),
        legend.position = "left") +
  labs(title = "Biological Process",
       subtitle = "Up-regulated Genes")

GO_bp_down <- ggplot(filter(full_GO, subontology == 'Biological_Process', Regulation == 'Down' ), #Description %in% bp_down_categories), 
                     aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
  geom_point(aes(size = Count, color = cell_type), position = position_jitter(height = 0.25)) +
  theme_dark(base_size = 12) +
  ylab(NULL) +
  scale_x_continuous(trans=reverselog_trans(10)) +
  scale_y_discrete(position = 'right') +
  guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(title = "Biological Process",
       subtitle = "Down-regulated Genes") 

tiff("GO_bp_unfiltered.tiff", units="in", width=16, height=35, res=500)
GO_bp_up+GO_bp_down
dev.off()         

#tiff('GO_bp_down.tiff', units = 'in', width = 10, height = 35, res=300)
#tiff('GO_bp_up.tiff', units = 'in', width = 10, height = 35, res=300)
GO_bp_up+GO_bp_down
dev.off()         
ggsave("GO_bp.tiff", plot = last_plot(), device = "tiff", dpi = 'retina')

bp_down_categories <- c('root morphogenesis','response to cadmium ion','sulfur coumpound biosynthetic process',
                        'response to jasmonic acid','purine-containing compound metabolic process','cell maturation',
                        'ribonucleotide metabolic process','cation homeostasis','steroid metabolic process','response to water',
                        'aromatic compound catabolic process','response to extracellular stimulus',
                        'cellular nitrogen compound catabolic process','cell wall biogenesis','developmental cell growth',
                        'transition metal ion transport','response to wounding','response to organonitrogen compound',
                        'cellular amino acid biosynthetic process','cell death','sulfur amino acid metabolic process',
                        'response to nutrient levels','regulation of immune system process','cellular polysaccharide biosynthetic process',
                        'response to starvation','inorganic anion transport','Golgi organization','protein localization to membrane',
                        'carbohydrate catabolic process','nitrate transport','response to auxin','metal ion homeostasis',
                        'secondary metabolite biosynthetic process','divalent metal ion transport','brassinosteroid metabolic process',
                        'Golgi vesicle transport','response to ethylene','hyperosmotic response','glucose metabolic process',
                        'flavonoid metabolic process','cysteine metabolic process','organophosphate catabolic process',
                        'iron ion transport','ATP metabolic process','calcium ion transport','response to cytokinin','water transport',
                        'amino acid transport','positive regulation of flavonoid biosynthetic process')
# 48 terms above

bp_up_categories <- c('response to toxic substance','response to oxidative stress','response to reactive oxygen species',
                      'response to organonitrogen compound','response to cadmium ion','response to heat','response to topologically incorrect protein',
                      'cell death','organic acid catabolic process','intracellular signal transduction','response to light intensity',
                      'protein localization to organelle','protein catabolic process','toxin metabolic process','nucleotide biosynthetic process',
                      'response to water','response to salicylic acid','response to chitin','detoxification',
                      'response to endoplasmic reticulum stress','cellular protein catabolic process','toxin catabolic process',
                      'aromatic compound catabolic process','response to jasmonic acid','oxidation-reduction process',
                      'plant-type hypersensitive response','antibiotic metabolic process','fatty acid metabolic process',
                      'cellular response to drug','regulation of immune response','response to ethylene','defense response to bacterium',
                      'response to high light intensity','protein folding','defense response to fungus','response to wounding',
                      'response to hydrogen peroxide','protein localization to membrane','cellular response to organic cyclic compound',
                      'ubiquitin-dependdent protein catabolic process','protein targeting to membrane','cellular response to lipid',
                      'lipid catabolic process','cellular response to antibiotic','reactive oxygen species metabolic process',
                      'Golgi vesicle transport','organic acid transport','vacuolar transport','autophagy')
# 48 terms above

########## Heatmaps ###########################################################################################

# Use each sig____all_info_filter to only include abs(log2FC) > 1 and RPKM > 6

heatmap_sigWER  <- sig_WER_all_info_filter %>% 
  filter(abs(log2FoldChange) > 1, (WER_0As_RPKM > 6 | WER_50As_RPKM > 6), padj < 0.0005)#, !is.na(symbol_list))

heatmap_sigCortex  <- sig_Cortex_all_info_filter %>% 
  filter(abs(log2FoldChange) > 1, (Cortex_0As_RPKM > 6 | Cortex_50As_RPKM > 6), padj < 0.0005)#, !is.na(symbol_list)) 

heatmap_sigSCR  <- sig_SCR_all_info_filter %>% 
  filter(abs(log2FoldChange) > 1, (SCR_0As_RPKM > 6 | SCR_50As_RPKM > 6), padj < 0.0005)#, !is.na(symbol_list)) 

heatmap_sigALL <- RPKM_avg %>%
  filter(gene_id %in% c(heatmap_sigWER$gene_id, heatmap_sigCortex$gene_id, heatmap_sigSCR$gene_id))

heatmap_sig_mat <- heatmap_sigALL %>% 
  column_to_rownames(var = 'gene_id') %>%
  as.matrix()  

# function to generate parameter for scaling each row; each RPKM observation - mean(RPKM) across row, div by sd 
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
heatmap_sig_mat_norm <- t(apply(heatmap_sig_mat, 1, cal_z_score))

# pheatmap(heatmap_sig_mat_norm)

my_sample_col <- data.frame(condition = c('0As', '50As', '0As', '50As','0As', '50As'))
row.names(my_sample_col) <- colnames(heatmap_sig_mat)


my_hclust_gene <- hclust(dist(heatmap_sig_mat_norm), method = "ward.D2")
as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)

my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 8) 
my_gene_col <- as.data.frame(my_gene_col)
my_gene_col <- dplyr::rename(my_gene_col, cluster = my_gene_col) 

# to extract genes w/in each cluster, save the previous line as csv file.  order of clusters isn't intuitive
overall_clusters <- my_gene_col %>% 
  rownames_to_column() %>% 
  dplyr::rename(gene_id = rowname) %>% 
  arrange(cluster)
  
dplyr::count(overall_clusters,cluster)

overall_clusters_match_heatmap <- overall_clusters %>% 
  mutate(heatmap_cluster = ifelse(cluster == 3 , 1, ifelse(cluster == 4 , 2, ifelse(cluster == 8 , 3, ifelse(cluster == 7 , 4,  
                                  ifelse(cluster == 2 , 5, ifelse(cluster == 5 , 6, ifelse(cluster == 1 , 7, ifelse(cluster == 6 , 8, 
                                                                                                                    NA)))))))))
figure_2A_cluster_genes <- overall_clusters_match_heatmap %>% 
  dplyr::select(gene_id, heatmap_cluster)

dplyr::count(figure_2A_cluster_genes,heatmap_cluster)

write_csv(figure_2A_cluster_genes, 'figure_2A_cluster_genes.csv', col_names = T) 

# adding RPKM values to the gene_cluster table
figure_2A_cluster_genes_table <- figure_2A_cluster_genes %>% 
  left_join(RPKM_avg_round, by = 'gene_id') %>% 
  arrange(heatmap_cluster)
# write_csv(figure_2A_cluster_genes_table, 'figure_2A_cluster_genes_table.csv', col_names = TRUE) 





mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)

my_color = list(cluster = c(cluster1='#a6cee3',cluster2='#1f78b4',cluster3='#b2df8a',cluster4='#33a02c',cluster5='#1a1a1a',
                             cluster6='#e31a1c',cluster7='#fdbf6f',cluster8='#ff7f00',cluster9='#cab2d6',
                             cluster10='#6a3d9a',cluster11='#ffff99',cluster12='#b15928',cluster13='#878787'))
summarise(as_tibble(heatmap_sig_mat_norm), n())


#tiff("pheatmap_RPKM_top2205.tiff", units="in", width=7, height=10, res=600)
tiff("pheatmap_RPKM_top7543.tiff", units="in", width=4.5, height=7, res=900)

pheatmap(heatmap_sig_mat_norm, 
         annotation_row = my_gene_col, 
         annotation_col = my_sample_col, 
         color = mypalette, 
         annotation_colors = my_color,
         annotation_legend = F,
         cellwidth = 30, 
         clustering_method = 'ward.D2',
         show_rownames=F,
         cutree_rows = 8,
         cutree_cols = 2,
         legend_labels = c(-3,-2,-1,0,1,2,3),
         #main = 'AsIII Induced Differentially Expressed Genes\nby Cell-type and Condition',
         angle_col = 45,
         labels_col = c('Epidermal_0As', 'Epidermal_50As', 'Cortex_0As', 'Cortex_50As', 'Endodermal_0As', 'Endodermal_50As'))
# scale is sd from mean(RPKM).  for each gene, mean(RPKM) of each cell-type/condition.  6 total  
dev.off()         
# sampleinfo = tibble(
#   sample = c('WER_minus_As', 'WER_plus_As', 'Cortex_minus_As', 'Cortex_plus_As', 'SCR_minus_As', 'SCR_plus_As'),
#   celltype = c('WER', 'WER', 'Cortex', 'Cortex', 'SCR', 'SCR'),
#   condition = c('0As', '50As', '0As', '50As','0As', '50As'))

############ heatmaps for different gene families
AtMATEid <- read_delim("~/R/edgeR_Todd_RNAseq/DEG_Families/MATEid/AtMATEid.csv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

AtMATEid <- dplyr::rename(AtMATEid, gene_id = Gene_ID)
# need to add some DTX genes manually to this list.  Some are missing. (No DTX for 20, 31, 36, 38, 39, 42, 43, 44, 45, 46, 47)
# (35 - AT4G25640, FRD3 = AT3G08040, MATE = AT1G51340, AT2G38330, AT4G38380, EDS5 = AT4G39030, EDS5H = AT2G21340, AT1G11670 )  
Gene_Name <- c('DTX35','FRD3','MATE','AT2G38330','AT4G38380','EDS5','EDS5H','AT1G11670')
gene_id <- c('AT4G25640','AT3G08040','AT1G51340','AT2G38330','AT4G38380','AT4G39030','AT2G21340','AT1G11670')
other_MATEs <- tibble(gene_id,Gene_Name)
AtMATEid <- AtMATEid %>% 
  dplyr::union(other_MATEs)
View(AtMATEid)

results_WER2 <- dplyr::rename(results_WER2, gene_id = GeneID)
results_Cortex2 <- dplyr::rename(results_Cortex2, gene_id = GeneID)
results_SCR2 <- dplyr::rename(results_SCR2, gene_id = GeneID)

sigMATE_WER <- results_WER2 %>%
  semi_join(AtMATEid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigMATE_Cortex <- results_Cortex2 %>%
  semi_join(AtMATEid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)
  
sigMATE_SCR <- results_SCR2 %>%
  semi_join(AtMATEid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

sigMATE <- sigMATE_WER %>% 
  left_join(sigMATE_Cortex, by='gene_id') %>% 
  left_join(sigMATE_SCR, by = 'gene_id') %>%
  left_join(AtMATEid, by = 'gene_id') %>% 
  left_join(RPKM_avg_round, by = 'gene_id') %>% 
  filter(WER_minus_As >= 2 | Cortex_minus_As >= 2 | SCR_minus_As >= 2 | WER_plus_As >= 2 | Cortex_plus_As >= 2 | SCR_plus_As >= 2) %>% 
  dplyr::select(Gene_Name,log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  column_to_rownames('Gene_Name')

my_hclust_MATE <- hclust(dist(sigMATE), method = "ward.D2")
as.dendrogram(my_hclust_MATE) %>%
  plot(horiz = TRUE)

my_MATE_col <- cutree(tree = as.dendrogram(my_hclust_MATE), k=3) 
my_MATE_col <- as.data.frame(my_MATE_col)
my_MATE_col <- dplyr::rename(my_MATE_col, cluster = my_MATE_col)
my_color2 = list(cluster = c(cluster1='#4daf4a',cluster2='#984ea3',cluster3='#bababa'))
breaksList = seq(-8, 8, by = .16)
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)

tiff("MATE_heatmap2.tiff", units="in", width=7, height=7, res=800)
pheatmap(sigMATE, 
         color = mypalette,
         annotation_row = my_MATE_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 3,
         show_rownames=T,
         breaks = breaksList,
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'))
dev.off()
MATE_RPKM <- AtMATEid %>% 
  inner_join(RPKM_avg_round, by = 'gene_id') %>% 
  filter(WER_minus_As >= 2 | Cortex_minus_As >= 2 | SCR_minus_As >= 2 | WER_plus_As >= 2 | Cortex_plus_As >= 2 | SCR_plus_As >= 2)
#write_csv(MATE_RPKM, 'MATE_RPKM.csv', col_names = T)


######### ABC #############################

AtABCid <- read_csv("~/R/edgeR_Todd_RNAseq/DEG_Families/ABCid/AtABCid.csv")
AtABCid <- dplyr::rename(AtABCid, gene_id = Gene_ID)

sigABC_WER <- results_WER2 %>%
  semi_join(AtABCid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigABC_Cortex <- results_Cortex2 %>%
  semi_join(AtABCid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigABC_SCR <- results_SCR2 %>%
  semi_join(AtABCid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

# sigABC <- sigABC_WER %>% 
#   left_join(sigABC_Cortex, by='gene_id') %>% 
#   left_join(sigABC_SCR, by = 'gene_id') %>%
#   left_join(AtABCid, by = 'gene_id') %>% 
#   column_to_rownames('Gene_Name') %>% 
#   drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
#   dplyr::select(-gene_id)
sigABC <- sigABC_WER %>% 
  left_join(sigABC_Cortex, by='gene_id') %>% 
  left_join(sigABC_SCR, by = 'gene_id') %>%
  left_join(AtABCid, by = 'gene_id') %>% 
  left_join(RPKM_avg_round, by = 'gene_id') %>% 
  filter(WER_minus_As > 2 | Cortex_minus_As > 2 | SCR_minus_As > 2 | WER_plus_As > 2 | Cortex_plus_As > 2 | SCR_plus_As > 2) %>% 
  dplyr::select(Gene_Name,log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  column_to_rownames('Gene_Name')

my_hclust_ABC <- hclust(dist(sigABC), method = "ward.D2")
as.dendrogram(my_hclust_ABC) %>%
  plot(horiz = TRUE)
my_ABC_col <- cutree(tree = as.dendrogram(my_hclust_ABC), k=4) 
my_ABC_col <- as.data.frame(my_ABC_col)
my_ABC_col <- dplyr::rename(my_ABC_col, cluster = my_ABC_col)
# to extract genes w/in each cluster, save the previous line as csv file.  order of clusters isn't intuitive
# write_csv(rownames_to_column(my_ABC_col), 'ABC_clusters.csv', col_names = T) 
my_color2 = list(cluster = c(cluster1='#4daf4a',cluster2='#984ea3',cluster3='#bababa', cluster4='#f1b6da'))
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)
breaksList = seq(-7, 7, by = .14)

tiff("ABC_heatmap_2.tiff", units="in", width=7, height=9.5, res=800)
pheatmap(sigABC, 
         color = mypalette, 
         annotation_row = my_ABC_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 4,
         show_rownames=T,
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10)
dev.off()

ABC_RPKM <- AtABCid %>% 
  inner_join(RPKM_avg_round, by = 'gene_id')
#write_csv(ABC_RPKM, 'ABC_RPKM.csv', col_names = T)

######## Aquaporins ##############################

AtAQUAid <- read_csv("~/R/edgeR_Todd_RNAseq/DEG_Families/AQUAPORINid/AtAQUAid.csv")
AtAQUAid <- dplyr::rename(AtAQUAid, gene_id = Gene_ID)
AtAQUAid[19, "Gene_Name"] <- 'NIP1;1'

sigAQUA_WER <- results_WER2 %>%
  semi_join(AtAQUAid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigAQUA_Cortex <- results_Cortex2 %>%
  semi_join(AtAQUAid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigAQUA_SCR <- results_SCR2 %>%
  semi_join(AtAQUAid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

sigAQUA <- sigAQUA_WER %>% 
  left_join(sigAQUA_Cortex, by='gene_id') %>% 
  left_join(sigAQUA_SCR, by = 'gene_id') %>%
  left_join(AtAQUAid, by = 'gene_id') %>% 
  left_join(RPKM_avg_round, by = 'gene_id') %>% 
  filter(WER_minus_As > 2 | Cortex_minus_As > 2 | SCR_minus_As > 2 | WER_plus_As > 2 | Cortex_plus_As > 2 | SCR_plus_As > 2) %>% 
  dplyr::select(Gene_Name,log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  column_to_rownames('Gene_Name')

my_hclust_AQUA <- hclust(dist(sigAQUA), method = "ward.D2")
as.dendrogram(my_hclust_AQUA) %>%
  plot(horiz = TRUE)

my_AQUA_col <- cutree(tree = as.dendrogram(my_hclust_AQUA), k=4) 
my_AQUA_col <- as.data.frame(my_AQUA_col)
my_AQUA_col <- dplyr::rename(my_AQUA_col, cluster = my_AQUA_col)
my_color2 = list(cluster = c(cluster1='#4daf4a',cluster2='#984ea3',cluster3='#bababa', cluster4='#f1b6da'))

tiff("AQUA_heatmap2.tiff", units="in", width=7, height=5.5, res=800)
pheatmap(sigAQUA, 
         color = mypalette, 
         annotation_row = my_AQUA_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 4,
         show_rownames=T,
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10)
dev.off()

AtAQUAid <- dplyr::rename(AtAQUAid, gene_id = Gene_ID)

Aqua_clusters <- AtAQUAid %>% 
  inner_join(RPKM_avg, by ='gene_id') %>% 
  mutate(avgRPKM_WER = round((WER_minus_As + WER_plus_As)/2,digits =0),
         avgRPKM_Cortex = round((Cortex_minus_As + Cortex_plus_As)/2,digits =0),
         avgRPKM_SCR = round((SCR_minus_As + SCR_plus_As)/2,digits =0)) %>% 
  dplyr::select(gene_id, Gene_Name, avgRPKM_WER, avgRPKM_Cortex, avgRPKM_SCR) 

Aqua_clusters <- AtAQUAid %>% 
  inner_join(RPKM_avg, by ='gene_id') %>% 
  purrr::modify_if(is.numeric, ~round(., 0)) 
  

aqua_cluster_1 <- c("TIP1;1","PIP2A","PIP1;4",'NIP1;1','TIP4;1','PIP2;5','NIP3;1','SIP2;1','NIP5;1','NIP6;1','PIP1B')
aqua_cluster_2 <- c("PIP2E","TIP2;2","TIP2",'PIP1;5','PIP2;4','RD28','DELTA-TIP','PIP1A','PIP2B','TIP2;3')
aqua_cluster_3 <- c("SIP1;2","BRX","NIP2;1")
aqua_cluster_4 <- c("BETA-TIP","NIP7;1","TIP1;3",'NIP1;2','PIP1C','PIP2;8','PIP3','SIP1A','TIP3;1','TIP5;1')

aqua_cl_1_RPKM <- Aqua_clusters %>%
  arrange(factor(Gene_Name, levels = aqua_cluster_1)) %>% 
  dplyr::slice(1:11)
aqua_cl_2_RPKM <- Aqua_clusters %>%
  arrange(factor(Gene_Name, levels = aqua_cluster_2)) %>% 
  dplyr::slice(1:10)
aqua_cl_3_RPKM <- Aqua_clusters %>%
  arrange(factor(Gene_Name, levels = aqua_cluster_3)) %>% 
  dplyr::slice(1:3)
aqua_cl_4_RPKM <- Aqua_clusters %>%
  arrange(factor(Gene_Name, levels = aqua_cluster_4)) %>% 
  dplyr::slice(1:10)

View(aqua_cl_4_RPKM)
grid.table(aqua_cl_1_RPKM, rows = NULL)
# save from R Graphics Device
grid.table(aqua_cl_2_RPKM, rows = NULL)
grid.table(aqua_cl_3_RPKM, rows = NULL)
grid.table(aqua_cl_4_RPKM, rows = NULL)

aqua_RPKM <- AtAQUAid %>% 
  inner_join(RPKM_avg_round, by = 'gene_id')
# write_csv(aqua_RPKM, 'aqua_RPKM.csv', col_names = T)
# format table in excel to line up with heatmap
########## antiporters ###########

AtANTIPORTERid <- read_csv("~/R/edgeR_Todd_RNAseq/DEG_Families/ANTIPORTERSid/AtANTIPORTERSid.csv")
AtANTIPORTERid <- dplyr::rename(AtANTIPORTERid, gene_id = Gene_ID)

sigANTIPORTER_WER <- results_WER2 %>%
  semi_join(AtANTIPORTERid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigANTIPORTER_Cortex <- results_Cortex2 %>%
  semi_join(AtANTIPORTERid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigANTIPORTER_SCR <- results_SCR2 %>%
  semi_join(AtANTIPORTERid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

sigANTIPORTER <- sigANTIPORTER_WER %>% 
  left_join(sigANTIPORTER_Cortex, by='gene_id') %>% 
  left_join(sigANTIPORTER_SCR, by = 'gene_id') %>%
  left_join(AtANTIPORTERid, by = 'gene_id') %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id)

my_hclust_ANTIPORTER <- hclust(dist(sigANTIPORTER), method = "ward.D2")

my_ANTIPORTER_col <- cutree(tree = as.dendrogram(my_hclust_ANTIPORTER), k=5) 
my_ANTIPORTER_col <- as.data.frame(my_ANTIPORTER_col)
my_ANTIPORTER_col <- dplyr::rename(my_ANTIPORTER_col, cluster = my_ANTIPORTER_col)
my_color2 = list(cluster = c(cluster1='#4daf4a',cluster2='#984ea3',cluster3='#bababa', cluster4='#f1b6da', cluster5='#fe9929'))

tiff("ANTIPORTER_heatmap.tiff", units="in", width=7, height=9, res=600)
pheatmap(sigANTIPORTER, 
         color = mypalette, 
         annotation_row = my_ANTIPORTER_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 5,
         show_rownames=T,
         #main = 'ANTIPORTER Expression (log2FC) after AsIII Stress',
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10)
dev.off()

ANTIPORTER_RPKM <- sigANTIPORTER_WER %>% 
  left_join(sigANTIPORTER_Cortex, by='gene_id') %>% 
  left_join(sigANTIPORTER_SCR, by = 'gene_id') %>%
  left_join(AtANTIPORTERid, by = 'gene_id') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  left_join(RPKM_avg_round, by ='gene_id') %>% 
  dplyr::filter(WER_minus_As > 2 | Cortex_minus_As > 2 | SCR_minus_As > 2 | WER_plus_As > 2 | Cortex_plus_As > 2 | SCR_plus_As > 2)

################# bhlh transcription factors

AtBHLHid <- read_csv("~/R/edgeR_Todd_RNAseq/DEG_Families/BHLHid/AtBHLHid.csv")
AtBHLHid <- dplyr::rename(AtBHLHid, gene_id = Gene_ID)

sigBHLH_WER <- results_WER2 %>%
  semi_join(AtBHLHid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigBHLH_Cortex <- results_Cortex2 %>%
  semi_join(AtBHLHid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigBHLH_SCR <- results_SCR2 %>%
  semi_join(AtBHLHid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

sigBHLH <- sigBHLH_WER %>% 
  left_join(sigBHLH_Cortex, by='gene_id') %>% 
  left_join(sigBHLH_SCR, by = 'gene_id') %>%
  left_join(AtBHLHid, by = 'gene_id') %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id)

my_hclust_BHLH <- hclust(dist(sigBHLH), method = "ward.D2")
as.dendrogram(my_hclust_BHLH) %>%
  plot(horiz = TRUE)
# 10 clusters looks best imo
my_BHLH_col <- cutree(tree = as.dendrogram(my_hclust_BHLH), k=10) 
my_BHLH_col <- as.data.frame(my_BHLH_col)
my_BHLH_col <- dplyr::rename(my_BHLH_col, cluster = my_BHLH_col)
my_color2 = list(cluster = c(cluster1='#a6cee3',cluster2='#1f78b4',cluster3='#b2df8a',cluster4='#33a02c',cluster5='#1a1a1a',
                             cluster6='#e31a1c',cluster7='#fdbf6f',cluster8='#ff7f00',cluster9='#cab2d6',
                             cluster10='#6a3d9a'))
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)

tiff("BHLH_heatmap_10.tiff", units="in", width=7, height=22, res=300)
pheatmap(sigBHLH, 
         color = mypalette, 
         annotation_row = my_BHLH_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 10,
         show_rownames=T,
         #main = 'BHLH Expression (log2FC) after AsIII Stress',
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10)
dev.off()

BHLH_RPKM <- AtBHLHid %>% 
  inner_join(RPKM_avg_round, by = 'gene_id')
#write_csv(BHLH_RPKM, 'BHLH_RPKM.csv', col_names = T)

################################### WRKY transcription factors
AtWRKYid <- read_csv("~/R/edgeR_Todd_RNAseq/DEG_Families/WRKYid/AtWRKYid.csv")
AtWRKYid <- dplyr::rename(AtWRKYid, gene_id = Gene_ID)

sigWRKY_WER <- results_WER2 %>%
  semi_join(AtWRKYid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigWRKY_Cortex <- results_Cortex2 %>%
  semi_join(AtWRKYid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigWRKY_SCR <- results_SCR2 %>%
  semi_join(AtWRKYid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

target <- c("Tom", "Lynn")
filter(dat, name %in% target)  # equivalently, dat %>% filter(name %in% target)


candidates <- c('WRKY4', 'WRKY6', 'WRKY11', 'WRKY15', 'WRKY20', 'WRKY21', 'WRKY32', 'WRKY33', 'WRKY30', 'WRKY39', 'WRKY75')
sigWRKY <- sigWRKY_WER %>% 
  left_join(sigWRKY_Cortex, by='gene_id') %>% 
  left_join(sigWRKY_SCR, by = 'gene_id') %>%
  left_join(AtWRKYid, by = 'gene_id') %>%
  dplyr::filter(Gene_Name %in% candidates) %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id) 



my_hclust_WRKY <- hclust(dist(sigWRKY), method = "ward.D2")
as.dendrogram(my_hclust_WRKY) %>%
  plot(horiz = TRUE)
# 7 clusters looks best imo
my_WRKY_col <- cutree(tree = as.dendrogram(my_hclust_WRKY), k=2) 
my_WRKY_col <- as.data.frame(my_WRKY_col)
my_WRKY_col <- dplyr::rename(my_WRKY_col, cluster = my_WRKY_col)
breaksList = seq(-5, 5, by = .1)
my_color2 = list(cluster = c(cluster1='#a6cee3',cluster2='#1f78b4',cluster3='#b2df8a',cluster4='#33a02c',cluster5='#1a1a1a',
                             cluster6='#e31a1c',cluster7='#fdbf6f'))
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)

tiff("WRKY_heatmap_7.tiff", units="in", width=7, height=3.2, res=600)
pheatmap(sigWRKY, 
         color = mypalette, 
         annotation_row = my_WRKY_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 2,
         show_rownames=T,
         #main = 'WRKY Expression (log2FC) after AsIII Stress',
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10,
         breaks = breaksList)
dev.off()

WRKY_RPKM <- sigWRKY_WER %>% 
  left_join(sigWRKY_Cortex, by='gene_id') %>% 
  left_join(sigWRKY_SCR, by = 'gene_id') %>%
  left_join(AtWRKYid, by = 'gene_id') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  left_join(RPKM_avg_round, by ='gene_id') %>% 
  dplyr::filter(Gene_Name %in% candidates) %>% 
  dplyr::filter(WER_minus_As > 2 | Cortex_minus_As > 2 | SCR_minus_As > 2 | WER_plus_As > 2 | Cortex_plus_As > 2 | SCR_plus_As > 2)

#write_csv(WRKY_RPKM, 'WRKY_RPKM.csv', col_names = T)

####################### MYB transcription factors

AtMYBid <- read_csv("~/R/edgeR_Todd_RNAseq/DEG_Families/MYBid/AtMYBid.csv")
AtMYBid <- dplyr::rename(AtMYBid, gene_id = Gene_ID)

sigMYB_WER <- results_WER2 %>%
  semi_join(AtMYBid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigMYB_Cortex <- results_Cortex2 %>%
  semi_join(AtMYBid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigMYB_SCR <- results_SCR2 %>%
  semi_join(AtMYBid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

sigMYB <- sigMYB_WER %>% 
  left_join(sigMYB_Cortex, by='gene_id') %>% 
  left_join(sigMYB_SCR, by = 'gene_id') %>%
  left_join(AtMYBid, by = 'gene_id') %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id)
candidates <- c('MYB13', 'MYB30', 'MYB38', 'MYB45', 'MYB51', 'MYB96', 'MYB109', 'MYB1')
sigMYB <- sigMYB_WER %>% 
  left_join(sigMYB_Cortex, by='gene_id') %>% 
  left_join(sigMYB_SCR, by = 'gene_id') %>%
  left_join(AtMYBid, by = 'gene_id') %>%
  dplyr::filter(Gene_Name %in% candidates) %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id) 

my_hclust_MYB <- hclust(dist(sigMYB), method = "ward.D2")
as.dendrogram(my_hclust_MYB) %>%
  plot(horiz = TRUE)
# 7 clusters looks best imo
my_MYB_col <- cutree(tree = as.dendrogram(my_hclust_MYB), k=1) 
my_MYB_col <- as.data.frame(my_MYB_col)
my_MYB_col <- dplyr::rename(my_MYB_col, cluster = my_MYB_col)
my_color2 = list(cluster = c(cluster1='#a6cee3',cluster2='#1f78b4',cluster3='#b2df8a',cluster4='#33a02c',cluster5='#1a1a1a',
                             cluster6='#e31a1c',cluster7='#fdbf6f', cluster8='#ff7f00'))
breaksList = seq(-5, 5, by = .1)
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)

tiff("MYB_heatmap_7.tiff", units="in", width=7, height=2.5, res=300)
pheatmap(sigMYB, 
         color = mypalette, 
         annotation_row = my_MYB_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 1,
         show_rownames=T,
         #main = 'MYB Expression (log2FC) after AsIII Stress',
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10,
         breaks = breaksList)
dev.off()

MYB_RPKM <- sigMYB_WER %>% 
  left_join(sigMYB_Cortex, by='gene_id') %>% 
  left_join(sigMYB_SCR, by = 'gene_id') %>%
  left_join(AtMYBid, by = 'gene_id') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  left_join(RPKM_avg_round, by ='gene_id') %>% 
  dplyr::filter(Gene_Name %in% candidates) %>% 
  dplyr::filter(WER_minus_As > 2 | Cortex_minus_As > 2 | SCR_minus_As > 2 | WER_plus_As > 2 | Cortex_plus_As > 2 | SCR_plus_As > 2)
write_csv(MYB_RPKM, 'MYB_RPKM.csv', col_names = T)

################ NAC transcription factors

AtNACid <- read_csv("~/R/edgeR_Todd_RNAseq/DEG_Families/NACid/AtNACid.csv")
AtNACid <- dplyr::rename(AtNACid, gene_id = Gene_ID)

sigNAC_WER <- results_WER2 %>%
  semi_join(AtNACid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigNAC_Cortex <- results_Cortex2 %>%
  semi_join(AtNACid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigNAC_SCR <- results_SCR2 %>%
  semi_join(AtNACid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

sigNAC <- sigNAC_WER %>% 
  left_join(sigNAC_Cortex, by='gene_id') %>% 
  left_join(sigNAC_SCR, by = 'gene_id') %>%
  left_join(AtNACid, by = 'gene_id') %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id)

my_hclust_NAC <- hclust(dist(sigNAC), method = "ward.D2")
as.dendrogram(my_hclust_NAC) %>%
  plot(horiz = TRUE)
# 3 clusters looks best imo
my_NAC_col <- cutree(tree = as.dendrogram(my_hclust_NAC), k=3) 
my_NAC_col <- as.data.frame(my_NAC_col)
my_NAC_col <- dplyr::rename(my_NAC_col, cluster = my_NAC_col)
my_color2 = list(cluster = c(cluster1='#a6cee3',cluster2='#1f78b4',cluster3='#b2df8a'))
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)

tiff("NAC_heatmap_3.tiff", units="in", width=7, height=7, res=300)
pheatmap(sigNAC, 
         color = mypalette, 
         annotation_row = my_NAC_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 3,
         show_rownames=T,
         #main = 'NAC Expression (log2FC) after AsIII Stress',
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10)
dev.off()

NAC_RPKM <- AtNACid %>% 
  inner_join(RPKM_avg_round, by = 'gene_id')
#write_csv(NAC_RPKM, 'NAC_RPKM.csv', col_names = T)

############################

Gene_ID <- c('AT5G43350','AT5G43370', 'AT5G43360', 'AT2G38940','AT2G32830','AT5G43340','AT3G54700','AT1G20860','AT1G76430')#,
           #  'AT2G21045','AT5G44070','AT1G03980','AT1G30220','AT4G16480')
Gene_Name <- c('PHT1;1','PHT1;2','PHT1;3','PHT1;4','PHT1;5','PHT1;6','PHT1;7','PHT1;8','PHT1;9')#,'HAC1','PCS1','PCS2','INT2','INT4')
AtPHTid <- tibble(Gene_ID,Gene_Name)
AtPHTid <- dplyr::rename(AtPHTid, gene_id = Gene_ID)

sigPHT_WER <- results_WER2 %>%
  semi_join(AtPHTid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigPHT_Cortex <- results_Cortex2 %>%
  semi_join(AtPHTid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigPHT_SCR <- results_SCR2 %>%
  semi_join(AtPHTid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

sigPHT <- sigPHT_WER %>% 
  left_join(sigPHT_Cortex, by='gene_id') %>% 
  left_join(sigPHT_SCR, by = 'gene_id') %>%
  left_join(AtPHTid, by = 'gene_id') %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id)

my_hclust_PHT <- hclust(dist(sigPHT), method = "ward.D2")
as.dendrogram(my_hclust_PHT) %>%
  plot(horiz = TRUE)
# 2 clusters looks best imo
my_PHT_col <- cutree(tree = as.dendrogram(my_hclust_PHT), k=2) 
my_PHT_col <- as.data.frame(my_PHT_col)
my_PHT_col <- dplyr::rename(my_PHT_col, cluster = my_PHT_col)
my_color2 = list(cluster = c(cluster1='#b2df8a',cluster2='#ff7f00'))
breaksList = seq(-5, 5, by = .1)
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)


tiff("PHT_heatmap_2.tiff", units="in", width=7, height=2.4, res=300)
pheatmap(sigPHT, 
         color = mypalette, 
         annotation_row = my_PHT_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 2,
         show_rownames=T,
         #main = 'PHT Expression (log2FC) after AsIII Stress',
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10,
         breaks = breaksList
         )
dev.off()

PHT_RPKM <- AtPHTid %>% 
  inner_join(RPKM_avg_round, by = 'gene_id')
write_csv(PHT_RPKM, 'PHT_RPKM.csv', col_names = T)

######################################## Random arsenic genes

# maybe add At2g43330 INT1, At2g35740 INT3

Gene_ID <- c('AT2G21045','AT5G44070', 'AT5G43350', 'AT2G38940', 'AT4G02050', 'AT5G26340', 'AT5G62890', 'AT5G41610', 'AT5G13490', 'AT4G28390','AT5G17400','AT4G02050', 'AT5G26340')
Gene_Name <- c('HAC1','PCS1', 'PHT1;1', 'PHT1;4', 'STP7', 'STP13', 'NAT6', 'CHX18', 'AT5G13490', 'AT4G28390', 'AT5G17400','PLT7', 'PLT13')
AtRandomid <- tibble(Gene_ID,Gene_Name)
AtRandomid <- dplyr::rename(AtRandomid, gene_id = Gene_ID)

sigRandom_WER <- results_WER2 %>%
  semi_join(AtRandomid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigRandom_Cortex <- results_Cortex2 %>%
  semi_join(AtRandomid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigRandom_SCR <- results_SCR2 %>%
  semi_join(AtRandomid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

candidates <- c('HAC1','PCS1', 'PHT1;1', 'PHT1;4', 'STP7', 'STP13', 'NAT6', 'CHX18', 'AT5G13490', 'AT4G28390', 'AT5G17400','PLT7', 'PLT13')

sigRandom <- sigRandom_WER %>% 
  left_join(sigRandom_Cortex, by='gene_id') %>% 
  left_join(sigRandom_SCR, by = 'gene_id') %>%
  left_join(AtRandomid, by = 'gene_id') %>%
  dplyr::filter(Gene_Name %in% candidates) %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id) 


my_hclust_Random <- hclust(dist(sigRandom), method = "ward.D2")
as.dendrogram(my_hclust_Random) %>%
  plot(horiz = TRUE)
# 1 clusters looks best imo
my_Random_col <- cutree(tree = as.dendrogram(my_hclust_Random), k=3) 
my_Random_col <- as.data.frame(my_Random_col)
my_Random_col <- dplyr::rename(my_Random_col, cluster = my_Random_col)
my_color2 = list(cluster = c(cluster1='#a6cee3',cluster2='#1f78b4',cluster3='#b2df8a'))
breaksList = seq(-5, 5, by = .1)
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)


tiff("Random_heatmap_2.tiff", units="in", width=7, height=3.5, res=600)
pheatmap(sigRandom, 
         color = mypalette, 
         annotation_row = my_Random_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 3,
         show_rownames=T,
         #main = 'Random Expression (log2FC) after AsIII Stress',
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10,
         breaks = breaksList)
dev.off()

Random_RPKM <- sigRandom_WER %>% 
  left_join(sigRandom_Cortex, by='gene_id') %>% 
  left_join(sigRandom_SCR, by = 'gene_id') %>%
  left_join(AtRandomid, by = 'gene_id') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  left_join(RPKM_avg_round, by ='gene_id') %>% 
  dplyr::filter(Gene_Name %in% candidates)  
write_csv(Random_RPKM, 'Random_RPKM.csv', col_names = T)

############################# Nucleobase Ascorbate Transporters
Gene_ID <- c('AT2G05760','AT2G34190','AT2G26510','AT1G49960','AT5G49990', 'AT5G62890','AT1G60030','AT1G10540','AT5G25420','AT1G65550','AT4G38050','AT2G27810')
Gene_Name <- c('NAT1','NAT2','NAT3','NAT4','NAT5','NAT6','NAT7','NAT8','NAT9','NAT10','NAT11','NAT12')
AtNATid <- tibble(Gene_ID,Gene_Name)
AtNATid <- dplyr::rename(AtNATid, gene_id = Gene_ID)

sigNAT_WER <- results_WER2 %>%
  semi_join(AtNATid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigNAT_Cortex <- results_Cortex2 %>%
  semi_join(AtNATid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigNAT_SCR <- results_SCR2 %>%
  semi_join(AtNATid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

sigNAT <- sigNAT_WER %>% 
  left_join(sigNAT_Cortex, by='gene_id') %>% 
  left_join(sigNAT_SCR, by = 'gene_id') %>%
  left_join(AtNATid, by = 'gene_id') %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id)

my_hclust_NAT <- hclust(dist(sigNAT), method = "ward.D2")
as.dendrogram(my_hclust_NAT) %>%
  plot(horiz = TRUE)
# 3 clusters looks best imo
my_NAT_col <- cutree(tree = as.dendrogram(my_hclust_NAT), k=3) 
my_NAT_col <- as.data.frame(my_NAT_col)
my_NAT_col <- dplyr::rename(my_NAT_col, cluster = my_NAT_col)
my_color2 = list(cluster = c(cluster1='#b2df8a',cluster2='#ff7f00',cluster3='#1a1a1a'))
breaksList = seq(-2.5, 2.5, by = .05)
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)


tiff("NAT_heatmap_3.tiff", units="in", width=7, height=3, res=300)
pheatmap(sigNAT, 
         color = mypalette, 
         annotation_row = my_NAT_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 3,
         show_rownames=T,
         #main = 'NAT Expression (log2FC) after AsIII Stress',
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10,
         breaks = breaksList
)
dev.off()

NAT_RPKM <- AtNATid %>% 
  inner_join(RPKM_avg_round, by = 'gene_id')
View(NAT_RPKM)

########################## Sucrose Transporters
Gene_ID <- c('AT1G71880','AT1G22710','AT2G02860','AT1G09960','AT1G71890', 'AT5G43610','AT1G66570','AT2G14670','AT5G06170')
Gene_Name <- c('SUC1','SUC2','SUC3','SUC4','SUC5','SUC6','SUC7','SUC8','SUC9')
AtSUCid <- tibble(Gene_ID,Gene_Name)
AtSUCid <- dplyr::rename(AtSUCid, gene_id = Gene_ID)

sigSUC_WER <- results_WER2 %>%
  semi_join(AtSUCid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigSUC_Cortex <- results_Cortex2 %>%
  semi_join(AtSUCid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigSUC_SCR <- results_SCR2 %>%
  semi_join(AtSUCid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

sigSUC <- sigSUC_WER %>% 
  left_join(sigSUC_Cortex, by='gene_id') %>% 
  left_join(sigSUC_SCR, by = 'gene_id') %>%
  left_join(AtSUCid, by = 'gene_id') %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id)

my_hclust_SUC <- hclust(dist(sigSUC), method = "ward.D2")
as.dendrogram(my_hclust_SUC) %>%
  plot(horiz = TRUE)
# 2 clusters looks best imo
my_SUC_col <- cutree(tree = as.dendrogram(my_hclust_SUC), k=2) 
my_SUC_col <- as.data.frame(my_SUC_col)
my_SUC_col <- dplyr::rename(my_SUC_col, cluster = my_SUC_col)
my_color2 = list(cluster = c(cluster1='#b2df8a',cluster2='#ff7f00'))
breaksList = seq(-2.5, 2.5, by = .05)
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)


tiff("SUC_heatmap_3.tiff", units="in", width=7, height=1.9, res=300)
pheatmap(sigSUC, 
         color = mypalette, 
         annotation_row = my_SUC_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 3,
         show_rownames=T,
         #main = 'SUC Expression (log2FC) after AsIII Stress',
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10,
         breaks = breaksList
)
dev.off()

SUC_RPKM <- AtSUCid %>% 
  inner_join(RPKM_avg_round, by = 'gene_id')
View(SUC_RPKM)

###############################  Monosaccharide (hexoses / pentoses)-H+ symporter family
Gene_ID <- c('AT1G11260','AT1G07340','AT5G61520','AT3G19930','AT1G34580', 'AT3G05960','AT4G02050','AT5G26250','AT1G50310','AT3G19940',
             'AT5G23270','AT4G21480','AT5G26340','AT1G77210')

Gene_Name <- c('STP1','STP2','STP3','STP4','STP5','STP6','STP7','STP8','STP9','STP10','STP11','STP12','STP13','STP14')
AtSTPid <- tibble(Gene_ID,Gene_Name)
AtSTPid <- dplyr::rename(AtSTPid, gene_id = Gene_ID)

sigSTP_WER <- results_WER2 %>%
  semi_join(AtSTPid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigSTP_Cortex <- results_Cortex2 %>%
  semi_join(AtSTPid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigSTP_SCR <- results_SCR2 %>%
  semi_join(AtSTPid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

sigSTP <- sigSTP_WER %>% 
  left_join(sigSTP_Cortex, by='gene_id') %>% 
  left_join(sigSTP_SCR, by = 'gene_id') %>%
  left_join(AtSTPid, by = 'gene_id') %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id)

my_hclust_STP <- hclust(dist(sigSTP), method = "ward.D2")
as.dendrogram(my_hclust_STP) %>%
  plot(horiz = TRUE)
# 4 clusters looks best imo
my_STP_col <- cutree(tree = as.dendrogram(my_hclust_STP), k=4) 
my_STP_col <- as.data.frame(my_STP_col)
my_STP_col <- dplyr::rename(my_STP_col, cluster = my_STP_col)
my_color2 = list(cluster = c(cluster1='#b2df8a',cluster2='#ff7f00', cluster3='#cab2d6',cluster4='#33a02c'))
breaksList = seq(-2.5, 2.5, by = .05)
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)


tiff("STP_heatmap_4.tiff", units="in", width=7, height=3.2, res=300)
pheatmap(sigSTP, 
         color = mypalette, 
         annotation_row = my_STP_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 4,
         show_rownames=T,
         #main = 'STP Expression (log2FC) after AsIII Stress',
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10,
         breaks = breaksList
)
dev.off()

STP_RPKM <- AtSTPid %>% 
  inner_join(RPKM_avg_round, by = 'gene_id')
View(STP_RPKM)

############################## H+-Symporter family for polyols and monosaccharides
Gene_ID <- c('AT1G11260','AT1G07340','AT5G61520','AT3G19930','AT1G34580', 'AT3G05960','AT4G02050','AT5G26250','AT1G50310','AT3G19940',
             'AT5G23270','AT4G21480','AT5G26340','AT1G77210')

Gene_Name <- c('PLT1','PLT2','PLT3','PLT4','PLT5','PLT6','PLT7','PLT8','PLT9','PLT10','PLT11','PLT12','PLT13','PLT14')
AtPLTid <- tibble(Gene_ID,Gene_Name)
AtPLTid <- dplyr::rename(AtPLTid, gene_id = Gene_ID)

sigPLT_WER <- results_WER2 %>%
  semi_join(AtPLTid, by ='gene_id') %>%
  dplyr::rename(log2FC_WER = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_WER)

sigPLT_Cortex <- results_Cortex2 %>%
  semi_join(AtPLTid, by ='gene_id') %>%
  dplyr::rename(log2FC_Cortex = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_Cortex)

sigPLT_SCR <- results_SCR2 %>%
  semi_join(AtPLTid, by ='gene_id') %>%
  dplyr::rename(log2FC_SCR = log2FoldChange) %>% 
  dplyr::select(gene_id, log2FC_SCR)

sigPLT <- sigPLT_WER %>% 
  left_join(sigPLT_Cortex, by='gene_id') %>% 
  left_join(sigPLT_SCR, by = 'gene_id') %>%
  left_join(AtPLTid, by = 'gene_id') %>% 
  column_to_rownames('Gene_Name') %>% 
  drop_na(log2FC_WER, log2FC_Cortex, log2FC_SCR) %>% 
  dplyr::select(-gene_id)

my_hclust_PLT <- hclust(dist(sigPLT), method = "ward.D2")
as.dendrogram(my_hclust_PLT) %>%
  plot(horiz = TRUE)
# 4 clusters looks best imo
my_PLT_col <- cutree(tree = as.dendrogram(my_hclust_PLT), k=4) 
my_PLT_col <- as.data.frame(my_PLT_col)
my_PLT_col <- dplyr::rename(my_PLT_col, cluster = my_PLT_col)
my_color2 = list(cluster = c(cluster1='#b2df8a',cluster2='#ff7f00', cluster3='#cab2d6',cluster4='#33a02c'))
breaksList = seq(-2.5, 2.5, by = .05)
mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)


tiff("PLT_heatmap_4.tiff", units="in", width=7, height=3.2, res=300)
pheatmap(sigPLT, 
         color = mypalette, 
         annotation_row = my_PLT_col, 
         annotation_legend = F,
         annotation_colors = my_color2,
         cluster_cols = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cutree_rows = 4,
         show_rownames=T,
         #main = 'PLT Expression (log2FC) after AsIII Stress',
         angle_col = 45,
         labels_col = c('Epidermal', 'Cortex', 'Endodermal'),
         fontsize=10,
         breaks = breaksList
)
dev.off()

PLT_RPKM <- AtPLTid %>% 
  inner_join(RPKM_avg_round, by = 'gene_id')
View(PLT_RPKM)

###############################
# Bar plot of WER, Cortex, and SCR FPKM units in control and AsIII conditions.  Mention in Discussion that WER and SCR are cell-type specific in both conditions, but Cortex is downregulated under AsIII.
# However, the GFP signal in the Cortex line after 24 hr AsIII was still visible and specific to cortex cells, albeit slightly weaker, and was therefor still capable of sorting cortex cells under AsIII conditions

SupFig2_se <- RPKM_complete %>% 
  dplyr::filter(gene_id %in% c('AT5G14750', 'AT1G09750', 'AT3G54220')) %>% 
  pivot_longer(-gene_id, names_to = 'sample2', values_to = 'RPKM') %>% 
  mutate(sample = case_when(grepl("Wer_minus", sample2) ~ "WER_0_As",
                            grepl("SCR_minus", sample2) ~ "SCR_0_As", 
                            grepl("Cortex_minus", sample2) ~ "Cortex_0_As",
                            grepl("Wer_plus", sample2) ~ "WER_50_As",
                            grepl("SCR_plus", sample2) ~ "SCR_50_As", 
                            grepl("Cortex_plus", sample2) ~ "Cortex_50_As")) %>% 
  dplyr::select(-sample2) %>% 
  group_by(gene_id, sample) %>%
  summarise(count = n(),
            value = 'RPKM',
            mean = mean(RPKM, na.rm = TRUE),
            sd = sd(RPKM, na.rm = TRUE), 
            se = sd/sqrt(count)) %>%
  ungroup() %>% 
  dplyr::rename(RPKM = mean)
  
View(SupFig2_se)

ggplot(data=SupFig2_se, aes(x=reorder(gene_id, -RPKM), y=RPKM, fill=sample, group=sample)) +
  geom_bar(colour="black", stat="identity", position=position_dodge(width=0.8), width=0.8) +
  scale_fill_brewer(palette="Paired") +
  geom_errorbar(aes(x=gene_id,group=sample, ymin = RPKM-se, ymax = RPKM+se), position=position_dodge(width=0.8), width = 0.6) +
  scale_x_discrete(labels=c("AT5G14750 \nWER::GFP","AT3G54220 \nSCR::GFP","AT1G09750 \nCortex::GFP"), name ="Gene") +
  theme(legend.position = c(0.9, 0.83),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))

ggsave("FACS_reporter_RPKM.tiff", plot = last_plot(), device = "tiff", dpi = 600)

###############################
# A. thaliana GO term annotation violin plots

MF_sort <- GO_term_id %>% 
  filter(X8 == "F") %>% 
  group_by(X5) %>% 
  summarise(Gene_Count = n_distinct(X1)) %>% 
  arrange(desc(Gene_Count)) %>% 
  mutate(Sub_ontology = "Molecular Function") %>% 
  dplyr::rename(Term = X5)
View(MF_sort)

BP_sort <- GO_term_id %>% 
  filter(X8 == "P") %>% 
  group_by(X5) %>% 
  summarise(Gene_Count = n_distinct(X1)) %>% 
  arrange(desc(Gene_Count)) %>% 
  mutate(Sub_ontology = "Biological Process") %>% 
  dplyr::rename(Term = X5)
View(BP_sort)


CC_sort <- GO_term_id %>% 
  filter(X8 == "C") %>% 
  group_by(X5) %>% 
  summarise(Gene_Count = n_distinct(X1)) %>% 
  arrange(desc(Gene_Count)) %>% 
  mutate(Sub_ontology = "Cellular Component") %>% 
  dplyr::rename(Term = X5)
View(CC_sort)

GO_term_count <- MF_sort %>% 
  bind_rows(BP_sort) %>% 
  bind_rows(CC_sort)
View(GO_term_count)

library(ggrepel)
GO_term_violin

label_terms <- c('biological_process', 'acetyl-CoA biosynthetic process from acetate', 'nucleus', 'eRF1 methyltransferase complex', 
                 'molecular_function', 'response to cadmium ion', 
                 'mitotic cell cycle')

ggplot(GO_term_count, aes(x=Sub_ontology, y=Gene_Count)) +
  geom_violin(aes(color=Sub_ontology), scale = 'count', size = 1.5) +
  geom_jitter(aes(color=Sub_ontology), alpha=0.4, position = position_jitter(width = 0.05)) +
  scale_y_log10(name = 'Gene Count') +
  scale_x_discrete(labels=c("Biological Process \nn=3,969","Cellular Component \nn=839","Molecular Function \nn=2,727"), name='Sub-Ontology') +
  theme(axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=11),
        legend.position = "none",
        title = element_text(size = 13, face = 'bold')) +
  labs(title = "A. thaliana GO term annotation",
       caption = "Data source: arabidopsis.org") +
  geom_label_repel(data = subset(GO_term_count, Term %in% label_terms), aes(label=Term, fill = Sub_ontology),
                   fontface='bold', color = 'white', size = 4.5,
                   box.padding = unit(0.25, "lines"),
                   point.padding = unit(1, "lines"),
                   segment.color = 'black', 
                   segment.size = 1,
                   force = 10, 
                   nudge_x = 0.15,
                   nudge_y = -0.15)

ggsave('Arabidopsis_GO_term_annotation.tiff', plot=last_plot(), device = 'tiff', dpi = 500)
