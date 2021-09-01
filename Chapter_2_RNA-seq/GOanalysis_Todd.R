# Run DESeq2_Dec_2019_2 script first

library(org.At.tair.db) # might need to click from Packages tab 

columns(org.At.tair.db)
keys <- head( keys(org.At.tair.db) )
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
  inner_join(FPKM_avg, by = 'gene_id') %>% 
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
                   WER_0As_FPKM = round(WER_minus_As, digits = 2), 
                   WER_50As_FPKM = round(WER_plus_As, digits = 2))

as.data.frame(sig_WER_all_info) %>%
  write.table("sig_WER_all_info")  # Saving all genes differentially expressed in WER

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
  inner_join(FPKM_avg, by = 'gene_id') %>% 
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
                   Cortex_0As_FPKM = round(Cortex_minus_As, digits = 2), 
                   Cortex_50As_FPKM = round(Cortex_plus_As, digits = 2))

as.data.frame(sig_Cortex_all_info) %>%
  write.table("sig_Cortex_all_info")  # Saving all genes differentially expressed in WER

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
  inner_join(FPKM_avg, by = 'gene_id') %>% 
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
                   SCR_0As_FPKM = round(SCR_minus_As, digits = 2), 
                   SCR_50As_FPKM = round(SCR_plus_As, digits = 2))

as.data.frame(sig_SCR_all_info) %>%
  write.table("sig_SCR_all_info")  # Saving all genes differentially expressed in WER



############################################## GO analysis
#BiocManager::install("DOSE")
library(DOSE)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)

OrgDb <-org.At.tair.db ##will need to use this later for GO analysis

sig_WER_all_info_filter <- sig_WER_all_info %>% 
  filter(abs(log2FoldChange) > 1.3, (WER_0As_FPKM > 8 | WER_50As_FPKM > 8))
sig_Cortex_all_info_filter <- sig_Cortex_all_info %>% 
  filter(abs(log2FoldChange) > 1.3, (Cortex_0As_FPKM > 8 | Cortex_50As_FPKM > 8))
sig_SCR_all_info_filter <- sig_SCR_all_info %>% 
  filter(abs(log2FoldChange) > 1.3, (SCR_0As_FPKM > 8 | SCR_50As_FPKM > 8))


geneList_WER <- sig_WER_all_info_filter$log2FoldChange  # filtered logFC > 1.3 & FPKM > 8
names(geneList_WER) <- sig_WER_all_info_filter$ENTREZID ##will use GO based on ENTREZ ID

geneList_Cortex <- sig_Cortex_all_info_filter$log2FoldChange  # filtered logFC > 1.3 & FPKM > 8
names(geneList_Cortex) <- sig_Cortex_all_info_filter$ENTREZID ##will use GO based on ENTREZ ID

geneList_SCR <- sig_SCR_all_info_filter$log2FoldChange  # filtered logFC > 1.3 & FPKM > 8
names(geneList_SCR) <- sig_SCR_all_info_filter$ENTREZID ##will use GO based on ENTREZ ID

ego <- clusterProfiler::enrichGO(gene          = names(geneList_SCR),
                                 OrgDb         = OrgDb,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05, 
                                 readable      = TRUE)

barplot(ego, drop=TRUE, showCategory=60)

clusterProfiler::dotplot(ego, showCategory=60) ##same but as doplot









##I prefer the above one for Go enrichemnt analysis.
# Group GO
ggo <- clusterProfiler::groupGO(gene     = names(geneList),
                                OrgDb    = OrgDb,
                                ont      = "BP",
                                level    = 3,
                                readable = TRUE)
plot(ggo)
barplot(ggo, drop=TRUE, showCategory=40)