###### GO term enrichment #######################################################################################
library(tidyverse)
library(patchwork) # used to patch plots together.  running 'GO_cc_up + GO_cc_down' will then render both plots side by side 
library(BiocManager)
library(biomaRt)
library(grid)
#BiocManager::install("org.At.tair.db")
#BiocManager::install("org.Hs.eg.db") # used to force the tair file to change subontology groups
library(org.Hs.eg.db)
library(org.At.tair.db)
library(DOSE)
library(clusterProfiler)


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
# if you want all differentially expressed genes in one table; both up and down reg together

# sig_WER_all_info_filter <- sig_WER_all_info %>% 
#   filter(abs(log2FoldChange) > 1.3, (WER_0As_FPKM > 8 | WER_50As_FPKM > 8))
# sig_Cortex_all_info_filter <- sig_Cortex_all_info %>% 
#   filter(abs(log2FoldChange) > 1.3, (Cortex_0As_FPKM > 8 | Cortex_50As_FPKM > 8))
# sig_SCR_all_info_filter <- sig_SCR_all_info %>% 
#   filter(abs(log2FoldChange) > 1.3, (SCR_0As_FPKM > 8 | SCR_50As_FPKM > 8))

# if you want up and down reg diff express genes seperate

sig_WER_down <- sig_WER_all_info %>% 
  filter(log2FoldChange < -1.3, WER_0As_FPKM > 8)
sig_Cortex_down <- sig_Cortex_all_info %>% 
  filter(log2FoldChange < -1.3, Cortex_0As_FPKM > 8)
sig_SCR_down <- sig_SCR_all_info %>% 
  filter(log2FoldChange < -1.3, SCR_0As_FPKM > 8)

sig_WER_up <- sig_WER_all_info %>% 
  filter(log2FoldChange > 1.3, WER_50As_FPKM > 8)
sig_Cortex_up <- sig_Cortex_all_info %>% 
  filter(log2FoldChange > 1.3, Cortex_50As_FPKM > 8)
sig_SCR_up <- sig_SCR_all_info %>% 
  filter(log2FoldChange > 1.3, SCR_50As_FPKM > 8)

# summary table of sig up/down genes per cell type
GO_sum_table<-tibble('Cell-Type'=c('Epidermal','Cortex','Endodermal'),'Up-Reg'=c(2564,2459,2145),'Down-Reg'=c(2465,2042,2454))
grid.table(GO_sum_table, rows = NULL)


# Make enrichResult objects for up/down, cell-type, & subontology categories (Biological process = BP, Molecular Function = MF, Cellular Components = CC)
# GO_WER_down_bp, GO_WER_down_mf, GO_WER_down_cc, GO_WER_up_bp, GO_WER_up_mf, GO_WER_up_cc, ...18 total objects

#### BUG #### In order to change ont ('MF', 'BP', 'CC'), you need to run the example org.Hs.eg.db below, then go back to our org.At.tair.db data

GO_SCR_down_mf <- clusterProfiler::enrichGO(
  gene          = sig_SCR_down$ENTREZID,
  OrgDb         = 'org.At.tair.db',
  ont           = "MF", # BP, MF, or CC 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.01, 
  readable      = T,
  pool = F)

# use this to force the tair package to change from BP to MF or CC.  Run this with 'CC', then go back and run the tair file with 'CC'
library(org.Hs.eg.db)
data(geneList, package = "DOSE")
de <- names(geneList)[1:100]

yy <- enrichGO(de, 'org.Hs.eg.db', ont="CC", pvalueCutoff=0.01)
yy <- as_tibble(yy)
View(yy)   

### the merge_results function doesn't seem to be working in the clusterProfiler package, so we need to compile the results manually
# 1) change each enrichResult object to a tibble 
# 2) create new columns for groups (cell_type = c(Epidermal,Cortex,Endodermal), subontology = c(Biological_Process, Molecular_Function, Cellular_Components), regulation = c(Up, Down))
# 3) paste tables together 
# 4) ggplot with facets

#Cortex down
GO_Cortex_down_bp <- as_tibble(GO_Cortex_down_bp) %>% 
  mutate(cell_type = 'Cortex', Regulation = 'Down', subontology = 'Biological_Process')
GO_Cortex_down_mf <- as_tibble(GO_Cortex_down_mf) %>% 
  mutate(cell_type = 'Cortex', Regulation = 'Down', subontology = 'Molecular_Function')
GO_Cortex_down_cc <- as_tibble(GO_Cortex_down_cc) %>% 
  mutate(cell_type = 'Cortex', Regulation = 'Down', subontology = 'Cellular_Component')
#Cortex up
GO_Cortex_up_bp <- as_tibble(GO_Cortex_up_bp) %>% 
  mutate(cell_type = 'Cortex', Regulation = 'Up', subontology = 'Biological_Process')
GO_Cortex_up_mf <- as_tibble(GO_Cortex_up_mf) %>% 
  mutate(cell_type = 'Cortex', Regulation = 'Up', subontology = 'Molecular_Function')
GO_Cortex_up_cc <- as_tibble(GO_Cortex_up_cc) %>% 
  mutate(cell_type = 'Cortex', Regulation = 'Up', subontology = 'Cellular_Component')

#WER down
GO_WER_down_bp <- as_tibble(GO_WER_down_bp) %>% 
  mutate(cell_type = 'Epidermal', Regulation = 'Down', subontology = 'Biological_Process')
GO_WER_down_mf <- as_tibble(GO_WER_down_mf) %>% 
  mutate(cell_type = 'Epidermal', Regulation = 'Down', subontology = 'Molecular_Function')
GO_WER_down_cc <- as_tibble(GO_WER_down_cc) %>% 
  mutate(cell_type = 'Epidermal', Regulation = 'Down', subontology = 'Cellular_Component')
#WER up
GO_WER_up_bp <- as_tibble(GO_WER_up_bp) %>% 
  mutate(cell_type = 'Epidermal', Regulation = 'Up', subontology = 'Biological_Process')
GO_WER_up_mf <- as_tibble(GO_WER_up_mf) %>% 
  mutate(cell_type = 'Epidermal', Regulation = 'Up', subontology = 'Molecular_Function')
GO_WER_up_cc <- as_tibble(GO_WER_up_cc) %>% 
  mutate(cell_type = 'Epidermal', Regulation = 'Up', subontology = 'Cellular_Component')

#SCR down
GO_SCR_down_bp <- as_tibble(GO_SCR_down_bp) %>% 
  mutate(cell_type = 'Endodermal', Regulation = 'Down', subontology = 'Biological_Process')
GO_SCR_down_mf <- as_tibble(GO_SCR_down_mf) %>% 
  mutate(cell_type = 'Endodermal', Regulation = 'Down', subontology = 'Molecular_Function')
GO_SCR_down_cc <- as_tibble(GO_SCR_down_cc) %>% 
  mutate(cell_type = 'Endodermal', Regulation = 'Down', subontology = 'Cellular_Component')
#SCR up
GO_SCR_up_bp <- as_tibble(GO_SCR_up_bp) %>% 
  mutate(cell_type = 'Endodermal', Regulation = 'Up', subontology = 'Biological_Process')
GO_SCR_up_mf <- as_tibble(GO_SCR_up_mf) %>% 
  mutate(cell_type = 'Endodermal', Regulation = 'Up', subontology = 'Molecular_Function')
GO_SCR_up_cc <- as_tibble(GO_SCR_up_cc) %>% 
  mutate(cell_type = 'Endodermal', Regulation = 'Up', subontology = 'Cellular_Component')

#bind all tibbles together to make 1 large dataset, .id = 'dataset' creates a column that tells you which dataset each row comes from
full_GO <- bind_rows(GO_Cortex_down_bp, GO_Cortex_down_cc, GO_Cortex_down_mf,
                     GO_Cortex_up_bp,   GO_Cortex_up_cc,   GO_Cortex_up_mf, 
                     GO_WER_down_bp,    GO_WER_down_cc,    GO_WER_down_mf, 
                     GO_WER_up_bp,      GO_WER_up_cc,      GO_WER_up_mf, 
                     GO_SCR_down_bp,    GO_SCR_down_cc,    GO_SCR_down_mf, 
                     GO_SCR_up_bp,      GO_SCR_up_cc,      GO_SCR_up_mf, 
                     .id = "dataset")

# need this function to make a reverse log scale x-axis
library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

GO_cc_up <- ggplot(
              filter(full_GO, subontology == 'Cellular_Component', Regulation == 'Up'), 
                    aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
                geom_point(aes(size = Count, color = cell_type)) +
                theme_dark(base_size = 12) +
                ylab(NULL) +
                scale_x_continuous(trans=reverselog_trans(10) ) +
                guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                  panel.grid.minor = element_blank()) +
                labs(title = "Cellular Component",
                  subtitle = "Up-regulated Genes")

GO_cc_down <- ggplot(
                filter(full_GO, subontology == 'Cellular_Component', Regulation == 'Down'), 
                     aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
                  geom_point(aes(size = Count, color = cell_type), position = position_jitter(height = 0.25)) +
                  theme_dark(base_size = 12) +
                  ylab(NULL) +
                  scale_x_continuous(trans=reverselog_trans(10)) +
                  scale_y_discrete(position = 'right') +
                  guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
                  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                    panel.grid.minor = element_blank(),
                    legend.position = "none") +
                  labs(title = "Cellular Component",
                    subtitle = "Down-regulated Genes") 

GO_cc_up+GO_cc_down
#ggsave('GO_cc.tiff', plot=last_plot(), device = 'tiff', dpi = 'retina')

GO_mf_up <- ggplot(
              filter(full_GO, subontology == 'Molecular_Function', Regulation == 'Up'), 
                   aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
                geom_point(aes(size = Count, color = cell_type),position = position_jitter(height = 0.15)) +
                theme_dark(base_size = 12) +
                ylab(NULL) +
                scale_x_continuous(trans=reverselog_trans(10) ) +
                guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                  panel.grid.minor = element_blank()) +
                  #legend.position = "none") +
                labs(title = "Molecular Function",
                  subtitle = "Up-regulated Genes")

GO_mf_down <- ggplot(
                filter(full_GO, subontology == 'Molecular_Function', Regulation == 'Down'), 
                     aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
                  geom_point(aes(size = Count, color = cell_type), position = position_jitter(height = 0.15)) +
                  theme_dark(base_size = 12) +
                  ylab(NULL) +
                  scale_x_continuous(trans=reverselog_trans(10)) +
                  scale_y_discrete(position = 'right') +
                  guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
                  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                    panel.grid.minor = element_blank()) +
                    #legend.position = "none") +
                  labs(title = "Molecular Function",
                    subtitle = "Down-regulated Genes") 

GO_mf_up+GO_mf_down
#ggsave('GO_mf.tiff', plot=last_plot(), device = 'tiff', dpi = 'retina')

GO_bp_up <- ggplot(
              filter(full_GO, subontology == 'Biological_Process', Regulation == 'Up'), 
                   aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
                geom_point(aes(size = Count, color = cell_type),position = position_jitter(height = 0.15)) +
                theme_dark(base_size = 12) +
                ylab(NULL) +
                scale_x_continuous(trans=reverselog_trans(10) ) +
                guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                  panel.grid.minor = element_blank()) +
                  #legend.position = "none") +
                labs(title = "Biological Process",
                  subtitle = "Up-regulated Genes")

GO_bp_down <- ggplot(
                filter(full_GO, subontology == 'Biological_Process', Regulation == 'Down'), 
                     aes(x = p.adjust, y = fct_reorder(Description, Count))) + 
                  geom_point(aes(size = Count, color = cell_type), position = position_jitter(height = 0.15)) +
                  theme_dark(base_size = 12) +
                  ylab(NULL) +
                  scale_x_continuous(trans=reverselog_trans(10)) +
                  scale_y_discrete(position = 'right') +
                  guides(color = guide_legend(override.aes = list(size=5))) + # increases the color legend size
                  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                    panel.grid.minor = element_blank()) +
                    #legend.position = "none") +
                  labs(title = "Biological Process",
                    subtitle = "Down-regulated Genes") 

# need higher resolution than ggsave "retina"
tiff("GO_bp.tiff", units="in", width=15, height=22, res=900)
GO_bp_up+GO_bp_down
dev.off()         
