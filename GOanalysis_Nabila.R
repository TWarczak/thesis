##all genes in file and converting to dataframe
GenesWTSig<- as.data.frame(WTSig) %>%
  rownames_to_column("gene_id") 
GenesWTSig
View(GenesWTSig)
summary( GenesWTSig)
##GenesWTSig here is my datarframe containing genes DE with log2fold and pvalues in the dataframe
## GENE ANNOTATION USING TAIR WEBSITE

library(AnnotationHub)
library(AnnotationDbi)
library(org.At.tair.db) # might need to click from Packages tab 
library(EnrichmentBrowser)
BiocManager::install('EnrichmentBrowser')
columns(org.At.tair.db)
n
#id_list <- GenesWTSig[,1] ## taking in column with gene ids
WER_id_list <- sig_WER[,1] ## taking in column with gene ids

#id_list
WER_id_list

#id_list <- as.character(id_list)
WER_id_list <- as.character(WER_id_list)
WER_id_list

#symbol_list <- mapIds(org.At.tair.db, keys= id_list, column=c("SYMBOL"), keytype="TAIR", multiVals="first")
symbol_list <- mapIds(org.At.tair.db, keys= WER_id_list, column=c("SYMBOL"), keytype="TAIR", multiVals="first")

#ENTREZID <- mapIds(org.At.tair.db, keys= id_list, column=c("ENTREZID"), keytype="TAIR", multiVals="first")
ENTREZID <- mapIds(org.At.tair.db, keys= WER_id_list, column=c("ENTREZID"), keytype="TAIR", multiVals="first")
ENTREZID

ENTREZID <- as.data.frame(ENTREZID)
ENTREZID
description_list <- mapIds(org.At.tair.db, keys= id_list, column=c("GENENAME"), keytype="TAIR", multiVals="first")

descriptions <- as.data.frame(description_list)

finalannoated_WT <- cbind(GenesWTSig, symbol_list, descriptions,ENTREZID)

View (finalannoated_WT)
rownames(finalannoated_WT)
as.data.frame(finalannoated_WT) %>%
  write.table("WT genes_LOG2FOLD_1")  # Saving all genes differentially expressed in WT

##GO analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DOSE")

library(DOSE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")


library(clusterProfiler)


OrgDb <-org.At.tair.db ##will need to use this later for GO analysis

##below just to check whats in the orgdb and what can we use for analysis
columns(org.At.tair.db)
keytypes(org.At.tair.db)
levels(org.At.tair.db)



geneList <- as.vector(finalannoated_WT$log2FoldChange)
View(geneList)

names(geneList) <- finalannoated_WT$ENTREZID ##will use GO based on ENTREZ ID

names(geneList) 

class(names(geneList))

ego <- clusterProfiler::enrichGO(gene          = names(geneList),
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
?groupGO
ggo <- clusterProfiler::groupGO(gene     = names(geneList),
                                OrgDb    = OrgDb,
                                ont      = "BP",
                                level    = 3,
                                readable = TRUE)
ggo


plot(ggo)

barplot(ggo, drop=TRUE, showCategory=40)
