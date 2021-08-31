WER_FPKM <- FPKM_avg %>% 
  filter(gene_id %in% c('AT5G14750')) 

x <- read.delim("clipboard")
x$AGI <- toupper(as.character(x$AGI)) 
# remove the following
y= c('AT1G15270', 'AT1G20696', 'AT1G27090', 'AT1G27310', 'AT1G31070', 'AT1G49600', 'AT1G52360', 'AT1G76010', 'AT1G80410', 'AT2G07360', 
     'AT2G20490', 'AT2G20750','AT2G31410', 'AT2G37600', 'AT2G45080','AT3G01390', 'AT3G03380', 'AT3G07050', 'AT3G07660', 'AT3G13460', 
     'AT3G15357', 'AT3G18060', 'AT3G20630','AT3G21350', 'AT3G49080', 'AT3G50140', 'AT3G52640', 'AT3G54580', 'AT3G60280', 'AT3G61110', 
     'AT3G63460', 'AT4G15410', 'AT4G16143', 'AT4G27640', 'AT4G32030', 'AT4G34290','AT4G34870', 'AT4G38600', 'AT5G13030', 'AT5G23140', 
     'AT5G24690', 'AT5G27430', 'AT5G35430', 'AT5G39850', 'AT5G53480', 'AT5G56320', 'AT5G56950', 'AT5G59840') 
Epiderm_genes <- FPKM_avg %>% 
  filter(gene_id %in% x$AGI, !(gene_id %in% y)) %>%
  bind_rows(WER_FPKM) %>% 
  modify_if(~is.numeric(.), ~round(., 0))

View(Epiderm_genes)


z <- read_delim(clipboard(), col_names = F, delim = ' ') %>% 
  dplyr::rename(AGI = X1) 
z$AGI <- toupper(as.character(z$AGI))
View(z)
Cortex_FPKM <- FPKM_avg %>% 
  filter(gene_id %in% c('AT1G09750'))
y = c('AT1G22740','AT1G51940','AT1G65060','AT3G10525','AT3G20660','AT3G56200','AT5G07190','AT5G42100','AT5G62210','AT5G66310')

Cortex_genes <- FPKM_avg %>% 
  filter(gene_id %in% z$AGI, Cortex_minus_As > 2, !(gene_id %in% y)) %>%
  modify_if(~is.numeric(.), ~round(., 0))
View(Cortex_genes)


SCR_FPKM <- FPKM_avg %>% 
  filter(gene_id %in% c('AT3G54220'))
w <- read_delim(clipboard(), col_names = F, delim = ' ') %>% 
  dplyr::rename(AGI = X1) 
w <- w %>% 
  add_row(AGI = c('At4g14010','At1g63520','At5g14510','AT4G21340'))
w$AGI <- toupper(as.character(w$AGI))
Endoderm_genes <- FPKM_avg %>% 
  filter(gene_id %in% w$AGI, SCR_minus_As > 2, !(gene_id %in% y)) %>%
  bind_rows(SCR_FPKM) %>% 
  modify_if(~is.numeric(.), ~round(., 0))
View(Endoderm_genes)
y= c('AT1G66230','AT2G30400','AT2G32120','AT3G10340','AT3G18280','AT4G27350','AT4G33550','AT4G39675','AT5G46590')

cell_type_sp_genes <- bind_rows(Epiderm_genes, Cortex_genes, Endoderm_genes)
write.table(cell_type_sp_genes, 'cell_type_sp_genes.csv', row.names = F)
cell_type_sp_genes <- read_delim("~/R/DESeq2_Todd_RNAseq/cell_type_sp_genes.csv", 
                                  " ", escape_double = FALSE, trim_ws = TRUE)

cell_type_sp_genes <- cell_type_sp_genes %>% 
  select(gene_id, WER_minus_As, Cortex_minus_As, SCR_minus_As)

heatmap_cell_type_sp_mat <- cell_type_sp_genes %>% 
  column_to_rownames(var = 'gene_id') %>%
  as.matrix()  

# function to generate parameter for scaling each row; each FPKM observation - mean(FPKM) across row, div by sd 
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
heatmap_cell_type_sp_norm <- t(apply(heatmap_cell_type_sp_mat, 1, cal_z_score))
row.names(my_sample_col) <- colnames(heatmap_cell_type_sp_mat)

my_hclust_gene <- hclust(dist(heatmap_cell_type_sp_norm), method = "ward.D2")
as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)

my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 3) 
my_gene_col <- as.data.frame(my_gene_col)
my_gene_col <- dplyr::rename(my_gene_col, cluster = my_gene_col)

my_sample_col <- data.frame(condition = c('0As', '0As', '0As'))

mypalette <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)

my_color = list(cluster = c(cluster1='#6a3d9a',cluster2='#ff7f00',cluster3='#b2df8a'))#cluster4='#33a02c',cluster5='#1a1a1a',
                            # cluster6='#e31a1c',cluster7='#fdbf6f',cluster8='#ff7f00',cluster9='#cab2d6',
                            # cluster10='#6a3d9a',cluster11='#ffff99',cluster12='#b15928',cluster13='#878787'))

tiff("heatmap_cell_type_sp.tiff", units="in", width=7, height=10, res=900)
pheatmap(heatmap_cell_type_sp_norm, 
         annotation_row = my_gene_col, 
         color = mypalette, 
         annotation_colors = my_color,
         annotation_legend = F,
         cellwidth = 45, 
         clustering_method = 'ward.D2',
         cluster_cols = F,
         show_rownames=T,
         cutree_rows = 3,
         angle_col = 45,
         labels_col = c('Epidermal_0As', 'Cortex_0As', 'Endodermal_0As'))

dev.off()         


