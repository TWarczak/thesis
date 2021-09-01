# thesis
thesis-work figures 





***

![.Figure1.1/figure1.1.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Figure1.1/figure1.1.png)

**Figure 1.1: Generalized Diagram of Arsenic Transport and Metabolism in Plants**

Modified from Punshon et al. (2017), Zhao et al. (2009), and Ma et al. (2007). 

***

![.Ch_2_figures/figure2.1A_B_C.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/figure2.1A_B_C.png)

![.Ch_2_figures/figure2.1D.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/figure2.1D.png)

**Figure 2.1: Unsupervised Clustering & Exploratory Data Analysis**

(A) Principal Component Analysis (PCA) plot showing the samples clustering based on the two principal components, cell-type and condition, which explains 14% and 79% of the variation, respectively, in the RNA-seq data.  (B) Hierarchical clustering heatmap showing similarities in global gene expression between both the biological replicates and the different sample groups.  Perfect correlation = 1.  (C) Heatmap of genes in control conditions that show cell-type specific expression in Dinneney et al. 2008 and this study.  The promoters of At5g14750, At3g54220, and At1g09750 were used to drive expression of GFP in the WER::GFP,  SCR::GFP, and Cortex::GFP marker lines.  Color scale based on difference between RPKM of gene in cell-type and mean RPKM of gene in all samples, divided by the standard deviation (x - mean(x)) / sd(x). (D) Volcano plots showing the fold changes relative to the adjusted p-values for all genes in epidermal, cortex, and endodermal cell-types.  Blue dots are statistically significant differentially expressed genes padj < 0.05, red dots are not significant.  

***

![.Ch_2_figures/figure2.2.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/figure2.2.png)

**Figure 2.2: Differentially Expressed Genes by Cell-Type**

(A) Hierarchical clustering of 7543 significant differentially expressed genes with abs(logFC) > 1 and RPKM > 6 in at least one cell-type of any condition.  8 clusters can be seen showing relative DEG expression in each cell-type and condition.  Color scale the difference between RPKM of gene in cell-type and mean RPKM of gene in all samples, divided by the standard deviation (x - mean(x)) / sd(x).  (B) Venn diagram of same 7543 DEGs differentially expressed in unique or combinations of cell-types.  (C) Density plots of the up (bottom row) and down (top row) DEG RPKM values in epidermal, cortex, or endodermal cells in either control (blue) or AsIII (red) conditions.  

***

![.Ch_2_figures/figure2.3A_B.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/figure2.3A_B.png)

![.Ch_2_figures/figure2.3C.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/figure2.3C.png)

**Figure 2.3: GO-Term Enrichment**

Molecular Function (A), Cellular Component (B), and Biological Process (C) subontology terms enriched in response to AsIII.  Terms for each subontology arranged by adjusted p-values for up and down differentially expressed genes.  Size corresponds to the number of genes (count) enriched for that term for each cell-type. Differentially expressed genes were filtered to include only abs(logFC) > 1.3, and RPKM > 8 in AsIII condition for up-regulated genes or RPKM > 8 in control condition for down-regulated genes. Cell-types separated by color. 

***

![.Ch_2_figures/figure2.4.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/figure2.4.png)

**Figure 2.4: AsIII Induced Expression by Cell-Type**

Expression heatmaps showing LogFC by cell-type in response to AsIII on left and the corresponding RPKM units arranged by cell-type and condition for each gene on right. Genes from the Major Intrinsic Proteins (MIP) superfamily of aquaporins (A), ATP-Binding Cassette (ABC) family (B), and Multidrug and Toxic Compound Extrusion (MATE) family (C).  Genes that lacked a RPKM value of > 2 in at least one sample were filtered out. Heatmaps clustered by logFC.

***

![.Ch_2_figures/sup_figure2.1.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_figure2.1.png)

**Supplemental Figure 2.1: MA Plots for Differentially Expressed Genes**

MA plots (log fold-change versus mean of normalized read counts) showing statistically significant, differentially expressed genes (red dots) in epidermal (A), cortex (B), and endodermal (C) cell-types.

***

![.Ch_2_figures/sup_figure2.2A_B_C_D_E.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_figure2.2A_B_C_D_E.png)

![.Ch_2_figures/sup_figure2.2F.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_figure2.2F.png)

**Supplemental Figure 2.2: GFP-Reporter Lines**

7d old GFP reporter lines after 24 hours exposure to 50 μM AsIII show GFP fluorescence exclusive to epidermal cells (A) in WER::GFP, to endodermal cells (B) in SCR::GFP, and to cortex cells (D-E) in Cortex::GFP.  GFP fluorescence of Cortex::GFP line in control conditions (C). Scale in each image is 50 μm.  (F) Bar plot of RPKM units for the three GFP-reporter line genes in both experimental conditions.  Promoters of AT5G14750, AT3G54220, AT1G09750 were cloned for WER::GFP, SCR::GFP, and Cortex::GFP, respectively.  Error bars represent standard error from n=3 biological replicates.  

***

![.Ch_2_figures/sup_figure2.3A_B.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_figure2.3A_B.png)

![.Ch_2_figures/sup_figure2.3C_D.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_figure2.3C_D.png)

![.Ch_2_figures/sup_figure2.3E_F.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_figure2.3E_F.png)

![.Ch_2_figures/sup_figure2.3G.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_figure2.3G.png)

**Supplemental Figure 2.3: GFP-marker lines used for FACS**

Examples of complete gating parameters and sorting statistics used in FACS for 7d old WT and GFP-reporter lines after 24 hours exposure to control or 50 μM AsIII conditions then root protoplasting.  Protoplasts filtered through gates P1-3 for size and shape that fell within gate P5 were collected as GFP positive events.  WT control protoplasts (A), WER::GFP control protoplasts (B), WER::GFP AsIII protoplasts (C), Cortex::GFP control protoplasts (D), Cortex::GFP AsIII protoplasts (E), SCR::GFP control protoplasts (F), SCR::GFP AsIII protoplasts (G). 

***

![.Ch_2_figures/sup_figure2.5.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_figure2.5.png)

**Supplemental Figure 2.5: ICP-MS by Cell-Type**

Bar plots of total ng arsenic (As), sulfur (S), or silicon (Si) per 50,000 sorted cortex, endodermal, or epidermal cells in both experimental conditions.

***

![.Ch_2_figures/sup_figure2.6.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_figure2.6.png)

**Supplemental Figure 2.6: AsIII Induced Expression by Cell-Type**

Expression heatmaps showing LogFC by cell-type in response to AsIII on left and the corresponding RPKM units arranged by cell-type and condition for each gene on right. Genes from the WRKY family transcription factors (A), MYB family transcription factors (B), and other genes of interest (C). Heatmaps include genes with logFC > 1 and a RPKM value of > 15 in at least one sample. Heatmaps clustered by logFC.

***

![.Ch_2_figures/sup_figure2.7.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_figure2.7.png)

**Supplemental Figure 2.7: A. thaliana GO Term Limitations**

Violin plots of Biological Process (red), Cellular Component (green), and Molecular Function (blue) GO terms.  Data points are unique GO terms with current (May 2020) counts of A. thaliana genes annotated with that term.  Plots highlight some limitations to typical GO enrichment analysis.  Terms with < 10 gene counts are too specific (eRF1 methyltransferase complex) and terms with ~ 10,000 gene counts are too vague (molecular_function).  Many terms such as response to cadmium ion are over-represented (Gene Count = 313) while many terms like mitotic cell cycle are under-represented (Gene Count = 55).

***

![.Ch_2_figures/sup_table2.1.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_table2.1.png)

**Supplemental Table 2.1: RNA Library Counts & Alignment**

24 RNA samples used in this study were sequenced, trimmed, filtered, and aligned to the A. thaliana reference genome (Araport11).  Sequencing was performed using the Illumina Nextseq500 system with 75bp fragment single-read libraries. Mapped read depth of samples ranged from ~ 11,000,000 – 33,000,000 per sample. 

***

![.Ch_2_figures/sup_table2.2.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_table2.2.png)

**Supplemental Table 2.2: Gene Expression (RPKM) Counts from Up-Regulated DEGs by Cell-Type and Condition**

Numbers in the table represent the number of up-regulated DEGs within the RPKM Range bin and the corresponding percentage in that cell-type and condition.

***

![.Ch_2_figures/sup_table2.3.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_2_figures/sup_table2.3.png)

**Supplemental Table 2.3: Gene Expression (RPKM) Counts from Down-Regulated DEGs by Cell-Type and Condition**

Numbers in the table represent the number of down-regulated DEGs within the RPKM Range bin and the corresponding percentage in that cell-type and condition. 

***






