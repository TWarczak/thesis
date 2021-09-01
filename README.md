## My thesis **"Characterizing genetic factors that determine arsenic tolerance in plants"** is available online at [ProQuest](https://www.proquest.com/openview/70234b9ebe5923fa1a93e889d61fabe0/1?pq-origsite=gscholar&cbl=18750&diss=y). 

## Download full thesis for free here [Warczak, 2020](https://raw.githubusercontent.com/TWarczak/thesis/master/TWarczak_thesis.pdf)

### *Will post more code after publications*

Chapter 1 is an introduction to arsenic as a non-essential, toxic metalloid that poses serious risks towards crop yields and human health. I discuss sources of dietary arsenic and solutions for mitigating consumption. 

Chapter 2 covers my RNA-seq pipeline for gene expression of 25000+ plant (*Arabidopsis thaliana*) genes. This includes cell-type specific gene expression, clustering (PCA, hierarchical clustering), regression (GLMs/ANOVA), exploratory data analysis, and causal inference techniques. 

Chapter 3 covers how I engineered a novel genome-wide association study (GWAS) that identified genes controlling arsenic tolerance in plant roots. I leveraged expert data cleaning and wrangling in R to summarise findings from millions of observations across thousounds of plant genomes. Modeling was performed in C/C++ on Ubuntu virtual machine and all analysis in R. I ultimately determined the plant gene AtNIP1;1 is the major genetic factor for tolerating arsenic in *A. thaliana*  root cells and identified regions of interest on multiple chromosomes for future studies. 

### Additional publications:

Nachman et al., 2018. **Opportunities and Challenges for Dietary Arsenic Intervention**. *Environ Health Perspect* Aug; 126(8): 084503.
[Link](https://raw.githubusercontent.com/TWarczak/thesis/master/Nachman_et_al_2018.pdf)

Punshon et al., 2017. **Understanding arsenic dynamics in agronomic systems to predict and prevent uptake by crop plants**. *Sci Total Environ* **581-582**: 209-220
[Link](https://raw.githubusercontent.com/TWarczak/thesis/master/Punshon_et_al_2017.pdf)

## For now, here are the notable figures from the thesis, chapters 1-3.  

***

![.Figure1.1/figure1.1.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Figure1.1/figure1.1.png)

**Figure 1.1: Generalized Diagram of Arsenic Transport and Metabolism in Plants**

Modified from Punshon et al. (2017), Zhao et al. (2009), and Ma et al. (2007). 
***

**Notable Figures for Chapter 2 of Thesis**

**Chapter 2: Cell-type specific expression profiles in Arabidopsis thaliana roots in response to AsIII**

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
***
***

**Notable Figures for Chapter 3 of Thesis**

**Chapter 3: Natural variation in NIP1;1 expression determines AsIII tolerance in A. thaliana accessions**

***

![.Ch_3_figures/figure3.1.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/figure3.1.png)

**Figure 3.1: MAGIC Line Root Scoring for Arsenic Tolerance**

(A) Roots for MAGIC lines #193 and #402 scored for arsenic tolerance. (B) Distribution of 492 MAGIC line root scores used for QTL mapping.  MAGIC lines which tolerate AsIII produced negative root scores; MAGIC lines sensitive to AsIII produce positive root scores.  Min = -1.94, 1st Qu. = -0.32, Median = 0.33, Mean = 0.30, 3rd Qu. = 0.85, Max = 3.18.  Red dots in (A) indicate where the main root tip of each plant ended and thus where each measurement occurred.

***

![.Ch_3_figures/figure3.2A.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/figure3.2A.png)
![.Ch_3_figures/figure3.2B_C.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/figure3.2B_C.png)

**Figure 3.2: GWAS for Arsenic Tolerance Found in MAGIC Population**

(A) Manhattan plot showing association of ~3.3 million sequence variants with 392 root scores for arsenic tolerance within the MAGIC population, mapped onto the five A. thaliana chromosomes.  Horizontal lines represent conventional statistical significance thresholds of -log10(1e-5) (green dash) and -log10(5*10^-8) (blue dash).  Red line represents the p < 0.05 threshold after Bonferroni correction, -log10(0.05/3.3M).  (B) Magnified view of major QTL on chromosome four with top seven variants associated with AsIII tolerance labelled.  (C) Location of 2nd, 3rd, and 4th most significant variants (rs2412782, and rs2412888, rs2412612) compared to NIP1;1 locus and nearby genes.  Segments: blue = UTR, pink = exon, green = intron. 

***

![.Ch_3_figures/figure3.3A_B.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/figure3.3A_B.png)
![.Ch_3_figures/figure3.3C_D.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/figure3.3C_D.png)


**Figure 3.3: MAGIC Line Haplotypes at Major and Minor QTL**

Boxplots of MAGIC line haplotypes at the NIP1;1 (A) and Chr3 QTL (C) loci showing relationship between arsenic sensitivity phenotype rank (Parental Rank), arranged from most sensitive (Can) to most tolerant (Edi), and arsenic root scores.  Each box represents the 25%, 50%, and 75% quantiles.  Red dotted line represents root score of 0.  Correlation of the average root scores for each MAGIC line haplotype found at the NIP1;1 (B) and Chr3 QTL (D) loci and phenotype rank of the 19 MAGIC line founders.  NIP1;1 locus correlates with Rank (R2 = 0.44).  Chr3 QTL locus does not correlate with Rank (R2 = 0.1).  

***

![.Ch_3_figures/figure3.4A_B_C_D.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/figure3.4A_B_C_D.png)
![.Ch_3_figures/figure3.4E_F.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/figure3.4E_F.png)

**Figure 3.4: MAGIC founder NIP1;1 expression**

(A) Founder NIP1;1 expression (qPCR), relative to Col-0 at 0dAs.  Samples 0dAs – 5dAs represent plant roots harvested either after 7d on control medium (0dAs) or after 1-5 days transferred to 10 µM AsIII media.  Founders are arranged by AsIII tolerance phenotype rank, most sensitive (Can) to most tolerant (Edi) to show correlation of phenotype with NIP1;1 expression.  (B) Founder NIP1;1 expression (qPCR) at only 0dAs (7d old plants).  Statistical significance was calculated using one-way ANOVA (*P < 0.05) with Tukey’s Test.  (C) Strong correlation between cumulative NIP1;1 expression (0dAs – 5dAs) and phenotype rank of the 19 MAGIC line founders.  (D) Strong correlation between NIP1;1 expression (0dAs only) and phenotype rank.  (E) Phenotype of transgenic lines (C1, C2, C3) expressing NIP1;1 from Can-0 with native promoter in the nip1;1 background.  (F) NIP1;1 expression, relative to Col-0, from root tissue of 7d plants on control medium.    

***

![.Ch_3_figures/sup_figure3.1.1.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/sup_figure3.1.1.png)
![.Ch_3_figures/sup_figure3.1.2.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/sup_figure3.1.2.png)

**Supplemental Figure 3.1: MAGIC Founder Lines; ordered by AsIII sensitivity rank**

MAGIC founder lines grown on control medium (half-strength MS) for 11 days, and 10 μM AsIII medium for 15 days. Founders are arranged from most AsIII sensitive (Can-0) to tolerant (Edi-0). Order used as AsIII sensitivity rank (1-19): Can-0, Oy-0, Po-0, Rsch-4, Zu-0, Tsu-0, Wu-0, Ct-1, Wil-2, Ws-0, Bur-0, Col-0, Ler-0, Mt-0, No-0, Hi-0, Sf-2, Kn-0, Edi-0.

***

![.Ch_3_figures/sup_figure3.2.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/sup_figure3.2.png)

**Supplemental Figure 3.2: Variants of Interest**

(A) Top 7 variants at minor QTL on chromosome 3.  (B) Table summary of variants and position on chromosome 3.

***

![.Ch_3_figures/sup_figure3.3.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/sup_figure3.3.png)

**Supplemental Figure 3.3: GWAS with Condensed SNP Library Showing Association with AsIII Tolerance**

(A) Genome scan performed with the 492 root scores showing association with the 1260 SNP markers used to generate condensed genome mosaics of MAGIC lines.  Two minor QTLs found on chromosome 3 and one on chromosome 4.  One major QTL found on chromosome 4.  A 2nd major QTL might be present upstream of the main peak.  Orange dots signify peak at minor QTL, green triangle signifies major QTL, blue asterisk signifies potential QTL.  -logP of 3.51 corresponds to a genome-wide p-value of <0.05.  (B)  Allele effect estimates for the 19 MAGIC founders at the peak SNP (MN4_10482087) of the major QTL on chromosome 4. MN4_10482087 is located 59569 bp upstream of NIP1;1 but is the closest SNP to NIP1;1.  The next closest SNPs to NIP1;1 are MASC02548 and MASC01526, located 120628 bp downstream and 190423 bp upstream of the NIP1;1 5’UTR, respectively.  Higher effect estimates associate with higher root scores and AsIII sensitivity, lower effect estimates associate with lower root scores and AsIII tolerance.

***

![.Ch_3_figures/sup_figure3.4.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/sup_figure3.4.png)

**Supplemental Figure 3.4: Root Scores of Minor & Major QTL Haplotype Combinations**

Heatmap detailing average root score of every MAGIC line founder haplotype combination at minor (QTL3) and major (QTL4) QTLs.  Combinations show an uneven distribution of the two QTL loci.  Many haplotype combinations are not present in the MAGIC population.  For example, Ler-0 haplotype at QTL4 is well represented and found in combination with 17 of 19 other haplotypes at QTL3.  The Po-0 haplotype at QTL4, in contrast, is underrepresented and only found in combination with 5 of 19 other haplotypes at QTL3, making it difficult to assess the impact of Po-0 NIP1;1 on root scores.

***

![.Ch_3_figures/sup_figure3.5.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/sup_figure3.5.png)

**Supplemental Figure 3.5: Counts of MAGIC line Founders Per Genome**

Histogram counts the number of MAGIC founder haplotypes contributing to each MAGIC line of the 391 MAGIC line genomes used in Figure 2 GWAS.  Compare to Figure 1 in Korver et al. (2009).  

***

![.Ch_3_figures/figure3.omit_1.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/figure3.omit_1.png)

**Not Included in Thesis**

***

![.Ch_3_figures/figure3.omit_2.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/figure3.omit_2.png)

**Not Included in Thesis**

***

![.Ch_3_figures/figure3.omit_3.png](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/figure3.omit_3.png)

**Not Included in Thesis**

***

![.Ch_3_figures/drake_meme.tif](https://raw.githubusercontent.com/TWarczak/thesis/master/Ch_3_figures/drake_meme.tif)

**Not Included in Thesis**

***

