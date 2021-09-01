library(tidyverse)

# MAGIC_genotypes_byparent <- read_tsv("http://mtweb.cs.ucl.ac.uk/mus/www/POOLING/ARABIDOPSIS/FOUNDER/GENOTYPES/mosaic10.lowcov.txt", col_names = T)
MAGIC_genotypes_byparent <- read_delim("~/R/MAGIC_Analysis/GWAS/Sept_2019/mosaic.txt", skip = 1, delim = ' ', 
                                       col_names = c('magic',	'chr', 'acc'	,'from.bp'	,'to.bp'	,'from.site'	,
                                                     'to.site'	,'len.bp'	,'sites'	,'errors'	,'error.site',	'error.bp'))
# Only want MAGIC lines 1-527
MAGIC_genotypes_byparent <- MAGIC_genotypes_byparent %>% 
  dplyr::slice(985:24652) %>% 
  separate(col = magic, into = c('Magic', 'magic_line')) %>% 
  select(-Magic)
  
MAGIC_genotypes_byparent[,c('magic_line','chr','from.bp','to.bp','from.site','to.site','len.bp','sites','errors')] <- lapply(
            MAGIC_genotypes_byparent[,c('magic_line','chr','from.bp','to.bp','from.site','to.site','len.bp','sites','errors')], as.integer)

MAGIC_genotypes_byparent <- filter(MAGIC_genotypes_byparent, magic_line < 600)
#write_csv(MAGIC_genotypes_byparent, "haplotypes_every_MAGIC_line.csv", col_names = TRUE)
n_distinct(MAGIC_genotypes_byparent$magic_line)

nip1_chr4_loc <- MAGIC_genotypes_byparent %>% 
  filter(chr == 4)

n_distinct(nip1_chr4_loc$magic_line)
527-392
# 135 MAGIC line genotypes missing 

nip1_chr4_loc2 <- nip1_chr4_loc %>% 
  filter(!to.bp < 10421000, !from.bp > 10422800)

n_distinct(nip1_chr4_loc2$magic_line)

nip1_chr4_loc2 %>% 
  count(magic_line) %>% 
  arrange(desc(n))

#magic_scores <- read_tsv('~/R/MAGIC_Analysis/GWAS/phenotype.txt', col_names = T)
magic_scores <- read_csv('~/R/MAGIC_Analysis/GWAS/Sept_2019/Sept19_updated_figs/phenotype.csv', col_names = T)


magic_scores2 <- magic_scores %>% 
  separate(col = SUBJECT.NAME, into = c('Magic', 'magic_line')) %>% 
  dplyr::select(-Magic) %>% 
  drop_na(score)

n_distinct(magic_scores2$magic_line)
n_distinct(nip1_chr4_loc2$magic_line)

magic_scores2$magic_line <- as.integer(magic_scores2$magic_line)


nip1_chr4_loc2 <- inner_join(nip1_chr4_loc2, magic_scores2, by = "magic_line")

############### Box plots of magic lines genotype ~ score
nip1_chr4_loc2$acc <- nip1_chr4_loc2$acc %>%
  str_replace_all(c('can-0'='Can',  'oy-0'='Oy', 'po-0'='Po', 'rsch-4'='Rsch', 'zu-0'='Zu', 'tsu-0'='Tsu', 'wu-0'='Wu', 'ct-1'='Ct',
                    'wil-2'='Wil', 'ws-0'='Ws', "bur-0"= "Bur", 'col-0'='Col', 'ler-0'='Ler', 'mt-0'='Mt', 'no-0'='No', 'hi-0'='Hi',
                    'sf-2'='Sf', 'kn-0'='Kn', 'edi-0'='Edi')) %>% 
  factor(levels = c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi"))

n_distinct(nip1_chr4_loc2$magic_line)
# 384 magic lines with genotypes and root scores at NIP1 locus

ggplot(nip1_chr4_loc2, aes(x=acc, y=score)) + 
  geom_jitter(position=position_jitter(width=0.1, height=0.1), aes(colour=factor(acc)), alpha=0.9) +
  geom_boxplot(alpha = 0.5, show.legend = F, aes(fill=factor(acc)), outlier.shape = NA) + 
  theme_minimal() +
  theme(legend.position = 'none', 
        axis.text.x = element_text(size=11, color="black", angle = 90, margin = margin(t=-15, b=10),hjust = 0.2),
        axis.text.y = element_text(size=10, margin = margin(r=10, l=8)),
        axis.title.x = element_text(size = 11, face = 'bold'),
        axis.title.y = element_text(size = 11, face = 'bold'),
        panel.grid.major.x = element_blank()) +
  ylim(-2.15,3.3) +
  scale_y_continuous(labels=waiver()) +
  labs(x = expression(paste("\nMAGIC line haplotype @", italic('NIP1;1'),' locus')), 
       y = expression(paste("Root Score  ", italic('(neg. = tolerant, pos. = sensitive)')))) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size = 1 ) 
  
ggsave("nip1_chr4_locus_genotypes.tiff", plot = last_plot(), device = "tiff",
       dpi = 1000)


################## Correlation between root scores and ranked genotype
# need avg root score per boxplot

chr4_qtl_avg_root_score <- nip1_chr4_loc2 %>%
  group_by(acc) %>%
  summarise(avg_root_score = mean(score))  


chr4_qtl_avg_root_score$acc <- factor(chr4_qtl_avg_root_score$acc, levels = c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi"))

chr4_qtl_avg_root_score$Rank <- 1:19

formula <- y ~ x
ggplot(chr4_qtl_avg_root_score, aes(Rank, avg_root_score)) +
  geom_text(aes(label=acc)) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(size=11, angle = 45, hjust=1),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=12, face='bold'),
        axis.line = element_line(colour = "black")) +
  geom_smooth(method= 'lm', formula = formula) +
  stat_regline_equation(aes(label =  paste(..eq.label..)),
                        label.x = 13, label.y = 1.0) +
  stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~','~")), label.y = 1.1, label.x = 13) +
  scale_x_continuous(breaks = round(seq(min(chr3_qtl_avg_root_score$Rank), max(chr3_qtl_avg_root_score$Rank), by = 1),1)) +
  scale_y_continuous(breaks = round(seq(min(chr3_qtl_avg_root_score$avg_root_score), max(chr3_qtl_avg_root_score$avg_root_score), by = 0.2),1)) +
  labs(x=("Parental Rank \nSensitive-Tolerant (1-19)"),
       y=('Avg Root Score'),
       title = ('Major QTL Chr4 Haplotype Scores'))  

ggsave("Major_QTL_Chr4_Genotype_Scores_corr.tiff", plot = last_plot(), device = "tiff",
       dpi = 600)

################# QTL chr 3 plot 
chr3_qtl <- MAGIC_genotypes_byparent %>% 
  filter(chr == 3) 

chr3_qtl$magic_line <- as.integer(chr3_qtl$magic_line)

chr3_qtl <- chr3_qtl %>% 
  filter(!to.bp < 4854000, !from.bp > 4857000)

n_distinct(chr3_qtl$magic_line)

chr3_qtl %>% 
  count(magic_line) %>% 
  arrange(desc(n))

magic_scores2$magic_line <- as.integer(magic_scores2$magic_line)
chr3_qtl <- inner_join(chr3_qtl, magic_scores2, by = "magic_line")

chr3_qtl$acc <- chr3_qtl$acc %>%
  str_replace_all(c('can-0'='Can',  'oy-0'='Oy', 'po-0'='Po', 'rsch-4'='Rsch', 'zu-0'='Zu', 'tsu-0'='Tsu', 'wu-0'='Wu', 'ct-1'='Ct',
                    'wil-2'='Wil', 'ws-0'='Ws', "bur-0"= "Bur", 'col-0'='Col', 'ler-0'='Ler', 'mt-0'='Mt', 'no-0'='No', 'hi-0'='Hi',
                    'sf-2'='Sf', 'kn-0'='Kn', 'edi-0'='Edi')) %>% 
  factor(levels = c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi"))

ggplot(chr3_qtl, aes(x=acc, y=score)) + 
  geom_jitter(position=position_jitter(width=0.1, height=0.1), aes(colour=factor(acc)), alpha=0.9) +
  geom_boxplot(alpha = 0.5, show.legend = F, aes(fill=factor(acc)), outlier.shape = NA) + 
  theme_minimal() +
  theme(legend.position = 'none', 
        axis.text.x = element_text(size=11, color="black", angle = 90, margin = margin(t=-15, b=10),hjust = 0.2),
        axis.text.y = element_text(size=10, margin = margin(r=10, l=8)),
        axis.title.x = element_text(size = 11, face = 'bold'),
        axis.title.y = element_text(size = 11, face = 'bold'),
        panel.grid.major.x = element_blank()) +
  ylim(-2.15,3.3) +
  scale_y_continuous(labels=waiver()) +
  labs(x = expression(paste("\nMAGIC line haplotype @ Chr3 QTL locus")), 
       y = expression(paste("Root Score  ", italic('(neg. = tolerant, pos. = sensitive)')))) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size = 1 ) 

ggsave("chr3_qtl_genotypes.tiff", plot = last_plot(), device = "tiff",
       dpi = 1000)

################## Correlation between root scores and ranked genotype QTL3
# need avg root score per boxplot


chr3_qtl_avg_root_score <- chr3_qtl %>%
  group_by(acc) %>%
  summarise(avg_root_score = mean(score))  


chr3_qtl_avg_root_score$acc <- factor(chr3_qtl_avg_root_score$acc, levels = c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi"))

chr3_qtl_avg_root_score$Rank <- 1:19

formula <- y ~ x
ggplot(chr3_qtl_avg_root_score, aes(Rank, avg_root_score)) +
  geom_text(aes(label=acc)) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(size=11, angle = 45, hjust=1),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=12, face='bold'),
        axis.line = element_line(colour = "black")) +
  geom_smooth(method= 'lm', formula = formula) +
  stat_regline_equation(aes(label =  paste(..eq.label..)),
                        label.x = 13, label.y = 1.0) +
  stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~','~")), label.y = 1.1, label.x = 13) +
  scale_x_continuous(breaks = round(seq(min(chr3_qtl_avg_root_score$Rank), max(chr3_qtl_avg_root_score$Rank), by = 1),1)) +
  scale_y_continuous(breaks = round(seq(min(chr3_qtl_avg_root_score$avg_root_score), max(chr3_qtl_avg_root_score$avg_root_score), by = 0.2),1)) +
  labs(x=("Parental Rank \nSensitive-Tolerant (1-19)"),
       y=('Avg Root Score'),  
       title = ('Minor QTL Chr3 Haplotype Scores'))  

ggsave("Minor_QTL_Chr3_Genotype_Scores_corr.tiff", plot = last_plot(), device = "tiff",
       dpi = 600)


############### Put NIP1 locus data and Chr3 QTL data together to tease out influence on eachother
qtls_3 <- chr3_qtl %>% 
  select(magic_line, acc, from.bp, to.bp) %>% 
  rename(chr3_acc = acc, chr3_from.bp = from.bp, chr3_to.bp = to.bp)
qtls_4 <- nip1_chr4_loc2 %>% 
  select(magic_line, acc, from.bp, to.bp, score) %>% 
  rename(chr4_acc = acc, chr4_from.bp = from.bp, chr4_to.bp = to.bp)

qtls_3_4 <- left_join(qtls_3, qtls_4, by = "magic_line")
qtls_3_wu_hi_sf <- filter(qtls_3_4, chr3_acc %in% c( 'Sf'))

no_hi_wu_sf <- filter(qtls_3_4, !chr3_acc %in% c('Hi', 'Wu', 'Sf'))
no_mt_bur <- filter(qtls_3_4, !chr3_acc %in% c('Mt', 'Bur'))




ggplot(qtls_3_wu_hi_sf, aes(x=chr4_acc, y=score)) + 
  geom_jitter(position=position_jitter(width=0.1, height=0.1), aes(colour=factor(chr4_acc)), alpha=0.9) +
  geom_boxplot(alpha = 0.5, show.legend = F, aes(fill=factor(chr4_acc)), outlier.shape = NA) + 
  theme_minimal() +
  theme(legend.position = 'none', 
        axis.text.x = element_text(size=11, color="black", angle = 90, margin = margin(t=-15, b=10),hjust = 1),
        axis.text.y = element_text(size=10, margin = margin(r=10, l=8)),
        panel.grid.major.x = element_blank()) +
  ylim(-2.15,3.3) +
  scale_y_continuous(labels=waiver()) +
  labs(x = expression(paste("\nMAGIC line genotype @", italic('NIP1;1'),' locus, only Hi, Wu, Sf QTL3')), 
       y = expression(paste("Arsenic score  ", italic('(neg. = tolerant, pos. = sensitive)')))) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size = 1 ) #+
  #geom_point(data= filter(qtls_3_4, chr3_acc %in% c('Hi', 'Wu', 'Sf')), shape=21, fill="red", color="darkred", size=3.5) 

ggsave(".tiff", plot = last_plot(), device = "tiff",
       dpi = 'retina')


acc <- c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi")

QTL3_effect_QTL4$acc <- acc
QTL3_effect_QTL4$chr4_QTL_genome <- factor(QTL3_effect_QTL4$chr4_QTL_genome, levels = c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct",
                                                                                        "Wil","Ws",'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi"))

QTL3_effect_QTL4$Rank <- 1:19
QTL3_effect_QTL4$`Edi effect` <- QTL3_effect_QTL4$`Edi effect` + QTL3_effect_QTL4$avg_score


formula <- y ~ x
ggplot(QTL3_effect_QTL4, aes(Rank, `Edi effect`)) +
  geom_text(aes(label=acc)) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(size=11, angle = 45, hjust=1),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=12, face='bold'),
        axis.line = element_line(colour = "black")) +
  geom_smooth(method= 'lm', formula = formula) +
  stat_regline_equation(aes(label =  paste(..eq.label..)),
                        label.x = 13, label.y = 1.0) +
  stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~','~")), digits = 3, label.y = 1.1, label.x = 13) +
  scale_x_continuous(breaks = round(seq(min(chr3_qtl_avg_root_score$Rank), max(chr3_qtl_avg_root_score$Rank), by = 1),1)) +
  scale_y_continuous(breaks = round(seq(min(chr3_qtl_avg_root_score$avg_root_score), max(chr3_qtl_avg_root_score$avg_root_score), by = 0.2),1)) +
  labs(x=("Parental Rank \nSensitive-Tolerant (1-19)"),
       y=('Avg Root Score'),  
       title = ('NIP1 QTL minus Mt & Bur QTL3'))  

ggsave("No_Mt_Bur_corr.tiff", plot = last_plot(), device = "tiff",
       dpi = 300)

qtls3_4_no_Mt <- qtls_3_4 %>% 
  filter(chr3_acc != "Mt" ) %>% 
  #filter(chr3_acc != "Mt" & chr3_acc != 'Bur') %>% 
  group_by(chr4_acc) %>% 
  summarise(avg_sc = round(mean(score),3)) %>% 
  mutate(Rank = c(1:19))

qtls3_4_no_Sf_Tsu <- qtls_3_4 %>% 
  #filter(chr3_acc != "Sf" ) %>% 
  filter(chr3_acc != "Sf" & chr3_acc != 'Tsu') %>% 
  group_by(chr4_acc) %>%
  summarise(avg_sc = round(mean(score),3)) %>% 
  mutate(Rank = c(1:19))

  
formula <- y ~ x
ggplot(qtls3_4_no_Sf_Tsu, aes(Rank, avg_sc)) +
  geom_text(aes(label=chr4_acc)) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(size=11, angle = 45, hjust=1),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=12, face='bold'),
        axis.line = element_line(colour = "black")) +
  geom_smooth(method= 'lm', formula = formula) +
  stat_regline_equation(aes(label =  paste(..eq.label..)),label.x = 13, label.y = 1.0) +
  stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~','~")), digits = 4, label.y = 1.1, label.x = 13) +
  scale_x_continuous(breaks = round(seq(min(chr3_qtl_avg_root_score$Rank), max(chr3_qtl_avg_root_score$Rank), by = 1),1)) +
  scale_y_continuous(breaks = round(seq(min(chr3_qtl_avg_root_score$avg_root_score), max(chr3_qtl_avg_root_score$avg_root_score), by = 0.2),1)) +
  labs(x=("Parental Rank \nSensitive-Tolerant (1-19)"),
       y=('Avg Root Score'),  
       title = ('NIP1 QTL minus Sf QTL3'))  


# made corr plots for root score vs rank without magic lines w/ QTL3 of each haplotype. Put data into excel and copy here
library(clipr)
corr_tbl_no_QTL3_haplo <- read_clip_tbl()
corr_tbl_no_QTL3_haplo <- corr_tbl_no_QTL3_haplo %>% 
  mutate(diff=(R2-0.4415))
# n magic lines taken out when each QTL3 haplotype taken out.  0 for original corr.  Then Can, Oy, Po, etc.
corr_tbl_no_QTL3_haplo$w.out_QTL3_acc <- factor(corr_tbl_no_QTL3_haplo$w.out_QTL3_acc, c("All","Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi"))
corr_tbl_no_QTL3_haplo <- corr_tbl_no_QTL3_haplo %>% 
  arrange(w.out_QTL3_acc)
corr_tbl_no_QTL3_haplo$n_QTL3_haplo <- c(0,15,22,9,21,28,8,22,26,21,15,15,22,18,21,29,23,24,32,13)
# arrange by diff
corr_tbl_no_QTL3_haplo <- corr_tbl_no_QTL3_haplo %>% 
  arrange(diff)
corr_tbl_no_QTL3_haplo$w.out_QTL3_acc <- factor(corr_tbl_no_QTL3_haplo$w.out_QTL3_acc, c("Mt","Wu","Kn","Bur","Wil","Ct","Po","Ws","Edi","Rsch",
                                                                                         "Col",'All',"Ler","Zu","Hi","No","Oy","Can","Tsu","Sf"))

# bar plot w/ All at 0, each parental haplotype w/ n=_ above/below bar
#corr_tbl_no_QTL3_haplo$w.out_QTL3_acc <- factor(corr_tbl_no_QTL3_haplo$w.out_QTL3_acc, levels = c('Kn','Wil','Hi','Col','Rsch','Can','Sf','Wu','Zu','Ler','All','Po','Ct','Tsu','Ws','No','Edi','Oy','Bur','Mt'))

corr_tbl_no_QTL3_haplo <- corr_tbl_no_QTL3_haplo %>% 
  mutate(net= ifelse(diff>0,'Positive', ifelse(diff<0,'Negative', 'None')))

write.csv(corr_tbl_no_QTL3_haplo,'corr_tbl_no_QTL3_haplo.csv')
# just use the above table

ggplot(corr_tbl_no_QTL3_haplo,aes(x=w.out_QTL3_acc, y=diff)) +
  geom_col(aes(fill = net)) +
  geom_text(data = . %>% filter(diff>0), aes(x=w.out_QTL3_acc, y=-0.008, label=paste('n =', n_QTL3_haplo), hjust=0.5, vjust=-1, angle=45 )) +
  geom_text(data = . %>% filter(diff<0), aes(x=w.out_QTL3_acc, y= 0.002, label=paste('n =', n_QTL3_haplo), hjust=0.01  , vjust= 0.5, angle=45 )) +
  theme(axis.text.x = element_text(size=11, angle = 45, hjust=1),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=12, face='bold'),
        axis.line = element_line(colour = "black"),
        legend.position = 'none') +
  labs(x= 'MAGIC Line Haplotype (@ QTL3) Removed',
       y= bquote(bold(Chang~to~R^~2))) +
  annotate("text", x = 14, y = -0.02, family = "Poppins", size = 4, color = "gray20",
           label = bquote(bold(R^{2}==0.4415))) +
  annotate("text", x = 14, y = -0.025, family = "Poppins", size = 3.5, color = "gray20",
           label = '(No Magic Lines Removed)') +
  annotate("text", x = 16, y = 0.04, family = "Poppins", size = 4, color = "gray20",
           label = bquote(bold(R^{2}==0.4922))) +
  annotate("text", x = 5.6, y = -0.055, family = "Poppins", size = 4, color = "gray20", #lineheight = .6, 
           label = bquote(bold(R^{2}==0.3745))) +
  geom_curve(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.07, "inch")), size = 0.8,
             color = "gray20", curvature = -0.4)

ggsave("R2_change_Haploytpes_removed.tiff", plot = last_plot(), device = "tiff",
       dpi = 800)

# had to tinker with these a lot for arrows
arrows <- tibble(x1 = c(12.4, 16.5, 5),
                 x2 = c(12.0, 19.4, 1.5),
                 y1 = c(-0.02, 0.042, -0.059), 
                 y2 = c(0, 0.05, -0.067))






################################################
# QTL3s effect on QTL4
avg_summary_QTL4 <- qtls_3_4 %>% 
  group_by(chr4_acc) %>% 
  summarise(avg_score=round(mean(score),3))
Parents<-c("Bur", "Can", "Col", "Ct", "Edi", "Hi", "Kn", "Ler", "Mt", "No", "Oy", "Po", "Rsch", "Sf", "Tsu", "Wil", "Ws", "Wu", "Zu")

# make table for avg root score when you take out each getotype from QTL on chr 3 
for (i in Parents){ 
  no <- filter(qtls_3_4, !chr3_acc %in% i)
  summary <- no %>% 
    group_by(chr4_acc) %>% 
    summarise(avg_sc = round(mean(score),3)) 
    colnames(summary)[2] <- paste('avg_no ',i)
  avg_summary_QTL4 <- left_join(avg_summary_QTL4, summary, by = 'chr4_acc')
}

QTL3_effect_QTL4 = data.frame(avg_summary_QTL4[,1:2], avg_summary_QTL4$avg_score - avg_summary_QTL4[-c(1:2)]) 

colnames(QTL3_effect_QTL4) <- c('chr4_QTL_genome', 'avg_score' ,"Can effect","Oy effect","Po effect","Rsch effect","Zu effect","Tsu effect",
                                "Wu effect","Ct effect","Wil effect","Ws effect",'Bur effect',"Col effect","Ler effect","Mt effect",
                                "No effect","Hi effect","Sf effect","Kn effect","Edi effect") 







QTL3_effect_QTL4_summary <- as_tibble(mapply(sum, QTL3_effect_QTL4[,-(1:2)]), rownames = NA) 

row.names(QTL3_effect_QTL4_summary) <-c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",
                                        'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi")

colnames(QTL3_effect_QTL4_summary) <- 'QTL3_effect_on_QTL4'

QTL3_effect_QTL4_summary <- rownames_to_column(QTL3_effect_QTL4_summary, 'Founder')

QTL3_effect_QTL4_summary$Founder <- factor(QTL3_effect_QTL4_summary$Founder, levels = c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",
                                                                                      'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi"))
z2 <- chr3_qtl %>%
  rename(Founder=acc)
z2 <- z2 %>%   
  group_by(Founder) %>% 
  tally()
QTL3_effect_QTL4_summary <- QTL3_effect_QTL4_summary %>%
  left_join(z2, by='Founder') 
QTL3_effect_QTL4_summary <- QTL3_effect_QTL4_summary %>% 
  mutate(avg_effect_QTL3_on_QTL4 = round(QTL3_effect_on_QTL4/n, 5))

ggplot(QTL3_effect_QTL4_summary, aes(x=Founder, y = avg_effect_QTL3_on_QTL4)) +
  geom_col()

ggsave("avg_QTL3_effect_on_QTL4.tiff", plot = last_plot(), device = "tiff",
       dpi = 'retina')



n_summary_QTL3 <- qtls_3_4 %>% 
  group_by(chr3_acc) %>% 
  tally()

for (i in Parents){ 
  no <- filter(qtls_3_4, !chr4_acc %in% i)
  summary <- no %>% 
    group_by(chr3_acc) %>% 
    tally()  
  colnames(summary)[2] <- paste(i, ' chr4') 
  n_summary_QTL3 <- left_join(n_summary_QTL3, summary, by = 'chr3_acc')
}
n_QTL3_no_QTL4 <-  data.frame(n_summary_QTL3[,1:2], n_summary_QTL3$n - n_summary_QTL3[-c(1:2)]) 



library(gridExtra)

grid.table(QTL3_effect_QTL4_summary, rows=NULL)
# save from R graphics device
#############################
# QTL4s effect on QTL3
avg_summary_QTL3 <- qtls_3_4 %>% 
  group_by(chr3_acc) %>% 
  summarise(avg_score=round(mean(score),3))

# make table for avg root score when you take out each getotype from QTL on chr 4 
# for (i in Parents){ 
#   no <- filter(qtls_3_4, !chr4_acc %in% i)
#   summary <- no %>% 
#     group_by(chr3_acc) %>% 
#     summarise(avg_sc = round(mean(score),3)) 
#   colnames(summary)[2] <- paste('avg_no ',i)
#   avg_summary_QTL3 <- left_join(avg_summary_QTL3, summary, by = 'chr3_acc')
# }
# 
# QTL4_effect_QTL3 = data.frame(avg_summary_QTL3[,1:2], avg_summary_QTL3$avg_score - avg_summary_QTL3[-c(1:2)]) 
# 
# colnames(QTL4_effect_QTL3) <- c('chr3_QTL_genome', 'avg_score' ,"Can effect","Oy effect","Po effect","Rsch effect","Zu effect","Tsu effect",
#                                 "Wu effect","Ct effect","Wil effect","Ws effect",'Bur effect',"Col effect","Ler effect","Mt effect",
#                                 "No effect","Hi effect","Sf effect","Kn effect","Edi effect") 
# 
# QTL4_effect_QTL3_summary <- as_tibble(mapply(sum, QTL4_effect_QTL3[,-(1:2)]), rownames = NA) 
# 
# row.names(QTL4_effect_QTL3_summary) <-c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",
#                                         'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi")
# 
# colnames(QTL4_effect_QTL3_summary) <- 'QTL4_effect_on_QTL3'
# 
# QTL4_effect_QTL3_summary <- rownames_to_column(QTL4_effect_QTL3_summary, 'Founder')
# 
# QTL4_effect_QTL3_summary$Founder <- factor(QTL4_effect_QTL3_summary$Founder, levels = c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",
#                                                                                       'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi"))
#  
# 
# z <- nip1_chr4_loc2 %>%
#   rename(Founder=acc)
# z <- z %>%   
#   group_by(Founder) %>% 
#   tally() 
# 
# QTL4_effect_QTL3_summary <- QTL4_effect_QTL3_summary %>%
#   left_join(z, by='Founder')
# QTL4_effect_QTL3_summary <- QTL4_effect_QTL3_summary %>% 
#   mutate(avg_effect_QTL4_on_QTL3 = round(QTL4_effect_on_QTL3/n, 5))
# 
# ggplot(QTL4_effect_QTL3_summary, aes(x=Founder, y = avg_effect_QTL4_on_QTL3)) +
#   geom_col()
# 
# ggsave("avg_QTL4_effect_on_QTL3.tiff", plot = last_plot(), device = "tiff",
#        dpi = 'retina')
# 
# grid.table(QTL4_effect_QTL3_summary, rows=NULL)
#################


##########################

# to view why bur is having such a big effect, look at what chr 3 genotypes are affected
bur_qtl4_effect <- QTL4_effect_QTL3 %>% 
  select(chr3_QTL_genome, avg_score, `Bur effect`) %>% 
  mutate(Founder=chr3_QTL_genome) %>% 
  left_join(z2, by='Founder') %>% 
  arrange(desc(`Bur effect`)) %>% 
  select(-Founder) %>% 
  rename(n_chr4 = n, `Bur_chr4 effect` = `Bur effect`, chr3_acc = chr3_QTL_genome) 

chr4bur_n <- n_QTL3_no_QTL4 %>% 
  select(chr3_acc, Bur..chr4)

bur_qtl4_effect <- bur_qtl4_effect %>% 
  left_join(chr4bur_n, by='chr3_acc')


grid.table(bur_qtl4_effect, rows=NULL)

ct_qtl4_effect <- QTL4_effect_QTL3 %>% 
  select(chr3_QTL_genome, avg_score, `Ct effect`) %>% 
  mutate(Founder=chr3_QTL_genome) %>% 
  left_join(z2, by='Founder') %>% 
  arrange(desc(`Ct effect`)) %>% 
  select(-Founder) %>% 
  rename(n_chr4 = n, `Ct_chr4 effect` = `Ct effect`)

grid.table(ct_qtl4_effect, rows=NULL)

# to view why Tsu & Mt are having such a big effect, look at what chr 4 genotypes are affected

tsu_qtl3_effect <- QTL3_effect_QTL4 %>% 
  select(chr4_QTL_genome, avg_score, `Tsu effect`) %>%
  mutate(Founder=chr4_QTL_genome) %>% 
  left_join(z, by='Founder') %>% 
  arrange(desc(`Tsu effect`)) %>% 
  select(-Founder) %>% 
  rename(n_chr3 = n, `Tsu_chr3 effect` = `Tsu effect`)

grid.table(bur_qtl4_effect, rows=NULL)

ct_qtl4_effect <- QTL4_effect_QTL3 %>% 
  select(chr3_QTL_genome, avg_score, `Ct effect`) %>% 
  mutate(Founder=chr3_QTL_genome) %>% 
  left_join(z2, by='Founder') %>% 
  arrange(desc(`Ct effect`)) %>% 
  select(-Founder) %>% 
  rename(n_chr4 = n, `Ct_chr4 effect` = `Ct effect`)

grid.table(ct_qtl4_effect, rows=NULL)


############## heatmap qtl effects
ggplot(QTL4_effect_QTL3, aes(x=x, y=chr3_QTL_genome, fill= colors)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") 

dat <- (QTL4_effect_QTL3[,3:21])
row.names(dat) <- QTL4_effect_QTL3$chr3_QTL_genome

library(pheatmap)
pheatmap(dat, display_numbers = T, cluster_cols=F, cluster_rows=F)

combos <- qtls_3_4 %>% 
  group_by(chr3_acc, chr4_acc) %>% 
  summarise(avg_sc = round(mean(score),4)) %>% 
  rename(QTL3_Founder=chr3_acc, QTL4_Founder = chr4_acc, avg_score=avg_sc)

ggplot(combos, aes(x=QTL3_Founder, y=QTL4_Founder)) + 
  geom_tile(aes(fill=avg_score)) +
  scale_fill_viridis_c(option = "plasma") +
  geom_text(aes(label = round(avg_score, 2)), color='white', size=2.5, fontface ='bold') 
  
effect_combos <- QTL4_effect_QTL3
library(viridis)
# med_summary <- qtls_3_4 %>% 
#   group_by(chr4_acc) %>% 
#   summarise(med_score=round(median(score),4))
# 
# for (i in Parents){ 
#   no <- filter(qtls_3_4, !chr3_acc %in% i)
#   summary <- no %>% 
#     group_by(chr4_acc) %>% 
#     summarise(med_sc = round(median(score),4)) 
#   colnames(summary)[2] <- paste('med_no ',i)
#   med_summary <- left_join(med_summary, summary, by = 'chr4_acc')
# }






# ###################################
# #create empty matrix
# qtl_mat <- matrix(nrow = 19, ncol = 19)
# colnames(qtl_mat) <- c('mean-Can_q3','mean-Oy_q3','mean-Po_q3','mean-Rsch_q3','mean-Zu_q3','mean-Tsu_q3','mean-Wu_q3','mean-Ct_q3',
#                       'mean-Wil_q3','mean-Ws_q3','mean-Bur_q3','mean-Col_q3','mean-Ler_q3','mean-Mt_q3','mean-No_q3','mean-Hi_q3','mean-Sf_q3',
#                       'mean-Kn_q3','mean-Edi_q3')
# rownames(qtl_mat) <- c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi")

 


###############################################################################
# horizontal bar plot showing magic_line number on y-axis, and range of DNA locus containing NIP1 gene from each magic_line spanning the x-axis
# with the actual NIP1;1 gene at x=0.  Every magic_line bar will cross x=0, but will vary in size and length in +/- direction.  log scale? 
library('scales')

first_50_magic <- nip1_chr4_loc2 %>% 
  mutate(upstream = (from.bp - 10421520), downstream = (to.bp - 10423508)) %>% 
  arrange(magic_line) %>% 
  slice(1:50) 

ggplot(first_50_magic) +
  geom_col(aes(x=as.factor(magic_line), y=upstream, fill='blue'), na.rm = TRUE, show.legend = F) +
  geom_col(aes(x=as.factor(magic_line), y=downstream, fill='red'), na.rm = TRUE, show.legend = F) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.title.y=element_text(size=12,face="bold"),
        axis.title.x=element_text(size=12,face="bold"),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size = 12, color = 'black', angle = 35, hjust=1),
        plot.title = element_text(size=16, margin = margin(t= -5))) +
  labs(y = expression(paste('Distance from ',italic('NIP1'),' locus (bp)')),
       x = expression(paste('MAGIC line ( # / genotype @ ',italic(' NIP1'),' locus )')),
       title = expression(paste("Ranges of Classified Genotypes @ ", italic('NIP1'), ' Locus, 1st 50 MAGIC Lines'))) +
  scale_y_continuous(breaks = c(-10000000,-7500000,-5000000,-2500000 ,0,2500000,5000000,7500000,10000000),
                     labels = comma) +
  geom_text(aes(x=as.factor(magic_line), y=upstream, label=magic_line, hjust=1.2 )) +
  geom_text(aes(x=as.factor(magic_line), y=downstream, label=acc, hjust=-0.2 )) 

ggsave("ranges_classified_genotypes_1st50.tiff", plot = last_plot(), device = "tiff", dpi = 'retina')

######################################

QTL <- nip1_chr4_loc2 %>% 
  mutate(upstream = (from.bp - 10421520), downstream = (to.bp - 10423508)) %>% 
  arrange(as.factor(acc)) %>% 
  filter(acc == 'Oy' )
# acc == 'Oy' | acc =='Po',
  
ggplot(QTL) +
  geom_col(aes(x=as.factor(magic_line), y=upstream, fill='blue', group = acc), show.legend = F) +
  geom_col(aes(x=as.factor(magic_line), y=downstream, fill='red', group = acc), show.legend = F) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.title.y=element_text(size=12,face="bold"),
        axis.title.x=element_text(size=12,face="bold"),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size = 12, color = 'black', angle = 35, hjust=1),
        plot.title = element_text(size=16, margin = margin(t= -5))) +
  labs(y = expression(paste('Distance from ',italic('NIP1'),' locus (bp)')),
       x = expression(paste('MAGIC line ( # / genotype @ ',italic(' NIP1'),' locus )')),
       title = expression(paste("Ranges of Classified Genotypes @ ", italic('NIP1'), ' Locus, OY'))) +
  scale_y_continuous(breaks = c(-10000000,-7500000,-5000000,-2500000,-1000000,-500000 ,0,1000000, 2500000,5000000,7500000,10000000),
                     labels = comma) +
  geom_text(aes(x=as.factor(magic_line), y=upstream, label=magic_line, hjust=1.2 )) +
  geom_text(aes(x=as.factor(magic_line), y=downstream, label=acc, hjust=-0.2 )) 

ggsave("ranges_classified_genotypes_rootscore_under_oy.tiff", plot = last_plot(), device = "tiff",
       dpi = 'retina')


