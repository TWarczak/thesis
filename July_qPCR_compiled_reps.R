library(tidyverse)

Rep3_0d_1d <- read_csv("~/R/qPCR/July2019/ToddW_MAGIC_Parentals_UBQ1_NIP1_July9_2019_0d10As_1d10As.csv")
Rep3_2d_3d <- read_csv("~/R/qPCR/July2019/ToddW_MAGIC_Parentals_UBQ1_NIP1_July9_2019_2d10As_3d10As.csv")
Rep3_4d_5d <- read_csv("~/R/qPCR/July2019/ToddW_MAGIC_Parentals_UBQ1_NIP1_July17_2019_4d10As_5d10As.csv")
Rep2_2d_3d <- read_csv("~/R/qPCR/July2019/ToddW_MAGIC_Parentals_UBQ1_NIP1_July6_2019_2d10As_3d10As.csv")
Rep2_4d_5d <- read_csv("~/R/qPCR/July2019/ToddW_MAGIC_Parentals_UBQ1_NIP1_July7_2019_4d10As_5d10As.csv")
Rep2_Rep3_leftovers <- read_csv(('~/R/qPCR/July2019/ToddW_MAGIC_Parentals_UBQ1_NIP1_July18_2019_leftovers.csv'))
Rep2_0d_1d <- read_csv("~/R/qPCR/April2019/ToddW_MAGIC_Parentals_UBQ1_NIP1_April24_2019_7dCon_1d10As.csv")
Rep2_0d_1d <- select(Rep2_0d_1d, c('Target','Sample','Cq'))
Rep_reruns <- read_csv('~/R/qPCR/July2019/ToddW_MAGIC_Parentals_UBQ1_NIP1_July24_2019_reruns.csv') 
Rep_reruns2 <- read_csv('~/R/qPCR/July2019/ToddW_MAGIC_Parentals_UBQ1_NIP1_Sept_9_2019_reruns.csv')
########################################################
# Rep 3

# bind_rows for all Rep 3 samples
Rep3_leftovers <- Rep2_Rep3_leftovers %>% 
                    filter(str_detect(Sample, "Bur3") | str_detect(Sample, 'Kn')) 

Rep3_leftovers$Sample <- Rep3_leftovers$Sample %>%
  str_replace_all("Bur3", "Bur")

Rep3 <- bind_rows(Rep3_0d_1d, Rep3_2d_3d, Rep3_4d_5d, Rep3_leftovers) %>%
# Separating the Sample variable into Parental_ID and Condition variables
  separate(col = Sample, into = c("Parental_ID", "Condition"), sep = " ")

Rep3$Target <- as.factor(Rep3$Target)
Rep3$Condition <- as.factor(Rep3$Condition)
Rep3$Parental_ID <- as.factor(Rep3$Parental_ID)  # for plot at end

Rep3$Condition <- factor(Rep3$Condition, levels = c("0dAs", "1dAs", "2dAs", "3dAs", "4dAs", "5dAs"))

# Need to remove Mt 2dAs UBQ rows 398 & 399 and Zu 0d 
# Rows can be dropped with negative indices:
Rep3 <- slice(Rep3, -c(398:399, 118:120, 103:105))

###########################################################
# Rep2

# bing_rows for all Rep 2 samples
Rep2_leftovers <- Rep2_Rep3_leftovers %>% 
  filter(!str_detect(Sample, "Bur3") & !str_detect(Sample, 'Kn') & !str_detect(Sample, 'nip1_1')) 

Rep2_leftovers$Sample <- Rep2_leftovers$Sample %>%
  str_replace_all(c("Bur2"= "Bur", 'Ct2'= 'Ct', 'CtE'= 'Ct'))

Rep2 <- bind_rows(Rep2_0d_1d, Rep2_2d_3d, Rep2_4d_5d, Rep2_leftovers) %>%
  # Separating the Sample variable into Parental_ID and Condition variables
  separate(col = Sample, into = c("Parental_ID", "Condition"), sep = " ")

Rep2$Target <- as.factor(Rep2$Target)
Rep2$Condition <- as.factor(Rep2$Condition)
Rep2$Parental_ID <- as.factor(Rep2$Parental_ID)  # for plot at end

Rep2$Condition <- factor(Rep2$Condition, levels = c("0dAs", "1dAs", "2dAs", "3dAs", "4dAs", "5dAs"))

# Need to remove Mt 1dAs NIP1 row 87 & Rsch 2dAs NIP1 row 227   clear pipetting errors
# Rows can be dropped with negative indices:
Rep2 <- slice(Rep2, -c(87,227))

#################################################################
# Rep1 

# Measure expression of NIP1 after 7d, then expose plants to arsenic, measure NIP1 expression every day for 5 days
MAGIC_10As_Assay <- read_csv("~/R/qPCR/March2019/MAGIC_10As_Assay.csv")

# Only want 3 variables and don't need Col 15As data
Rep1 <- MAGIC_10As_Assay %>%
  select(Target, Sample, Cq) %>%  
  filter(!grepl('15', Sample))    

# Col data from Feb was in different Excel format; matching to the new qPCR data
Rep1$Sample[649:654]="Col-Con"
Rep1$Sample[655:660]="Col-1dAs"
Rep1$Sample[661:666]="Col-2dAs"
Rep1$Sample[667:672]="Col-3dAs"
Rep1$Sample[673:678]="Col-4dAs"
Rep1$Sample[679:684]="Col-5dAs"

# Separating the Sample variable into Parental_ID and Condition variables
Rep1 <- separate(Rep1, col = Sample, into = c("Parental_ID", "Condition"), sep = "-")

Rep1$Target <- as.factor(Rep1$Target)
Rep1$Condition <- as.factor(Rep1$Condition)
Rep1$Parental_ID <- as.factor(Rep1$Parental_ID)  # for plot at end

levels(Rep1$Condition)[levels(Rep1$Condition)=="Con"] <- "0dAs"
levels(Rep1$Target)[levels(Rep1$Target)=="UBQ1"] <- "UBQ"   # Used UBQ1 in the Col data

# Assign levels to Condition factor
Rep1$Condition <- factor(Rep1$Condition, levels = c("0dAs", "1dAs", "2dAs", "3dAs", "4dAs", "5dAs"))
# Data Cleaning: removing obvious outlier reps, averaging two remaining reps
Rep1 <- slice(Rep1, -c(481, 550, 526, 589, 619, 639, 520, 346, 325, 353, 389, 538, 561, 294))


###############################################################
# Ruruns

Reruns <- Rep_reruns %>%
  # Separating the Sample variable into Parental_ID and Condition variables
  separate(col = Sample, into = c("Parental_ID", "Condition"), sep = " ")

Reruns$Target <- as.factor(Reruns$Target)
Reruns$Condition <- as.factor(Reruns$Condition)
Reruns$Parental_ID <- as.factor(Reruns$Parental_ID)  # for plot at end

Reruns$Condition <- factor(Reruns$Condition, levels = c("0dAs", "1dAs", "2dAs", "3dAs", "4dAs", "5dAs"))

# Added Rep column in Reruns, Rep 1-3, so add Rep 4,5,6 to others. bad planning
Rep1 <- mutate(Rep1, Rep = 4)
Rep2 <- mutate(Rep2, Rep = 5)
Rep3 <- mutate(Rep3, Rep = 6)

###############################################################
# Ruruns2
Reruns2 <- Rep_reruns2 %>%
  # Separating the Sample variable into Parental_ID and Condition variables
  separate(col = Sample, into = c("Rep","Parental_ID", "Condition"), sep = " ")

Reruns2$Target <- as.factor(Reruns2$Target)
Reruns2$Parental_ID <- as.factor(Reruns2$Parental_ID)  # for plot at end
Reruns2$Condition <- factor(Reruns2$Condition, levels = c("0dAs", "1dAs", "2dAs", "3dAs", "4dAs", "5dAs"))
Reruns2$Rep <- as.numeric(Reruns2$Rep)
# Added Rep column in Reruns, Rep 1-3, so add Rep 4,5,6 to others. bad planning

###############################################################
# All data cleaned and combined for final plot
Reruns2
All_Reps <- bind_rows(Rep1, Rep2, Rep3, Reruns, Reruns2)
Wil <- filter(All_Reps, Parental_ID == 'Wil', Condition == '2dAs')
# Averaging the Cq values for each technical rep
All_Reps_avg <- arrange(All_Reps, Rep, Parental_ID, Condition, Target, Cq) %>%
  group_by(Parental_ID, Rep, Condition, Target) %>%
  summarise(avg_Cq = mean(Cq))

# Want to compare UBQ and NIP values
All_Reps_spread <- spread(All_Reps_avg, Target, avg_Cq)

# Calculate fold change of NIP1 expression
All_Reps_fc <- mutate(All_Reps_spread,NIP1_Exp = 2^(UBQ - NIP1)) %>%
  arrange(Parental_ID, Condition, Rep) 

# Find the value of average Col 0dAs NIP1 expression
Col_avg_exp_0dAs <- mean(All_Reps_fc$NIP1_Exp[42:44])

# Normalize each NIP1 exp to Col 0d, then group the reps, average the normalized expression of the reps
All_Reps_norm <- All_Reps_fc %>%
  mutate(NIP1_Exp_Norm = NIP1_Exp / Col_avg_exp_0dAs) %>%
  group_by(Parental_ID, Condition) %>%
  summarise(avg_NIP1_Exp = mean(NIP1_Exp_Norm), sd = sd(NIP1_Exp_Norm), se = sd/sqrt(3))

# For samples found in se_reruns below, se needs to be sd/sqrt(6))
se_reruns <- bind_rows(Reruns, Reruns2)
se_reruns <- arrange(se_reruns, Rep, Parental_ID, Condition, Target, Cq) %>%
  group_by(Parental_ID, Rep, Condition, Target) %>%
  summarise(avg_Cq = mean(Cq))
se_reruns <- spread(se_reruns, Target, avg_Cq)
se_reruns <- mutate(se_reruns,NIP1_Exp = 2^(UBQ - NIP1)) %>%
  arrange(Parental_ID, Condition, Rep) 
se_reruns <- se_reruns %>%
  mutate(NIP1_Exp_Norm = NIP1_Exp / Col_avg_exp_0dAs) %>%
  group_by(Parental_ID, Condition) %>%
  summarise(avg_NIP1_Exp = mean(NIP1_Exp_Norm), sd = sd(NIP1_Exp_Norm), se = sd/sqrt(3)) %>% 
  filter(is.na(se)== F) %>% 
  select(Parental_ID, Condition)
# I'm positive there's an easier way to do this
R <- semi_join(All_Reps_norm, se_reruns, by = c('Parental_ID', 'Condition'))
R$se <- R$sd/sqrt(6)
# data.table method for replacing values in 1 dt with values in another dt
library(data.table)
setDT(All_Reps_norm)
setDT(R)

All_Reps_norm[R, on = c("Parental_ID", 'Condition'), se := i.se]

# Order levels of Parental_ID by sensitivity to arsenic 
All_Reps_norm$Parental_ID <- factor(All_Reps_norm$Parental_ID, levels = c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi"))

#write_csv(All_Reps_norm, 'All_Reps_norm.csv', col_names = T)
#write_csv(All_Reps_fc, 'All_Reps_fc.csv', col_names = T)
#write_csv(All_Reps_avg, 'All_Reps_avg.csv', col_names = T)


ggplot(All_Reps_norm, aes(x = Condition, y = avg_NIP1_Exp, color = Parental_ID)) +
  geom_line(aes(group = Parental_ID), size=1) +
  geom_point() +
  facet_wrap(. ~ Parental_ID) +
  ylim(0.3, 1.5) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=12,face="bold"),
        axis.title.x=element_text(size=12,face="bold"),
        axis.text.y=element_text(size = 11),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1),
        plot.title = element_text(size=16, margin = margin(t= -5))) +
  geom_hline(yintercept=1, 
             linetype="dashed", 
             color = "red") +
  labs(x = "\nDays Exposed to arsenic", 
       y = expression(paste(italic('NIP1 '), 'Exp.  ', italic('(rel. to Col 0dAs)'),'\n')),
       title = expression(paste("MAGIC Parental ", italic('NIP1 '), 'Expression'))) +
  geom_errorbar(aes(ymin=avg_NIP1_Exp-se, ymax=avg_NIP1_Exp+se), width=.13)  # Width of the error bars 

# ggsave("MAGIC_Parent_10As_qPCR_TimePoint_Sept30.tiff", plot = last_plot(), device = "tiff",
#        dpi = 'retina')

###############################################################
# Coorelation between Parental NIP1 Expr (0-5dAs) and phenotype rank (sensitive to tolerant 1-19)

nip1_exp_byparent <- All_Reps_norm %>%
  group_by(Parental_ID) %>%
  summarise(total_nip_exp = sum(avg_NIP1_Exp))  

nip1_exp_byparent$Parental_ID <- factor(nip1_exp_byparent$Parental_ID, levels = c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi"))

nip1_exp_byparent$Rank <- 1:19

formula <- y ~ x
ggplot(nip1_exp_byparent, aes(Rank, total_nip_exp, label = Parental_ID)) +
  geom_text() +
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
                        label.x = 13, label.y = 7.1) +
  stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~','~")), label.y = 7.5, label.x = 13) +
  scale_x_continuous(breaks = round(seq(min(nip1_exp_byparent$Rank), max(nip1_exp_byparent$Rank), by = 1),1)) +
  scale_y_continuous(breaks = round(seq(min(nip1_exp_byparent$total_nip_exp), max(nip1_exp_byparent$total_nip_exp), by = 0.5),1)) +
  labs(x=("\nParental Rank \nSensitive-Tolerant (1-19)"),
      y=(expression(paste("Cumulative  ", italic('NIP1 '), 'Exp. (rel. to Col 0dAs)'))))  

ggsave("MAGIC_Parent_qPCR_corr_Aug3.tiff", plot = last_plot(), device = "tiff",
       dpi = 600)


library(agricolae)

subset <- All_Reps_fc %>%
  mutate(NIP1_Exp_Norm = NIP1_Exp / Col_avg_exp_0dAs) %>% 
  ungroup()

subset2 <- subset %>% 
  mutate(both = str_c(subset$Parental_ID, subset$Condition)) %>% 
  filter(Condition=='0dAs')

subset2$both <- as.factor(subset2$both)

nip.lm <- lm(NIP1_Exp_Norm~both, data = subset2)
nip.aov <- aov(nip.lm)

summary(nip.aov)
tukey.test <- HSD.test(nip.aov, trt = 'both')
tukey.test

tukey_t <- read.table("clipboard", sep="", header = T)
tukey_t <- rownames_to_column(tukey_t)
tukey_t <- tukey_t %>%
  separate(col = rowname, into = c("Founder", "Condition"), sep = -4) %>% 
  rename(NIP1_Exp = NIP1_Exp_Norm) 

tukey_t[,3] <-round(tukey_t[,3],3) #the "-1" excludes column 1
tukey_t

All_Reps_norm_0dAs <- All_Reps_norm %>% 
  filter(Condition == '0dAs') %>% 
  left_join(tukey_t, by = c('Parental_ID'='Founder', 'Condition'))

All_Reps_norm_0dAs$Parental_ID <- factor(All_Reps_norm_0dAs$Parental_ID, levels = c("Can","Oy","Po","Rsch","Zu","Tsu","Wu","Ct","Wil","Ws",'Bur',"Col","Ler","Mt","No","Hi","Sf","Kn","Edi"))
All_Reps_norm_0dAs <- All_Reps_norm_0dAs %>% 
  arrange(Parental_ID)
All_Reps_norm_0dAs$Rank <- 1:19

############################################################
# Bar plot of OdAs NIP1;1 expression, with tukey test results

ggplot(All_Reps_norm_0dAs, aes(x= Parental_ID, y= avg_NIP1_Exp, fill=viridis_pal(option = "D")(19))) + 
  geom_col(show.legend = F) +
  coord_cartesian(ylim=c(0.3,1.5)) +
  geom_errorbar(aes(ymin=avg_NIP1_Exp-se, ymax=avg_NIP1_Exp+se), width=.15) +
  theme(axis.title.y=element_text(size=12,face="bold"),
        axis.title.x=element_text(size=12,face="bold"),
        axis.text.y=element_text(size = 11),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1)) +
  geom_hline(size = 1,yintercept=1, linetype="dashed", color = "red") +
  labs(y= expression(paste(italic('NIP1;1')~Expression~' 0d'~As^{III}~'(7d), '~"rel. to Col")),
       x= "MAGIC Founder Line") +
  geom_text(aes(x= Parental_ID, y= (avg_NIP1_Exp+se+.1), label=groups), hjust=0.5, vjust=1.8, size=5)

ggsave("MAGIC_Parent_qPCR_0dAs_Aug3.tiff", plot = last_plot(), device = "tiff",
       dpi = 600)

######################################
# Correlation plot
library(ggpubr)

formula <- y ~ x
ggplot(All_Reps_norm_0dAs, aes(Rank, NIP1_Exp)) +
  geom_text(aes(label=Parental_ID)) +
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
                        label.x = 14, label.y = 1.18) +
  stat_cor(aes(label=paste(..rr.label.., ..p.label.., sep = "~','~")), label.y = 1.25, label.x = 14) +
  scale_x_continuous(breaks = round(seq(min(nip1_exp_byparent$Rank), max(nip1_exp_byparent$Rank), by = 1),1)) +
  #scale_y_continuous(breaks = round(seq(min(nip1_exp_byparent$total_nip_exp), max(nip1_exp_byparent$total_nip_exp), by = 0.5),1)) +
  labs(x=("Parental Rank \nSensitive-Tolerant (1-19)"),
       y=(expression(paste(italic('NIP1;1')~Exp~'0d'~As^{III}~"(rel. to Col)"))))

ggsave("MAGIC_Parent_qPCR_0dAs_corr_Aug3.tiff", plot = last_plot(), device = "tiff",
       dpi = 600)

####################################
# overexpression lines 
OE_Can_NIP1 <- read.table("clipboard", sep="", header = T)
write_csv(OE_Can_NIP1, 'OE_Can_NIP1.csv', col_names = T)

OE_Reps_avg <- arrange(OE_Can_NIP1, Rep, Target, Sample, Cq) %>%
  group_by(Sample, Target, Rep) %>%
  summarise(avg_Cq = mean(Cq))

# Want to compare UBQ and NIP values
OE_Reps_spread <- spread(OE_Reps_avg, Target, avg_Cq)

# Calculate fold change of NIP1 expression
OE_Reps_fc <- mutate(OE_Reps_spread,NIP1_Exp = 2^(UBQ - NIP1.1)) %>%
  arrange(Sample, Rep) 

# Find the value of average Col 0dAs NIP1 expression
Col_avg_exp <- mean(OE_Reps_fc$NIP1_Exp[7:8])

# Normalize each NIP1 exp to Col 0d, then group the reps, average the normalized expression of the reps
OE_Reps_norm <- OE_Reps_fc %>%
  mutate(NIP1_Exp_Norm = NIP1_Exp / Col_avg_exp) %>%
  group_by(Sample) %>%
  summarise(avg_NIP1_Exp = mean(NIP1_Exp_Norm), sd = sd(NIP1_Exp_Norm), se = sd/sqrt(2)) # need one one rep for each
OE_Reps_norm$Sample <- factor(OE_Reps_norm$Sample, levels = c('Col', 'nip1;1', 'Can_NIP1_C1', 'Can_NIP1_C2', 'Can_NIP1_C3'))
library(viridis)
library(scales)
show_col(viridis_pal()(10))

ggplot(OE_Reps_norm, aes(x= Sample, y= avg_NIP1_Exp)) + 
  geom_col(aes(fill= Sample),show.legend = F) +
  scale_fill_manual(values=c("#482677FF", "#238A8DFF", "#73D055FF", "#73D055FF", "#73D055FF")) +
  scale_x_discrete(labels=c("Can_NIP1_C1" = "nip1;1/\nCan_NIP1_C1", "Can_NIP1_C2" = "nip1;1/\nCan_NIP1_C2", "Can_NIP1_C3" = "nip1;1/\nCan_NIP1_C3")) +
  geom_errorbar(aes(ymin=avg_NIP1_Exp-se, ymax=avg_NIP1_Exp+se), width=.15) +
  theme(axis.title.y=element_text(size=12,face="bold"),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size = 11),
        axis.text.x=element_text(size = 11, angle = 35, hjust = 1)) +
  geom_hline(size = 1,yintercept=1, linetype="dashed", color = "red") +
  labs(y= expression(paste(italic('NIP1;1')~Expression~" (rel. to Col)"))) 

ggsave("Can_NIP1_nip1_comp.tiff", plot = last_plot(), device = "tiff",
       dpi = 800)
