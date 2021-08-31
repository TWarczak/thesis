library(agricolae)
# data("PlantGrowth")
# plant.lm <- lm(weight ~ group, data = PlantGrowth)
# plant.av <- aov(plant.lm)

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
library(gridExtra)

tbl <- tableGrob(tukey_t, rows=NULL)

grid.table(tukey_t, rows=NULL)
# save from R graphics device

###################
# t-tests 0dAs to 1dAs

t_tests <- All_Reps_fc %>% 
  ungroup() %>% 
  select(Parental_ID, Condition, NIP1_Exp)  %>% 
  filter(Condition=='0dAs' | Condition=='1dAs') %>% 
  group_by(Condition) %>% 
  data.frame()

Bur_Zu <- data.frame(Parental_ID=c('Bur', 'Zu'), Condition=c('1dAs', '0dAs'), NIP1_Exp=c(mean(0.06844044, 0.09489514), mean(0.09735566, 0.11900133)))
t_tests <- rbind(t_tests, Bur_Zu) %>% 
  group_by(Parental_ID, Condition)

View(count(t_tests))


Bur <- t_tests %>% 
  filter(Parental_ID=='Bur')
  
Bur_t <- t.test(NIP1_Exp~Condition, data=Bur,
                var.equal=T,
                conf.level=0.95)

Can <- t_tests %>% 
  filter(Parental_ID=='Can')

Can_t <- t.test(NIP1_Exp~Condition, data=Can,
                var.equal=T,
                conf.level=0.95)

Col <- t_tests %>% 
  filter(Parental_ID=='Col')

Col_t <- t.test(NIP1_Exp~Condition, data=Col,
                var.equal=T,
                conf.level=0.95)

Ct <- t_tests %>% 
  filter(Parental_ID=='Ct')

Ct_t <- t.test(NIP1_Exp~Condition, data=Ct,
               var.equal=T,
               conf.level=0.95)

Edi <- t_tests %>% 
  filter(Parental_ID=='Edi')

Edi_t <- t.test(NIP1_Exp~Condition, data=Edi,
                var.equal=T,
                conf.level=0.95)

Hi <- t_tests %>% 
  filter(Parental_ID=='Hi')

Hi_t <- t.test(NIP1_Exp~Condition, data=Hi,
               var.equal=T,
               conf.level=0.95)

Kn <- t_tests %>% 
  filter(Parental_ID=='Kn')

Kn_t <- t.test(NIP1_Exp~Condition, data=Kn,
               var.equal=T,
               conf.level=0.95)

Ler <- t_tests %>% 
  filter(Parental_ID=='Ler')

Ler_t <- t.test(NIP1_Exp~Condition, data=Ler,
               var.equal=T,
               conf.level=0.95)

Mt <- t_tests %>% 
  filter(Parental_ID=='Mt')

Mt_t <- t.test(NIP1_Exp~Condition, data=Mt,
               var.equal=T,
               conf.level=0.95)

No <- t_tests %>% 
  filter(Parental_ID=='No')

No_t <- t.test(NIP1_Exp~Condition, data=No,
               var.equal=T,
               conf.level=0.95)

Oy <- t_tests %>% 
  filter(Parental_ID=='Oy')

Oy_t <- t.test(NIP1_Exp~Condition, data=Oy,
               var.equal=T,
               conf.level=0.95)

Po <- t_tests %>% 
  filter(Parental_ID=='Po')

Po_t <- t.test(NIP1_Exp~Condition, data=Po,
               var.equal=T,
               conf.level=0.95)

Rsch <- t_tests %>% 
  filter(Parental_ID=='Rsch')

Rsch_t <- t.test(NIP1_Exp~Condition, data=Rsch,
               var.equal=T,
               conf.level=0.95)

Sf <- t_tests %>% 
  filter(Parental_ID=='Sf')

Sf_t <- t.test(NIP1_Exp~Condition, data=Sf,
               var.equal=T,
               conf.level=0.95)

Tsu <- t_tests %>% 
  filter(Parental_ID=='Tsu')

Tsu_t <- t.test(NIP1_Exp~Condition, data=Tsu,
               var.equal=T,
               conf.level=0.95)

Wil <- t_tests %>% 
  filter(Parental_ID=='Wil')

Wil_t <- t.test(NIP1_Exp~Condition, data=Wil,
               var.equal=T,
               conf.level=0.95)

Ws <- t_tests %>% 
  filter(Parental_ID=='Ws')

Ws_t <- t.test(NIP1_Exp~Condition, data=Ws,
               var.equal=T,
               conf.level=0.95)

Wu <- t_tests %>% 
  filter(Parental_ID=='Wu')

Wu_t <- t.test(NIP1_Exp~Condition, data=Wu,
               var.equal=T,
               conf.level=0.95)

Zu <- t_tests %>% 
  filter(Parental_ID=='Zu')

Zu_t <- t.test(NIP1_Exp~Condition, data=Zu,
               var.equal=T,
               conf.level=0.95)

Bur_t # not sig
Can_t # not sig
Col_t # sig
Ct_t  # not sig
Edi_t # not sig
Hi_t  # not sig
Kn_t  # sig
Ler_t # not sig
Mt_t  # not sig
No_t  # not sig
Oy_t  # not sig
Po_t  # not sig
Rsch_t # not sig
Sf_t  # not sig
Tsu_t # not sig
Wil_t # not sig
Ws_t  # not sig
Wu_t  # not sig
Zu_t  # not sig

