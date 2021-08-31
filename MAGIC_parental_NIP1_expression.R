MAGIC_Parentals<-c("Bur", "Can", "Col", "Ct", "Edi", "Hi", "Kn", "Ler", "Mt", "No", "Oy", "Po", "Rsch", "Sf", "Tsu", "Wil", "Ws", "Wu", "Zu")
phenotype<-c(8, 1, 12, 6, 19, 13, 18, 17, 16, 15, 2, 5, 4, 14, 9, 11, 7, 10, 3) #rank from most sensitive (1) to most tolerant (19)

NIP1_Expression_As10_14d<-c(2, 19, 4, 16, 1, 6, 5, 9, 17, 3, 18, 13, 8, 7, 12, 10, 11, 15, 14)
NIP1_Expression_Control_10d<-c(8, 6, 5, 10, 2, 9, 3, 19, 11, 4, 13, 15, 17, 7, 16, 18, 14, 12, 1)
NIP1_Expression_Control_7d <- c(16, 19, 11, 17, 1, 9, 4, 12, 7, 3, 14, 8, 13, 2, 18, 10, 5, 15, 6)
NIP1_Expression_As50_1d<- c(13, 8, 4, 3, 14, 9, 19, 10, 11, 6, 12, 18, 5, 17, 7, 16, 15, 2, 1)


MAGIC_matrix <- matrix(c(phenotype, NIP1_Expression_As10_14d), nrow=19, byrow = FALSE, dimnames = list(MAGIC_Parentals) )
MAGIC_matrix
colnames(MAGIC_matrix)<- c("Phenotype_Rank", "NIP1_Exp_As10")

rownames(MAGIC_matrix)<-MAGIC_Parentals
MAGIC_table<-data.frame(MAGIC_matrix)
MAGIC_table
 
MAGIC_table2 <- cbind(MAGIC_table, NIP1_Expression_Control_10d, NIP1_Expression_Control_7d, NIP1_Expression_As50_1d)
MAGIC_table2

colnames(MAGIC_table2)<- c("Phenotype_Rank", "NIP1_Exp_As10", "NIP1_Exp_Con", "NIP1_Exp_Con_7d", "NIP1_Exp_As50_1d")

str(MAGIC_table2)

SpearmanRank1 <-cor.test(MAGIC_table2$Phenotype_Rank, MAGIC_table2$NIP1_Exp_As10,  method = "spearman")
SpearmanRank2 <-cor.test(MAGIC_table2$Phenotype_Rank, MAGIC_table2$NIP1_Exp_Con,  method = "spearman")
SpearmanRank3 <-cor.test(MAGIC_table2$Phenotype_Rank, MAGIC_table2$NIP1_Exp_Con_7d,  method = "spearman")
SpearmanRank4 <-cor.test(MAGIC_table2$Phenotype_Rank, MAGIC_table2$NIP1_Exp_As50_1d,  method = "spearman")
SpearmanRank2



KendallRank <-cor.test(MAGIC_table2$Phenotype_Rank, MAGIC_table2$NIP1_Exp_As10,  method = "kendall")
KendallRank

#install.packages("ggpubr")
library("ggpubr")

ggscatter(MAGIC_table2, x = "Phenotype_Rank", y = "NIP1_Exp_As10", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coef.coord = c(10, 19), cor.method = "spearman",
          xlab = "Phenoptype Rank \n(Sensitive-Tolerant)", ylab = "NIP1_Expression", label = MAGIC_Parentals)

ggscatter(MAGIC_table2, x = "Phenotype_Rank", y = "NIP1_Exp_Con", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coef.coord = c(10, 19), cor.method = "spearman",
          xlab = "Phenoptype Rank \n(Sensitive-Tolerant)", ylab = "NIP1_Expression", label = MAGIC_Parentals)

ggscatter(MAGIC_table2, x = "Phenotype_Rank", y = "NIP1_Exp_Con_7d", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coef.coord = c(10, 19), cor.method = "spearman",
          xlab = "Phenoptype Rank \n(Sensitive-Tolerant)", ylab = "NIP1_Expression", label = MAGIC_Parentals)

ggscatter(MAGIC_table2, x = "Phenotype_Rank", y = "NIP1_Exp_As50_1d", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coef.coord = c(10, 19), cor.method = "spearman",
          xlab = "Phenoptype Rank \n(Sensitive-Tolerant)", ylab = "NIP1_Expression", label = MAGIC_Parentals)


ggscatter(MAGIC_table2, x = "Phenotype_Rank", y = "NIP1_Exp_As10", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.coef.coord = c(10, 19), cor.method = "kendall",
          xlab = "Phenoptype Rank \n(Sensitive-Tolerant)", ylab = "NIP1_Expression", label = MAGIC_Parentals)
?ggscatter

