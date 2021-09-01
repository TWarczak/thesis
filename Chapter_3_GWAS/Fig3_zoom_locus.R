# Chromosome 4 QTL, top markers near NIP1 plot
library(tidyverse)
#NIP1 segments

segment start end
5utr 10423508 10423409
exon1 10423408 10423241
intron1 10423240 10423129
exon2 10423128 10422904
intron2 10422903 10422480
exon3 10422479 10422261
intron3 10422260 10422110
exon4 10422109 10422048
intron4 10422047 10421944
exon5 10421943 10421727
3utr 10421726 10421520

NIP1_bp<- read.table(file = "clipboard", sep = " ", header=T)
NIP1_bp

nip <- data.frame(x1=NIP1_bp$start, x2=NIP1_bp$end, 
             y1=c( 0.15, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.15), 
             y2=c(-0.15,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.15), 
             segment=c('UTR','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon', 'UTR'))
lcr7
segment start end
5utr 10429000 10428976
exon1 10428975 10428915
intron1 10428914 10428001
exon2 10428000 10427807
3utr 10427806 10427615

lcr7_bp<- read.table(file = "clipboard", sep = " ", header=T)
lcr7_bp

lcr7 <- data.frame(x1=lcr7_bp$start, x2=lcr7_bp$end, 
                   y1=c( 0.15, 0.2, 0.06, 0.2, 0.15), 
                   y2=c(-0.15,-0.2,-0.06,-0.2,-0.15), 
                   segment=c('UTR','Exon','Intron','Exon','UTR'))

lcr15
segment start end
5utr 10431043 10430960
exon1 10430959 10430899
intron1 10430898 10430268
exon2 10430267 10430074
3utr 10430073 10429871

lcr15_bp<- read.table(file = "clipboard", sep = " ", header=T)
lcr15_bp

lcr15 <- data.frame(x1=lcr15_bp$start, x2=lcr15_bp$end, 
                    y1=c( 0.15, 0.2, 0.06, 0.2, 0.15), 
                    y2=c(-0.15,-0.2,-0.06,-0.2,-0.15), 
                    segment=c('UTR','Exon','Intron','Exon','UTR'))

edr2
segment start end
5utr 10437566 10437174
exon1 10437173 10437038
intron1 10437037 10436910
exon2 10436909 10436834
intron2 10436833 10436660
exon3 10436659 10436608
intron3 10436607 10436504
exon4 10436503 10436440
intron4 10436439 10436314
exon5 10436313 10436212
intron5 10436211 10435718
exon6 10435717 10435654
intron6 10435653 10435572
exon7 10435571 10435456
intron7 10435455 10435090
exon8 10435089 10435050
intron8 10435049 10434948
exon9 10434947 10434846
intron9 10434845 10434556
exon10 10434555 10434462
intron10 10434461 10434168
exon11 10434167 10434104
intron11 10434103 10434012
exon12 10434011 10433940
intron12 10433939 10433840
exon13 10433839 10433686
intron13 10433685 10433562
exon14 10433561 10433336
intron14 10433335 10433257
exon15 10433256 10433201
intron15 10433200 10433126
exon16 10433125 10432978
intron16 10432977 10432886
exon17 10432885 10432786
intron17 10432785 10432634
exon18 10432633 10432589
intron18 10432588 10432498
exon19 10432497 10432359
intron19 10432358 10432272
exon20 10432271 10432176
intron20 10432175 10432094
exon21 10432093 10431995
intron21 10431994 10431902
exon22 10431901 10431801
3utr 10431800 10431516

edr2_bp<- read.table(file = "clipboard", sep = " ", header=T)
edr2_bp

edr2 <- data.frame(x1=edr2_bp$start, x2=edr2_bp$end, 
                    y1=c( 0.15, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2,
                          0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.15), 
                    y2=c(-0.15,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,
                         -0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.15), 
                    segment=c('UTR','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon',
                              'Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon',
                              'Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','UTR'))

mob1
segment start end
5utr 10439901 10439789
exon1 10439788 10439766
intron1 10439765 10439429
exon2 10439428 10439380
intron2 10439379 10439289
exon3 10439288 10439180
intron3 10439179 10439108
exon4 10439107 10439014
intron4 10439013 10438774
exon5 10438773 10438641
intron5 10438640 10438559
exon6 10438558 10438396
intron6 10438395 10438288
exon7 10438287 10438214
3utr 10438213 10438024

mob1_bp<- read.table(file = "clipboard", sep = " ", header=T)

mob1 <- data.frame(x1=mob1_bp$start, x2=mob1_bp$end, 
                    y1=c(-0.25,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.25), 
                    y2=c(-0.55,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.55), 
                    segment=c('UTR','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','UTR'))

cmt2
segment start end
5utr 10414428 10414526
exon1 10414527 10415088
intron1 10415089 10415196
exon2 10415197 10415233
intron2 10415197 10415342
exon3 10415343 10416378
intron3 10416379 10416560
exon4 10416561 10416593
intron4 10416594 10416704
exon5 10416705 10416803
intron5 10416804 10416980
exon6 10416981 10417084
intron6 10417085 10417170
exon7 10417171 10417288
intron7 10417289 10417418
exon8 10417419 10417508
intron8 10417509 10417688
exon9 10417689 10417872
intron9 10417873 10417982
exon10 10417983 10418045
intron10 10418046 10418128
exon11 10418129 10418375
intron11 10418376 10418541
exon12 10418542 10418600
intron12 10418601 10418676
exon13 10418677 10418740
intron13 10418741 10418870
exon14 10418871 10419198
intron14 10419199 10419272
exon15 10419273 10419340
intron15 10419341 10419432
exon16 10419433 10419518
intron16 10419519 10419594
exon17 10419595 10419681
intron17 10419682 10419798
exon18 10419799 10419916
intron18 10419917 10419992
exon19 10419993 10420094
intron19 10420095 10420281
exon20 10420282 10420324
intron20 10420325 10420408
exon21 10420409 10420472
intron21 10420473 10420562
exon22 10420563 10420664
intron22 10420665 10420770
exon23 10420771 10420935
3utr 10420936 10421210

cmt2_bp<- read.table(file = "clipboard", sep = " ", header=T)

cmt2 <- data.frame(x1=cmt2_bp$start, x2=cmt2_bp$end, 
                   y1=c(-0.25,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,
                        -0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.25), 
                   y2=c(-0.55,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,
                        -0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.55), 
                   segment=c('UTR','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon',
                             'Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon',
                             'Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','UTR'))

At4g19010
segment start end
5utr 10414324 10414222
exon1 10414221 10413154
intron1 10413153 10412928
exon2 10412927 10412738
intron2 10412737 10412453
exon3 10412452 10412308
intron3 10412307 10412231
exon4 10412230 10412164
intron4 10412163 10412031
exon5 10412030 10411929
intron5 10411928 10411841
exon6 10411840 10411716
3utr 10411715 10411467

At4g19010_bp<- read.table(file = "clipboard", sep = " ", header=T)

At4g19010 <- data.frame(x1=At4g19010_bp$start, x2=At4g19010_bp$end, 
                   y1=c( 0.15, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.06, 0.2, 0.15), 
                   y2=c(-0.15,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.06,-0.2,-0.15), 
                   segment=c('UTR','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','UTR'))

At4g19006
segment start end
5utr 10411467 10411349
exon1 10411348 10411199
intron1 10411198 10410895
exon2 10410894 10410449
intron2 10410448 10410254
exon3 10410253 10410174
intron3 10410173 10410080
exon4 10410079 10409884
intron4 10409883 10409794
exon5 10409793 10409694
intron5 10409693 10409535
exon6 10409534 10409350
3utr 10409349 10409149

At4g19006_bp<- read.table(file = "clipboard", sep = " ", header=T)

At4g19006 <- data.frame(x1=At4g19006_bp$start, x2=At4g19006_bp$end, 
                        y1=c(-0.25,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.34,-0.2,-0.25), 
                        y2=c(-0.55,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.46,-0.6,-0.55), 
                        segment=c('UTR','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','Intron','Exon','UTR'))

ggplot() + 
  theme_grey() +
  labs(x="\nbp (Chromosome 4)") + 
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(face = 'bold'),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size = 10),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_x_continuous(breaks = seq(10410000,10440000, 2500), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 0.65), expand = c(0,0)) +
  geom_segment(aes(x=c(10423949,10427795,10430319), y=c(-1,-1,-1), xend=c(10423949, 10427795,10430319), yend=c(-0.68,-0.68,-0.68)), size = 1) +
  geom_rect(data=nip, mapping=aes(xmin=x2, xmax=x1, ymin=y2, ymax=y1, fill=segment), color="black", alpha=0.5) +
  geom_rect(data=lcr7, mapping=aes(xmin=x2, xmax=x1, ymin=y2, ymax=y1, fill=segment), color="black", alpha=0.5) +
  geom_rect(data=lcr15, mapping=aes(xmin=x2, xmax=x1, ymin=y2, ymax=y1, fill=segment), color="black", alpha=0.5) +
  geom_rect(data=edr2, mapping=aes(xmin=x2, xmax=x1, ymin=y2, ymax=y1, fill=segment), color="black", alpha=0.5) +
  geom_rect(data=mob1, mapping=aes(xmin=x2, xmax=x1, ymin=y2, ymax=y1, fill=segment), color="black", alpha=0.5) +
  geom_rect(data=cmt2, mapping=aes(xmin=x2, xmax=x1, ymin=y2, ymax=y1, fill=segment), color="black", alpha=0.5) +
  geom_rect(data=At4g19010, mapping=aes(xmin=x2, xmax=x1, ymin=y2, ymax=y1, fill=segment), color="black", alpha=0.5) +
  geom_rect(data=At4g19006, mapping=aes(xmin=x2, xmax=x1, ymin=y2, ymax=y1, fill=segment), color="black", alpha=0.5) +
  #geom_segment(aes(x=10423508, y=0.3, xend=10421520, yend=0.3),
   # lineend = 'butt', linejoin = 'mitre',
    #size = 1, arrow = arrow(length = unit(0.1, "inches"))) +
  geom_text(aes(x=c(10423508,10421580, 10422534,10422534,10424049,10427895,10430419,10428305,10428305,10430495,10430495,10434500,10434500, 
                    10439000,10439000,10418000, 10418000, 10413000, 10410400,10409240,10411447,10411537,10414324,10414448,10421190,
                    10427665,10429000,10429931,10431083,10431576,10437566,10438084,10439881), 
                y=c(0.27,0.27,0.45,0.3,-0.6,-0.6,-0.6,0.45,0.3,0.45,0.3,0.45,0.3,0.05,-0.1,0.05,-0.1,0.32,-0.08,     
                    -0.66,-0.66,0.27,0.27,-0.66,-0.66,-0.25,-0.25,-0.25,-0.25,0.27,0.27,-0.66,-0.66), 
                label=c("5'","3'","At4g19030","NIP1;1","rs2412612","rs2412782","rs2412888","At4g19035","LCR7","At4g19038","LCR15", 
                        "At4g19040","EDR2","At4g19045","MOB1","At4g19020","CMT2","At4g19010","At4g19006",
                        "3'","5'","3'","5'","5'","3'","3'","5'","3'","5'","3'","5'","3'","5'")), size = 3.2) 
  