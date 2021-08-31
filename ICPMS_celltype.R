library("readxl")
ICPMS_celltype <- read_excel("C:/Users/Todd/Desktop/Todd Lab-Stuff/ICPMS/CellType_ICPMS_Arsenic_Jan_9_2018.xlsx")

ICPMS_celltype <- ICPMS_celltype %>% 
  dplyr::rename(sample = '...1') %>% 
  filter(!is.na(sample)) %>% 
  dplyr::na_if('BDL') %>%
  mutate_at(vars(-sample), funs(as.double(.))) %>% 
  purrr::modify_if(is.numeric, ~round(., 2)) %>% 
  mutate_if(is.character, funs(str_replace(., "\\d", ""))) %>% 
  mutate_if(is.character, funs(str_replace(., '[+]', " 50"))) %>%
  mutate_if(is.character, funs(str_replace(., '[-]', " 0"))) %>%
  separate(sample, into= c('celltype', 'condition'), remove = F)
View(ICPMS_celltype)

As_sd <- group_by(ICPMS_celltype, sample) %>%
  summarise(count = n(),
            element = 'As',
            mean = mean(As, na.rm = TRUE),
            sd = sd(As, na.rm = TRUE), 
            se = sd/sqrt(count)) %>%
  ungroup() %>% 
  separate(sample, into= c('celltype', 'condition'), remove = F)

Si_sd <- group_by(ICPMS_celltype, sample) %>%
  summarise(count = n(),
            element = 'Si',
            mean = mean(Si, na.rm = TRUE),
            sd = sd(Si, na.rm = TRUE), 
            se = sd/sqrt(count)) %>%
  ungroup() %>% 
  separate(sample, into= c('celltype', 'condition'), remove = F)

S_sd <- group_by(ICPMS_celltype, sample) %>%
  summarise(count = n(),
            element = 'S',
            mean = mean(S, na.rm = TRUE),
            sd = sd(S, na.rm = TRUE), 
            se = sd/sqrt(count)) %>%
  ungroup() %>% 
  separate(sample, into= c('celltype', 'condition'), remove = F)

View(Si_sd)
As_Si_S <- bind_rows(As_sd, Si_sd, S_sd)
View(As_Si_S)

Si_bar <- ggplot(Si_sd, aes(x=sample, y= mean, fill= condition, group = celltype)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2) +
  theme_light(base_size = 12) +
  scale_fill_manual(values=c('cyan3','firebrick2')) +
  scale_x_discrete(labels=c("Cortex", "Cortex", "Endoderm", "Endoderm", "Epiderm", "Epiderm")) +
  scale_y_continuous(expand = c(0, 0.5)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 12 ),
        axis.text.y = element_text(size = 14),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank()) +
  labs(title = "Si",
       subtitle = "per 50K sorted cells")

S_bar <- ggplot(S_sd, aes(x=sample, y= mean, fill= condition, group = celltype)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2) +
  theme_light(base_size = 12) +
  scale_fill_manual(values=c('cyan3','firebrick2')) +
  scale_x_discrete(labels=c("Cortex", "Cortex", "Endoderm", "Endoderm", "Epiderm", "Epiderm")) +
  scale_y_continuous(expand = c(0, 0.5)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 12 ),
        axis.text.y = element_text(size = 14),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "none") +
  labs(title = "S",
       subtitle = "per 50K sorted cells")

As_bar <- ggplot(As_sd, aes(x=sample, y= mean, fill= condition, group = celltype)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2) +
  theme_light(base_size = 12) +
  scale_fill_manual(values=c('cyan3','firebrick2')) +
  scale_x_discrete(labels=c("Cortex", "Cortex", "Endoderm", "Endoderm", "Epiderm", "Epiderm")) +
  scale_y_continuous(expand = c(0, 0.5)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 12 ),
        axis.text.y = element_text(size = 14),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = 14, face = 'bold'),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "none") +
  labs(title = "As",
       subtitle = "per 50K sorted cells", 
       y = 'Total (ng)\n')

As_bar+S_bar+Si_bar

ggsave('ICPMS_celltype_As_Si_S.tiff', plot = last_plot(), dpi = 600)
View(ICPMS_celltype)
res.aov <- aov(As ~ sample, data = ICPMS_celltype)
summary(res.aov)
TukeyHSD(res.aov)

# Tukey_a <- WER0,Cortex0,SCR0, Cortex50
# Tukey_b <- SCR50
# Tukey_c <- WER50
