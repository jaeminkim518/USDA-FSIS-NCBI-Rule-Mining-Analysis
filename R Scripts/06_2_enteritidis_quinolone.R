## 06_2_enteritidis_quinolone.R
## Quinolone resistance determinants in quinolone-resistant Enteritidis

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/USDA/Rule Mining Work")

isolates <- read_csv("Output/isolates_clean.csv",
                     col_select = c(BioSample, specimen_type, sample_type, Year, serotype))
isolate_classes <- read_csv("Output/isolate_classes.csv")
amr_long <- read_csv("Output/amr_genes_long.csv")
gene_map <- read_csv("Output/gene_class_map.csv")

# Genes and mutations in the quinolone class
qui_genes <- gene_map %>% filter(grepl("QUINOLONE", class))
qui_genes$gene

# Quinolone-resistant SE from chicken and CDC clinical
ent <- isolates %>%
  filter(serotype == "Enteritidis") %>%
  filter((specimen_type != "other" & sample_type == "Chicken") |
         (specimen_type == "other" & sample_type == "Clinical")) %>%
  left_join(isolate_classes, by = "BioSample") %>%
  filter(QUINOLONE == 1)

ent$group <- paste(ent$specimen_type, ent$sample_type)
ent$group[ent$group == "other Clinical"] <- "CDC clinical"

d87y <- amr_long %>% filter(gene == "gyrA_D87Y")
ent$gyrA_D87Y <- ifelse(ent$BioSample %in% d87y$BioSample, 1, 0)

qui_long <- amr_long %>%
  filter(BioSample %in% ent$BioSample, gene %in% qui_genes$gene)

profiles <- qui_long %>%
  group_by(BioSample) %>%
  summarise(profile = paste(sort(gene), collapse = " + "))

ent <- ent %>% left_join(profiles, by = "BioSample")
sum(is.na(ent$profile))

# 1. Commonness of each profile
profile_table <- ent %>%
  count(group, profile) %>%
  group_by(group) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(group, desc(n))

print(profile_table, n = 50)
write_csv(profile_table, "Output/enteritidis_quinolone_profiles.csv")

# 2. What genes do clinical isolates without gyrA D87Y have
no_d87y <- ent %>% filter(group == "CDC clinical", gyrA_D87Y == 0)
nrow(no_d87y)

qui_long %>%
  filter(BioSample %in% no_d87y$BioSample) %>%
  count(gene, sort = TRUE)

# 3. % of quinolone-resistant isolates with gyrA D87Y, by year
d87y_year <- ent %>%
  group_by(group, Year) %>%
  summarise(N = n(), pct_d87y = round(100 * mean(gyrA_D87Y), 1)) %>%
  filter(N >= 20)

print(d87y_year, n = 50)
write_csv(d87y_year, "Output/enteritidis_d87y_share.csv")

# 4. Quinolone resistance split into w/ and w/out gyrA D87Y, as % of all SE
ent_all <- isolates %>%
  filter(serotype == "Enteritidis") %>%
  filter((specimen_type != "other" & sample_type == "Chicken") |
           (specimen_type == "other" & sample_type == "Clinical")) %>%
  left_join(isolate_classes, by = "BioSample")

ent_all$group <- paste(ent_all$specimen_type, ent_all$sample_type)
ent_all$group[ent_all$group == "other Clinical"] <- "CDC clinical"

ent_all$gyrA_D87Y <- ifelse(ent_all$BioSample %in% d87y$BioSample, 1, 0)
ent_all$QUI_other <- ifelse(ent_all$QUINOLONE == 1 & ent_all$gyrA_D87Y == 0, 1, 0)

qui_split <- ent_all %>%
  group_by(group, Year) %>%
  summarise(N = n(),
            `With gyrA D87Y` = 100 * mean(gyrA_D87Y),
            `Without gyrA D87Y` = 100 * mean(QUI_other))


split_long <- qui_split %>%
  pivot_longer(cols = c(`With gyrA D87Y`, `Without gyrA D87Y`),
               names_to = "type", values_to = "percent")

ggplot(split_long, aes(x = Year, y = percent, fill = type)) +
  geom_col() +
  facet_wrap(~ group) +
  scale_fill_manual(values = c("#B4645C", "#9AA0A6")) +
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  labs(y = "% of Enteritidis isolates", fill = "Quinolone resistant",
       title = "Quinolone resistance in Enteritidis, with and without gyrA D87Y")

ggsave("Output/Figures/enteritidis_quinolone_split.png", width = 10, height = 6)