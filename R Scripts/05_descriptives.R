## 05_descriptives.R
## Descriptive summaries: isolate counts, serotypes, AMR gene prevalence / HPS frequency, unique AMR patterns

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/USDA/Rule Mining Work")

isolates <- read_csv("Output/isolates_clean.csv",
                     col_select = c(BioSample, specimen_type, sample_type, Year, serotype))
amr_long <- read_csv("Output/amr_genes_long.csv")
hps <- read_csv("Output/hps.csv")
isolate_classes <- read_csv("Output/isolate_classes.csv")

counts <- isolates %>%
  count(specimen_type, sample_type, Year) %>%
  pivot_wider(names_from = Year, values_from = n, values_fill = 0)

counts
write_csv(counts, "Output/counts_specimen_sample_year.csv")


# Serotypes
serotype_counts <- isolates %>%
  count(specimen_type, sample_type, serotype) %>%
  group_by(specimen_type, sample_type) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  slice_max(n, n = 10) %>%
  ungroup()

serotype_counts
write_csv(serotype_counts, "Output/top_serotypes.csv")

# % of each top serotype by year
serotype_year <- isolates %>%
  count(specimen_type, sample_type, Year, serotype) %>%
  group_by(specimen_type, sample_type, Year) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  semi_join(serotype_counts, by = c("specimen_type", "sample_type", "serotype")) %>%
  select(-n) %>%
  arrange(Year) %>%
  pivot_wider(names_from = Year, values_from = pct, values_fill = 0) %>%
  arrange(specimen_type, sample_type)

write_csv(serotype_year, "Output/serotypes_by_year.csv")

# Infantis over time
infantis <- isolates %>%
  filter(sample_type %in% c("Beef", "Chicken", "Pork", "Turkey", "Clinical")) %>%
  group_by(specimen_type, sample_type, Year) %>%
  summarise(N = n(), pct_infantis = round(100 * mean(serotype %in% "Infantis"), 1)) %>%
  ungroup()

infantis %>% filter(sample_type == "Chicken")

ggplot(infantis, aes(x = Year, y = pct_infantis, color = specimen_type)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ sample_type) +
  labs(y = "% Infantis", title = "Salmonella Infantis by year")
ggsave("Output/infantis_by_year.png", width = 9, height = 5)


# AMR gene prevalence
gene_prev <- amr_long %>%
  left_join(isolates, by = "BioSample") %>%
  count(specimen_type, sample_type, gene)

group_sizes <- isolates %>% count(specimen_type, sample_type, name = "N")

gene_prev <- gene_prev %>%
  left_join(group_sizes, by = c("specimen_type", "sample_type")) %>%
  mutate(pct = round(100 * n / N, 2),
         point_mutation = grepl("_", gene)) %>%
  arrange(specimen_type, sample_type, desc(pct))

head(gene_prev, 20)
write_csv(gene_prev, "Output/amr_gene_prevalence.csv")

# point mutations
gene_prev %>% filter(point_mutation) %>% head(20)


# HPS frequency
hps_freq <- hps %>%
  left_join(isolates, by = "BioSample") %>%
  group_by(specimen_type, sample_type) %>%
  summarise(N = n(), n_HPS = sum(HPS == "HPS"), pct_HPS = round(100 * mean(HPS == "HPS"), 1)) %>%
  ungroup()

hps_freq
write_csv(hps_freq, "Output/hps_frequency.csv")

hps_year <- hps %>%
  left_join(isolates, by = "BioSample") %>%
  group_by(specimen_type, sample_type, Year) %>%
  summarise(pct_HPS = round(100 * mean(HPS == "HPS"), 1)) %>%
  ungroup() %>%
  arrange(Year) %>%
  pivot_wider(names_from = Year, values_from = pct_HPS) %>%
  arrange(specimen_type, sample_type)

write_csv(hps_year, "Output/hps_by_year.csv")


# Unique AMR patterns
narms_classes <- c("AMINOGLYCOSIDE", "BETA-LACTAM", "FOLATE-PATHWAY-INHIBITOR", "MACROLIDE",
                   "PHENICOL", "COLISTIN", "QUINOLONE", "TETRACYCLINE")
short_names <- c("AMI", "BLA", "FPI", "MAC", "PHE", "COL", "QUI", "TET")

isolate_classes$pattern <- apply(isolate_classes[, narms_classes], 1, function(x) {
  paste(short_names[x == 1], collapse = "; ")
})

patterns <- isolate_classes %>%
  filter(n_classes >= 2) %>%
  left_join(isolates, by = "BioSample") %>%
  count(pattern, specimen_type, sample_type, name = "frequency") %>%
  arrange(desc(frequency))

head(patterns, 20)
write_csv(patterns, "Output/amr_patterns.csv")

patterns %>% count(specimen_type, sample_type, name = "n_unique_patterns")
