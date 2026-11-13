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
ggsave("Output/Figures/infantis_by_year.png", width = 9, height = 5)

# MDR within Infantis in chicken/clinical only by year
infantis <- isolates %>% filter(serotype == "Infantis")
infantis <- infantis %>%
  filter((specimen_type != "other" & sample_type == "Chicken") |
           (specimen_type == "other" & sample_type == "Clinical"))
infantis <- infantis %>% left_join(isolate_classes, by = "BioSample")

infantis$MDR <- ifelse(infantis$amr_category == "MDR", 1, 0)

infantis_mdr <- infantis %>%
  group_by(specimen_type, Year) %>%
  summarise(N = n(), pct_MDR = 100 * mean(MDR))

# skip years with fewer than 20 isolates
infantis_mdr <- infantis_mdr %>% filter(N >= 20)

infantis_mdr
write_csv(infantis_mdr, "Output/infantis_mdr_by_year.csv")

ggplot(infantis_mdr, aes(x = Year, y = pct_MDR, color = specimen_type)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = 2015:2026) +   # show every year on the x-axis
  labs(y = "% MDR among Infantis", title = "MDR within Salmonella Infantis by year")

ggsave("Output/Figures/infantis_mdr_by_year.png", width = 8, height = 5)

# AMR within Infantis by year (chicken stages + CDC clinical)
inf <- isolates %>%
  filter(serotype == "Infantis") %>%
  filter((specimen_type != "other" & sample_type == "Chicken") |
           (specimen_type == "other" & sample_type == "Clinical")) %>%
  left_join(isolate_classes, by = "BioSample") %>%
  left_join(hps, by = "BioSample")

inf$group <- paste(inf$specimen_type, inf$sample_type)
inf$group[inf$group == "other Clinical"] <- "CDC clinical"

inf$MDR <- ifelse(inf$amr_category == "MDR", 1, 0)

d87y <- amr_long %>% filter(gene == "gyrA_D87Y")
inf$gyrA_D87Y <- ifelse(inf$BioSample %in% d87y$BioSample, 1, 0)

# Other = resistant to at least one of the other classes
inf$OTHER <- ifelse(inf$AMINOGLYCOSIDE == 1 | inf$`FOLATE-PATHWAY-INHIBITOR` == 1 |
                      inf$TETRACYCLINE == 1 | inf$MACROLIDE == 1 | inf$COLISTIN == 1, 1, 0)

inf_year <- inf %>%
  group_by(group, Year) %>%
  summarise(N = n(),
            QUINOLONE = 100 * mean(QUINOLONE),
            `BETA-LACTAM` = 100 * mean(`BETA-LACTAM`),
            PHENICOL = 100 * mean(PHENICOL),
            OTHER = 100 * mean(OTHER),
            gyrA_D87Y = 100 * mean(gyrA_D87Y),
            MDR = 100 * mean(MDR))

# skip years with fewer than 20 Infantis isolates
inf_year <- inf_year %>% filter(N >= 20)
inf_year

# check the 5 virulence genes prevalence
inf %>%
  group_by(group) %>%
  summarise(invA = 100 * mean(invA), avrA = 100 * mean(avrA), lpfB = 100 * mean(lpfB),
            sspH2 = 100 * mean(sspH2), sseK2 = 100 * mean(sseK2))

inf_plot <- inf_year %>%
  pivot_longer(cols = c(QUINOLONE, `BETA-LACTAM`, PHENICOL, OTHER, gyrA_D87Y, MDR),
               names_to = "measure", values_to = "percent")

inf_levels <- c("QUINOLONE", "BETA-LACTAM", "PHENICOL", "OTHER", "gyrA_D87Y", "MDR")
inf_labels <- c("Class: quinolone", "Class: beta-lactam", "Class: phenicol",
                "Class: other", "Mutation: gyrA D87Y", "Summary: MDR (3+ classes)")

inf_plot$measure <- factor(inf_plot$measure, levels = inf_levels, labels = inf_labels)

inf_colors <- c("Class: quinolone" = "#B4645C", "Class: beta-lactam" = "#C8A165",
                "Class: phenicol" = "#A0786A", "Class: other" = "#9AA0A6",
                "Mutation: gyrA D87Y" = "black", "Summary: MDR (3+ classes)" = "black")

inf_lines <- c(rep("solid", 4), "dashed", "dotted")
names(inf_lines) <- inf_labels

ggplot(inf_plot, aes(x = Year, y = percent, color = measure, linetype = measure)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.1) +
  facet_wrap(~ group) +
  scale_color_manual(values = inf_colors) +
  scale_linetype_manual(values = inf_lines) +
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  labs(y = "% of Infantis isolates", color = NULL, linetype = NULL,
       title = "AMR in Salmonella Infantis by year",
       subtitle = "Virulence: invA, avrA, lpfB, sspH2, sseK2 present in nearly all isolates",
       caption = "Class: other = aminoglycoside, folate pathway inhibitor, tetracycline, macrolide or colistin")

ggsave("Output/Figures/infantis_amr_by_year.png", width = 10, height = 6)

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