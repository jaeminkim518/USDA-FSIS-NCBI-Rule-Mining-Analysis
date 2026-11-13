## 06_enteritidis.R
## Salmonella Enteritidis over time in chicken

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/USDA/Rule Mining Work")

isolates <- read_csv("Output/isolates_clean.csv",
                     col_select = c(BioSample, specimen_type, sample_type, Year, serotype, source))
isolate_classes <- read_csv("Output/isolate_classes.csv")
hps <- read_csv("Output/hps.csv")
amr_long <- read_csv("Output/amr_genes_long.csv")

narms_classes <- c("AMINOGLYCOSIDE", "BETA-LACTAM", "FOLATE-PATHWAY-INHIBITOR", "MACROLIDE",
                   "PHENICOL", "COLISTIN", "QUINOLONE", "TETRACYCLINE")
hps_genes <- c("invA", "avrA", "gtgA", "lpfB", "sspH2", "sseK2", "sseK3", "spvD")

ent <- isolates %>%
  filter(serotype == "Enteritidis") %>%
  filter((specimen_type != "other" & sample_type == "Chicken") |
         (specimen_type == "other" & sample_type == "Clinical"))

ent <- ent %>%
  left_join(isolate_classes, by = "BioSample") %>%
  left_join(hps, by = "BioSample")

ent$group <- paste(ent$specimen_type, ent$sample_type)
ent$group[ent$group == "other Clinical"] <- "CDC clinical"

ent$MDR <- ifelse(ent$amr_category == "MDR", 1, 0)

# pulled out as gyrA D87Y is the main quinolone mutation in Enteritidis
d87y <- amr_long %>% filter(gene == "gyrA_D87Y")
ent$gyrA_D87Y <- ifelse(ent$BioSample %in% d87y$BioSample, 1, 0)

table(ent$group, ent$Year)

ent_year <- ent %>%
  group_by(group, Year) %>%
  summarise(N = n(),
            across(c(all_of(narms_classes), gyrA_D87Y, MDR, all_of(hps_genes)),
                   ~ round(100 * mean(.x), 1))) %>%
  ungroup()

# skip years with fewer than 20 isolates
ent_year <- ent_year %>% filter(N >= 20)
ent_year
write_csv(ent_year, "Output/enteritidis_by_year.csv")


# 1. AMR by year
amr_plot <- ent_year %>%
  pivot_longer(cols = c(all_of(narms_classes), gyrA_D87Y, MDR),
               names_to = "measure", values_to = "percent")

amr_levels <- c("QUINOLONE", "BETA-LACTAM", "TETRACYCLINE", "AMINOGLYCOSIDE",
                "FOLATE-PATHWAY-INHIBITOR", "PHENICOL", "MACROLIDE", "COLISTIN",
                "gyrA_D87Y", "MDR")
amr_labels <- c("Class: quinolone", "Class: beta-lactam", "Class: tetracycline",
                "Class: aminoglycoside", "Class: folate pathway inhibitor",
                "Class: phenicol", "Class: macrolide", "Class: colistin",
                "Mutation: gyrA D87Y", "Summary: MDR (3+ classes)")

amr_plot$measure <- factor(amr_plot$measure, levels = amr_levels, labels = amr_labels)

amr_colors <- c("Class: quinolone" = "#B4645C", "Class: beta-lactam" = "#C8A165",
                "Class: tetracycline" = "#8E7CC3", "Class: aminoglycoside" = "#6FA287",
                "Class: folate pathway inhibitor" = "#5B7FA6", "Class: phenicol" = "#A0786A",
                "Class: macrolide" = "#C98BB0", "Class: colistin" = "#9AA0A6",
                "Mutation: gyrA D87Y" = "black", "Summary: MDR (3+ classes)" = "black")

amr_lines <- c(rep("solid", 8), "dashed", "dotted")
names(amr_lines) <- amr_labels

ggplot(amr_plot, aes(x = Year, y = percent, color = measure, linetype = measure)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.1) +
  facet_wrap(~ group) +
  scale_color_manual(values = amr_colors) +
  scale_linetype_manual(values = amr_lines) +
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  labs(y = "% of Enteritidis isolates", color = NULL, linetype = NULL,
       title = "AMR in Salmonella Enteritidis by year")

ggsave("Output/Figures/enteritidis_amr_by_year.png", width = 10, height = 6)

# 2. HPS genes by year
hps_plot <- ent_year %>%
  pivot_longer(cols = all_of(hps_genes), names_to = "gene", values_to = "percent")

ggplot(hps_plot, aes(x = Year, y = percent, color = gene)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ group) +
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  labs(y = "% of Enteritidis isolates", title = "HPS genes in Salmonella Enteritidis by year")

ggsave("Output/Figures/enteritidis_hps_genes_by_year.png", width = 10, height = 6)

# 3. AMR next to virulence
sync <- ent_year %>%
  pivot_longer(cols = c(QUINOLONE, gyrA_D87Y, spvD, sseK3),
               names_to = "measure", values_to = "percent")

sync$measure <- factor(sync$measure,
                       levels = c("QUINOLONE", "gyrA_D87Y", "spvD", "sseK3"),
                       labels = c("AMR: quinolone", "AMR: gyrA D87Y",
                                  "Virulence: spvD", "Virulence: sseK3"))

ggplot(sync, aes(x = Year, y = percent, color = measure, linetype = measure)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ group) +
  scale_color_manual(values = c("AMR: quinolone" = "#B4645C", "AMR: gyrA D87Y" = "#B4645C",
                                "Virulence: spvD" = "#5B7FA6", "Virulence: sseK3" = "#5B7FA6")) +
  scale_linetype_manual(values = c("AMR: quinolone" = "solid", "AMR: gyrA D87Y" = "dashed",
                                   "Virulence: spvD" = "solid", "Virulence: sseK3" = "dashed")) +
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  labs(y = "% of Enteritidis isolates", color = NULL, linetype = NULL,
       title = "AMR vs virulence genes in Salmonella Enteritidis")

ggsave("Output/Figures/enteritidis_amr_vs_virulence.png", width = 10, height = 6)

# 4. spvB vs spvD
spvb <- read_csv("Output/vir_matrix.csv.gz", col_select = c(BioSample, spvB))
ent <- left_join(ent, spvb, by = "BioSample")

table(ent$spvD, ent$spvB)

spv_year <- ent %>%
  group_by(group, Year) %>%
  summarise(N = n(),
            spvD = 100 * mean(spvD),
            spvB = 100 * mean(spvB))

spv_year <- spv_year %>% filter(N >= 20)
spv_year
write_csv(spv_year, "Output/enteritidis_spv_by_year.csv")

spv_long <- spv_year %>%
  pivot_longer(cols = c(spvD, spvB), names_to = "gene", values_to = "percent")

ggplot(spv_long, aes(x = Year, y = percent, color = gene)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ group) +
  scale_color_manual(values = c("spvD" = "#5B7FA6", "spvB" = "#B4645C")) +
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  labs(y = "% of Enteritidis isolates",
       title = "spvD and spvB in Salmonella Enteritidis by year")

ggsave("Output/Figures/enteritidis_spv_by_year.png", width = 10, height = 6)

# 5. spvD dip in non-cecal chicken
nc <- ent %>% filter(group == "non-cecal Chicken")
nc %>% count(source, sort = TRUE)
table(nc$source, nc$Year)

spvd_source <- nc %>%
  group_by(source, Year) %>%
  summarise(N = n(), spvD = 100 * mean(spvD))

spvd_source <- spvd_source %>% filter(N >= 20)
spvd_source
write_csv(spvd_source, "Output/enteritidis_spvd_by_source.csv")

# shorter names, ordered along the production chain
spvd_source$source_short <- spvd_source$source
spvd_source$source_short[spvd_source$source == "chicken - young chicken carcass rinse (pre-evisceration)"] <- "Rinse (pre-evisc)"
spvd_source$source_short[spvd_source$source == "chicken carcass"] <- "Carcass"
spvd_source$source_short[spvd_source$source == "chicken - young chicken carcass rinse (post-chill)"] <- "Rinse (post-chill)"
spvd_source$source_short[spvd_source$source == "raw intact chicken"] <- "Raw intact"
spvd_source$source_short[spvd_source$source == "Nonintact Chicken"] <- "Nonintact"
spvd_source$source_short[spvd_source$source == "comminuted chicken"] <- "Comminuted"

source_levels <- c("Rinse (pre-evisc)", "Carcass", "Rinse (post-chill)",
                   "Raw intact", "Nonintact", "Comminuted")

spvd_source$source_short <- factor(spvd_source$source_short, levels = source_levels)

source_colors <- c("Rinse (pre-evisc)" = "#C8A165", "Carcass" = "#B4645C",
                   "Rinse (post-chill)" = "#A0786A", "Raw intact" = "#6FA287",
                   "Nonintact" = "#5B7FA6", "Comminuted" = "#8E7CC3")

ggplot(spvd_source, aes(x = Year, y = spvD, color = source_short)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = source_colors) +
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  labs(y = "% with spvD", color = "Isolation source",
       title = "spvD in non-cecal chicken Enteritidis by isolation source",
       caption = "Source categories as submitted to NCBI / carcass and post-chill rinse categories change around 2022")

ggsave("Output/Figures/enteritidis_spvd_by_source.png", width = 10, height = 6)