## 05_2_serotype_composition_mdr.R
## Serotypes making up MDR in chicken

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/USDA/Rule Mining Work")

isolates <- read_csv("Output/isolates_clean.csv",
                     col_select = c(BioSample, specimen_type, sample_type, Year, serotype))
isolate_classes <- read_csv("Output/isolate_classes.csv")

narms_classes <- c("AMINOGLYCOSIDE", "BETA-LACTAM", "FOLATE-PATHWAY-INHIBITOR", "MACROLIDE",
                   "PHENICOL", "COLISTIN", "QUINOLONE", "TETRACYCLINE")
short_names <- c("AMI", "BLA", "FPI", "MAC", "PHE", "COL", "QUI", "TET")

chicken <- isolates %>%
  filter(sample_type == "Chicken", specimen_type != "other") %>%
  left_join(isolate_classes, by = "BioSample")

top_serotypes <- c("Infantis", "Kentucky", "Enteritidis", "Typhimurium", "Schwarzengrund")
chicken$serotype_group <- ifelse(chicken$serotype %in% top_serotypes, chicken$serotype, "Other")

chicken$pattern <- apply(chicken[, narms_classes], 1, function(x) {
  paste(short_names[x == 1], collapse = "; ")
})

# 1. Share of chicken MDR isolates by serotype per year
mdr_share <- chicken %>%
  filter(amr_category == "MDR") %>%
  count(Year, serotype_group) %>%
  group_by(Year) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

mdr_share %>%
  select(-n) %>%
  pivot_wider(names_from = Year, values_from = pct, values_fill = 0)

write_csv(mdr_share, "Output/chicken_mdr_share_by_serotype.csv")

# 2. Resistance of serotypes
serotype_amr <- chicken %>%
  group_by(serotype_group) %>%
  summarise(N = n(),
            pct_2plus = round(100 * mean(n_classes >= 2), 1),
            pct_MDR = round(100 * mean(amr_category == "MDR"), 1),
            n_patterns = n_distinct(pattern[n_classes >= 2])) %>%
  arrange(desc(pct_MDR))

serotype_amr
write_csv(serotype_amr, "Output/chicken_serotype_amr.csv")

# most common pattern within each serotype
chicken %>%
  filter(n_classes >= 2) %>%
  count(serotype_group, pattern) %>%
  group_by(serotype_group) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  slice_max(n, n = 1) %>%
  ungroup()

# plotting chicken serotype share
mdr_share <- read_csv("Output/chicken_mdr_share_by_serotype.csv")

mdr_share$serotype_group <- factor(mdr_share$serotype_group,
                                levels = c("Infantis", "Typhimurium", "Kentucky", 
                                           "Enteritidis", "Schwarzengrund", "Other"))

ggplot(mdr_share, aes(x = Year, y = pct, fill = serotype_group)) +
  geom_col() + scale_fill_manual(values = c("Infantis" = "#B4645C", "Typhimurium" = "#C8A165",
                               "Kentucky" = "#6FA287", "Enteritidis" = "#5B7FA6",
                               "Schwarzengrund" = "#8E7CC3", "Other" = "#9AA0A6")) + 
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  labs(y = "% of MDR chicken isolates", fill = "Serotype", title = "Serotypes making up MDR in chicken")

ggsave("Output/Figures/chicken_mdr_share_by_serotype.png", width = 8, height = 5)