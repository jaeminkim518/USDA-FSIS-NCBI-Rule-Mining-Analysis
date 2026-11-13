## 06_3_enteritidis_ssek3.R
## sseK3 in Enteritidis examination

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("~/USDA/Rule Mining Work")

isolates <- read_csv("Output/isolates_clean.csv",
                     col_select = c(BioSample, specimen_type, sample_type, Year, serotype,
                                    Location, SNP_cluster, min_same, min_diff))
hps <- read_csv("Output/hps.csv")

ent <- isolates %>%
  filter(serotype == "Enteritidis") %>%
  filter((specimen_type != "other" & sample_type == "Chicken") |
           (specimen_type == "other" & sample_type == "Clinical")) %>%
  left_join(hps, by = "BioSample")

ent$group <- paste(ent$specimen_type, ent$sample_type)
ent$group[ent$group == "other Clinical"] <- "CDC clinical"

# min_diff = SNP distance to the nearest isolate of a different isolation type
ent$min_diff <- as.numeric(ent$min_diff)
ent$min_same <- as.numeric(ent$min_same)

# 1. How many isolates carry sseK3 per group
ent %>%
  group_by(group) %>%
  summarise(N = n(), n_sseK3 = sum(sseK3), pct_sseK3 = round(100 * mean(sseK3), 1))

# 2. SNP clusters of the sseK3-positive clinical isolates
clin <- ent %>% filter(group == "CDC clinical")

sum(is.na(clin$SNP_cluster))

ssek3_clusters <- clin %>%
  filter(sseK3 == 1, !is.na(SNP_cluster)) %>%
  count(SNP_cluster, sort = TRUE)

head(ssek3_clusters, 20)
write_csv(ssek3_clusters, "Output/enteritidis_ssek3_clusters.csv")

nrow(ssek3_clusters)
sum(head(ssek3_clusters$n, 5)) / sum(ssek3_clusters$n)

# 3. Are sseK3-positive isolates further from non-clinical isolates
clin %>%
  group_by(sseK3) %>%
  summarise(N = n(),
            median_min_diff = median(min_diff, na.rm = TRUE),
            median_min_same = median(min_same, na.rm = TRUE),
            pct_within_20 = round(100 * mean(min_diff <= 20, na.rm = TRUE), 1))

ggplot(clin, aes(x = factor(sseK3), y = min_diff)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = "sseK3 present", y = "SNP distance to nearest non-clinical isolate",
       title = "sseK3-positive clinical Enteritidis closeness to food isolates")

ggsave("Output/Figures/enteritidis_ssek3_min_diff.png", width = 7, height = 5)

# 4. sseK3 over time in clinical
ssek3_year <- clin %>%
  group_by(Year) %>%
  summarise(N = n(), pct_sseK3 = round(100 * mean(sseK3), 1))

ssek3_year
write_csv(ssek3_year, "Output/enteritidis_ssek3_by_year.csv")

# 5. Do any chicken isolates share a cluster with the sseK3-positive clinical isolates
chicken_clusters <- ent %>%
  filter(group != "CDC clinical", !is.na(SNP_cluster)) %>%
  count(SNP_cluster, name = "n_chicken")

shared <- ssek3_clusters %>%
  left_join(chicken_clusters, by = "SNP_cluster") %>%
  arrange(desc(n))

head(shared, 20)
write_csv(shared, "Output/enteritidis_ssek3_shared_clusters.csv")

# 6. US or foreign test on lineages (uses raw file as 01_R takes away non US)
raw <- read_tsv("Data/isolates.tsv.gz", col_types = cols(.default = "c"))
raw <- raw %>% rename(SNP_cluster = `SNP cluster`)
raw$country <- str_remove(raw$Location, ":.*")
raw$country[is.na(raw$country)] <- "unknown"
raw$country[tolower(raw$country) %in% c("not provided", "missing", "not collected", "not applicable")] <- "unknown"

# US share of each cluster with at least 20 sseK3 cases
check <- ssek3_clusters$SNP_cluster[ssek3_clusters$n >= 20]

us_share <- raw %>%
  filter(SNP_cluster %in% check) %>%
  group_by(SNP_cluster) %>%
  summarise(total = n(),
            n_known = sum(country != "unknown"),
            pct_us = round(100 * sum(country == "USA") / n_known, 1)) %>%
  arrange(pct_us)

us_share
write_csv(us_share, "Output/enteritidis_ssek3_us_share.csv")

ggplot(us_share, aes(x = pct_us, y = reorder(SNP_cluster, pct_us))) +
  geom_point(color = "#B4645C") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey50") +
  labs(x = "% of global isolates from the USA", y = NULL,
       title = "SNP clusters carrying sseK3-positive clinical Enteritidis")

ggsave("Output/Figures/enteritidis_ssek3_us_share.png", width = 8, height = 9)

# foreign cluster shares
foreign <- us_share$SNP_cluster[us_share$pct_us < 50]
clin$foreign <- clin$SNP_cluster %in% foreign
clin %>% filter(sseK3 == 1) %>% count(foreign)
ssek3_foreign_year <- clin %>%
  filter(sseK3 == 1) %>%
  group_by(Year) %>%
  summarise(N = n(), pct_foreign = round(100 * mean(foreign), 1))

ssek3_foreign_year

ggplot(ssek3_foreign_year, aes(x = Year, y = pct_foreign)) +
  geom_line(color = "#B4645C") +
  geom_point(color = "#B4645C") +
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  coord_cartesian(ylim = c(0, 25)) +
  labs(y = "% in mostly foreign lineages",
       title = "sseK3-positive clinical Enteritidis in foreign lineages by year")

ggsave("Output/Figures/enteritidis_ssek3_foreign_by_year.png", width = 8, height = 5)