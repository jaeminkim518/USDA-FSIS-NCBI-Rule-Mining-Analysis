## 04_table4.R
## Table 4: % of isolates with resistance genes for each NARMS class, by year
## From Glass et al.

library(readr)
library(dplyr)
library(tidyr)

setwd("~/USDA/Rule Mining Work")

isolates <- read_csv("Output/isolates_clean.csv",
                     col_select = c(BioSample, specimen_type, sample_type, Year))
isolate_classes <- read_csv("Output/isolate_classes.csv")

df <- isolates %>% left_join(isolate_classes, by = "BioSample")
nrow(df)

narms_classes <- c("AMINOGLYCOSIDE", "BETA-LACTAM", "FOLATE-PATHWAY-INHIBITOR", "MACROLIDE",
                   "PHENICOL", "COLISTIN", "QUINOLONE", "TETRACYCLINE")

df$MDR <- ifelse(df$amr_category == "MDR", 1, 0)

# keep the main sample types
df <- df %>%
  filter((specimen_type != "other" & sample_type %in% c("Beef", "Chicken", "Pork", "Turkey")) |
           (specimen_type == "other" & sample_type == "Clinical"))

table(df$specimen_type, df$sample_type)

table4 <- df %>%
  group_by(specimen_type, sample_type, Year) %>%
  summarise(N = n(),
            across(all_of(c(narms_classes, "MDR")), ~ round(100 * mean(.x), 1))) %>%
  ungroup()

table4$Year <- as.character(table4$Year)

# drop years with fewer than 20 isolates
table4 %>% filter(N < 20)
table4 <- table4 %>% filter(N >= 20)

# Mean and SD across years per sample
mean_rows <- table4 %>%
  group_by(specimen_type, sample_type) %>%
  summarise(Year = "Mean",
            N = sum(N),
            across(all_of(c(narms_classes, "MDR")), ~ round(mean(.x), 1))) %>%
  ungroup()

sd_rows <- table4 %>%
  group_by(specimen_type, sample_type) %>%
  summarise(Year = "SD",
            N = NA,
            across(all_of(c(narms_classes, "MDR")), ~ round(sd(.x), 1))) %>%
  ungroup()

table4_all <- bind_rows(table4, mean_rows, sd_rows) %>%
  arrange(specimen_type, sample_type, Year)

for (sp in unique(table4_all$specimen_type)) {
  t <- table4_all <- bind_rows(table4, mean_rows, sd_rows) %>%
    select(sample_type, specimen_type, everything()) %>%
    arrange(sample_type, specimen_type, Year)
  # One file per sample type 
  for (st in unique(table4_all$sample_type)) {
    t <- table4_all %>% filter(sample_type == st)
    print(st)
    print(t, n = 20)
    write_csv(t, paste0("Output/table4_", st, ".csv"))
  }
}
