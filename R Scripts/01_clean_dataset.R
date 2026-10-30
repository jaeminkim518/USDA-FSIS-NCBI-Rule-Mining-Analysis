## 01_clean_dataset.R
## Cleaning the NCBI Salmonella isolates data / US only, isolates collected 2015 onward.
## Data downloaded from NCBI Isolates Browser on 10/23/2026

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

setwd("~/USDA/Rule Mining Work")
dir.create("Output", showWarnings = FALSE)

isolates <- read_tsv("Data/isolates.tsv.gz", col_types = cols(.default = "c"))
nrow(isolates)
names(isolates)

isolates <- isolates %>%
  rename(BioProject = `#BioProject`,
         collection_date = `Collection date`,
         source = `Isolation source`,
         isolation_type = `Isolation type`,
         amr_genotypes = `AMR genotypes`,
         vir_genotypes = `Virulence genotypes`,
         computed_types = `Computed types`,
         SNP_cluster = `SNP cluster`,
         min_same = `Min-same`,
         min_diff = `Min-diff`)


# Keep the 4 BioProjects
# PRJNA242847 = FSIS HACCP, PRJNA292666 = FSIS NARMS, PRJNA292661 = FDA NARMS, PRJNA230403 = CDC
isolates <- isolates %>%
  filter(BioProject %in% c("PRJNA242847", "PRJNA292666", "PRJNA292661", "PRJNA230403"))

nrow(isolates)
table(isolates$BioProject)


# Specimen Type
isolates$specimen_type <- NA
isolates$specimen_type[isolates$BioProject == "PRJNA242847"] <- "non-cecal"
isolates$specimen_type[isolates$BioProject == "PRJNA292666"] <- "cecal"
isolates$specimen_type[isolates$BioProject == "PRJNA292661"] <- "retail"
isolates$specimen_type[isolates$BioProject == "PRJNA230403"] <- "other"

table(isolates$specimen_type)


# Sample Type
src <- tolower(isolates$source)
src[is.na(src)] <- ""

isolates$sample_type <- "Other"
isolates$sample_type[str_detect(src, "pork|swine|porcine|\\bhog|\\bsow\\b|\\bpig")] <- "Pork"
isolates$sample_type[str_detect(src, "beef|cattle|bovine|\\bcow|steer|heifer|\\bbull\\b|calf|veal")] <- "Beef"
isolates$sample_type[str_detect(src, "chicken|gallus|broiler")] <- "Chicken"
isolates$sample_type[str_detect(src, "turkey|meleagris")] <- "Turkey"
isolates$sample_type[isolates$isolation_type %in% "clinical"] <- "Clinical"

table(isolates$specimen_type, isolates$sample_type)

isolates %>%
  filter(sample_type == "Other", specimen_type != "other") %>%
  count(source, sort = TRUE) %>%
  head(20)


# Antigen and Serotype
isolates$serotype <- str_remove(isolates$computed_types, ".*serotype=")
isolates$antigen <- str_remove(isolates$computed_types, ",serotype=.*")
isolates$antigen <- str_remove(isolates$antigen, "antigen_formula=")

head(isolates %>% select(computed_types, antigen, serotype))
sort(table(isolates$serotype), decreasing = TRUE)[1:10]


# Year
# collection date is like "2019" or "2019-05-14"
isolates$Year <- as.numeric(substr(isolates$collection_date, 1, 4))

table(isolates$Year, useNA = "ifany")

# US only, 2015 onward
isolates <- isolates %>%
  filter(is.na(Location) | Location == "not provided" | str_starts(Location, "USA")) %>%
  filter(Year >= 2015)

nrow(isolates)
table(isolates$specimen_type, isolates$Year)

# Save
isolates_clean <- isolates %>%
  select(BioProject, Strain, BioSample, specimen_type, sample_type, Year,
         collection_date, Location, source, isolation_type, antigen, serotype,
         SNP_cluster, min_same, min_diff, amr_genotypes, vir_genotypes)

write_csv(isolates_clean, "Output/isolates_clean.csv")