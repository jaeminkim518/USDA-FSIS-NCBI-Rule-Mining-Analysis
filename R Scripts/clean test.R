
## 01_clean_dataset_test.R
## Test version of 01
## Same cleaning but isolates with incomplete serotype calls in an Infantis SNP cluster are relabeled as Infantis.

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

setwd("~/USDA/Rule Mining Work")
dir.create("Output", showWarnings = FALSE)

isolates <- read_tsv("Data/isolates.tsv.gz", col_types = cols(.default = "c"))

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
isolates <- isolates %>%
  filter(BioProject %in% c("PRJNA242847", "PRJNA292666", "PRJNA292661", "PRJNA230403"))


# Specimen Type
isolates$specimen_type <- NA
isolates$specimen_type[isolates$BioProject == "PRJNA242847"] <- "non-cecal"
isolates$specimen_type[isolates$BioProject == "PRJNA292666"] <- "cecal"
isolates$specimen_type[isolates$BioProject == "PRJNA292661"] <- "retail"
isolates$specimen_type[isolates$BioProject == "PRJNA230403"] <- "other"


# Sample Type
src <- tolower(isolates$source)
src[is.na(src)] <- ""

isolates$sample_type <- "Other"
isolates$sample_type[str_detect(src, "pork|swine|porcine|\\bhog|\\bsow\\b|\\bpig")] <- "Pork"
isolates$sample_type[str_detect(src, "beef|cattle|bovine|\\bcow|steer|heifer|\\bbull\\b|calf|veal")] <- "Beef"
isolates$sample_type[str_detect(src, "chicken|gallus|broiler")] <- "Chicken"
isolates$sample_type[str_detect(src, "turkey|meleagris")] <- "Turkey"
isolates$sample_type[isolates$isolation_type %in% "clinical"] <- "Clinical"


# Antigen and Serotype
isolates$serotype <- str_remove(isolates$computed_types, ".*serotype=")
isolates$antigen <- str_remove(isolates$computed_types, ",serotype=.*")
isolates$antigen <- str_remove(isolates$antigen, "antigen_formula=")


# Year
isolates$Year <- as.numeric(substr(isolates$collection_date, 1, 4))

# US only, 2015 onward
isolates <- isolates %>%
  filter(is.na(Location) | Location == "not provided" | str_starts(Location, "USA")) %>%
  filter(Year >= 2015)

nrow(isolates)


# Incomplete serotype calls -> Infantis
isolates$cluster_id <- str_remove(isolates$SNP_cluster, "\\..*")

cluster_main <- isolates %>%
  filter(!is.na(cluster_id)) %>%
  count(cluster_id, serotype) %>%
  group_by(cluster_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cluster_id, cluster_serotype = serotype)

isolates <- isolates %>% left_join(cluster_main, by = "cluster_id")

partial_infantis <- c("I -:r:1,5", "I -:-:-", "- 7:r:1,5", "I 7:r:-", "I 7:-:1,5",
                      "- -:r:1,5", "I -:-:1,5")

isolates$serotype_raw <- isolates$serotype
fix <- isolates$serotype %in% partial_infantis & isolates$cluster_serotype %in% "Infantis"

sum(fix)
table(isolates$serotype_raw[fix], isolates$specimen_type[fix])

isolates$serotype[fix] <- "Infantis"


# Save
isolates_clean <- isolates %>%
  select(BioProject, Strain, BioSample, specimen_type, sample_type, Year,
         collection_date, Location, source, isolation_type, antigen, serotype, serotype_raw,
         SNP_cluster, min_same, min_diff, amr_genotypes, vir_genotypes)

write_csv(isolates_clean, "Output/isolates_clean_test.csv")