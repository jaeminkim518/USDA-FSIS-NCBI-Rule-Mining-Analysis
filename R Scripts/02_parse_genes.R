## 02_parse_genes.R
## Parse AMR and virulence genes, unique gene lists, presence/absence matrices
## Present = any status except MISTRANSLATION

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

setwd("~/USDA/Rule Mining Work")

isolates <- read_csv("Output/isolates_clean.csv", col_types = cols(.default = "c"))
nrow(isolates)


# AMR genes
amr_long <- isolates %>%
  select(BioSample, amr_genotypes) %>%
  filter(!is.na(amr_genotypes)) %>%
  separate_rows(amr_genotypes, sep = ",")

# gene is before the last "=", status is after
amr_long$gene <- str_remove(amr_long$amr_genotypes, "=[^=]*$")
amr_long$status <- str_extract(amr_long$amr_genotypes, "[^=]*$")

table(amr_long$status)

amr_long <- amr_long %>%
  filter(status != "MISTRANSLATION") %>%
  distinct(BioSample, gene)

amr_gene_list <- amr_long %>% count(gene, sort = TRUE)
nrow(amr_gene_list)
head(amr_gene_list, 20)

write_csv(amr_gene_list, "Output/amr_gene_list.csv")
write_csv(amr_long, "Output/amr_genes_long.csv")


# Virulence genes
vir_long <- isolates %>%
  select(BioSample, vir_genotypes) %>%
  filter(!is.na(vir_genotypes)) %>%
  separate_rows(vir_genotypes, sep = ",")

vir_long$gene <- str_remove(vir_long$vir_genotypes, "=[^=]*$")
vir_long$status <- str_extract(vir_long$vir_genotypes, "[^=]*$")

table(vir_long$status)

vir_long <- vir_long %>%
  filter(status != "MISTRANSLATION") %>%
  distinct(BioSample, gene)

vir_gene_list <- vir_long %>% count(gene, sort = TRUE)
nrow(vir_gene_list)
head(vir_gene_list, 20)

write_csv(vir_gene_list, "Output/vir_gene_list.csv")


# Presence/absence matrices (1 = present, 0 = absent)
amr_matrix <- amr_long %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = gene, values_from = present, values_fill = 0L)

amr_matrix <- isolates %>%
  select(BioSample) %>%
  left_join(amr_matrix, by = "BioSample") %>%
  mutate(across(-BioSample, ~ replace_na(.x, 0L)))

dim(amr_matrix)
write_csv(amr_matrix, "Output/amr_matrix.csv.gz")


vir_matrix <- vir_long %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = gene, values_from = present, values_fill = 0L)

vir_matrix <- isolates %>%
  select(BioSample) %>%
  left_join(vir_matrix, by = "BioSample") %>%
  mutate(across(-BioSample, ~ replace_na(.x, 0L)))

dim(vir_matrix)
write_csv(vir_matrix, "Output/vir_matrix.csv.gz")

# HPS = at least 5 of these 8 virulence genes present
hps_genes <- c("invA", "avrA", "gtgA", "lpfB", "sspH2", "sseK2", "sseK3", "spvD")
hps_genes[!hps_genes %in% names(vir_matrix)]   # should be empty

hps <- vir_matrix %>% select(BioSample, all_of(hps_genes))
hps$n_hps_genes <- rowSums(hps[, hps_genes])
hps$HPS <- ifelse(hps$n_hps_genes >= 5, "HPS", "non-HPS")

table(hps$HPS)
write_csv(hps, "Output/hps.csv")
