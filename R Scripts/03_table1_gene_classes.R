## 03_table1_gene_classes.R
## Match AMR genes to the Reference Gene Catalog and make Table 1
## Reference Gene Catalog version 2026-08-07.1, downloaded 10/23/2026

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

setwd("~/USDA/Rule Mining Work")

refgenes <- read_tsv("Data/refgenes.tsv", col_types = cols(.default = "c"))
amr_gene_list <- read_csv("Output/amr_gene_list.csv")
amr_long <- read_csv("Output/amr_genes_long.csv")
isolates <- read_csv("Output/isolates_clean.csv", col_select = c(BioSample, specimen_type))

refgenes <- refgenes %>%
  rename(allele = `#Allele`, gene_family = `Gene family`)


# Match genes to the catalog
# look in the allele column first, then gene family
allele_lookup <- refgenes %>%
  filter(!is.na(allele)) %>%
  distinct(allele, .keep_all = TRUE) %>%
  select(gene = allele, gene_family, class = Class, subclass = Subclass)

family_lookup <- refgenes %>%
  arrange(!is.na(allele)) %>%
  distinct(gene_family, .keep_all = TRUE) %>%
  select(gene_family, family_class = Class, family_subclass = Subclass)

gene_map <- amr_gene_list %>%
  rename(n_isolates = n) %>%
  left_join(allele_lookup, by = "gene")

sum(is.na(gene_map$class))

no_match <- is.na(gene_map$class)
gene_map$gene_family[no_match] <- gene_map$gene[no_match]

# mutations not in the catalog (ramR_...) use the part before "_"
mutation <- no_match & !(gene_map$gene %in% family_lookup$gene_family)
gene_map$gene_family[mutation] <- str_remove(gene_map$gene[mutation], "_.*$")

gene_map <- gene_map %>% left_join(family_lookup, by = "gene_family")
gene_map$class[no_match] <- gene_map$family_class[no_match]
gene_map$subclass[no_match] <- gene_map$family_subclass[no_match]

gene_map <- gene_map %>% select(gene, gene_family, class, subclass, n_isolates)

gene_map$gene[is.na(gene_map$class)]


# Genes set by hand (broad names not in the catalog)
by_hand <- c("dfrA51"      = "TRIMETHOPRIM",
             "dfr"         = "TRIMETHOPRIM",
             "aac(3)-IV"   = "AMINOGLYCOSIDE",
             "aac(3)-II"   = "AMINOGLYCOSIDE",
             "aac(6')-Ie"  = "AMINOGLYCOSIDE",
             "aph(3'')"    = "AMINOGLYCOSIDE",
             "aph(2'')-I"  = "AMINOGLYCOSIDE",
             "aph(2'')-Ia" = "AMINOGLYCOSIDE",
             "ant(3'')"    = "AMINOGLYCOSIDE",
             "mph"         = "MACROLIDE",
             "erm"         = "MACROLIDE",
             "mle"         = "MACROLIDE",
             "vat"         = "STREPTOGRAMIN",
             "tet"         = "TETRACYCLINE")

fix <- is.na(gene_map$class) & gene_map$gene %in% names(by_hand)
gene_map$class[fix] <- by_hand[gene_map$gene[fix]]
gene_map$subclass[fix] <- gene_map$class[fix]
gene_map$gene_family[fix] <- gene_map$gene[fix]

# should be empty
gene_map %>% filter(is.na(class))

table(gene_map$class)
write_csv(gene_map, "Output/gene_class_map.csv")


# Table 1
table1 <- amr_long %>%
  inner_join(gene_map, by = "gene") %>%
  left_join(isolates, by = "BioSample") %>%
  group_by(class, subclass) %>%
  summarise(genes = paste(sort(unique(gene_family)), collapse = ", "),
            n_genes = n_distinct(gene),
            n_isolates = n_distinct(BioSample),
            in_dataset = paste(sort(unique(specimen_type)), collapse = "; ")) %>%
  ungroup() %>%
  arrange(class, desc(n_isolates))

table1
write_csv(table1, "Output/table1_gene_classes.csv")
