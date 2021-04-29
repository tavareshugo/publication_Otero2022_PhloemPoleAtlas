library(tidyverse)
library(patchwork)


# Download TF data --------------------------------------------------------

# source: PlantTFDB 4.0

# Transcription factor list
# https://academic.oup.com/nar/article/45/D1/D1040/2290936
if(!file.exists("data/external/transcription_factors/Ath_TF_list.gz")){
  download.file("http://planttfdb.cbi.pku.edu.cn/download/TF_list/Ath_TF_list.txt.gz",
                destfile = "data/external/transcription_factors/Ath_TF_list.gz")
}

# Read data
tfs <- read_tsv(gzfile("data/external/transcription_factors/Ath_TF_list.gz")) %>%
  rename_all(tolower)

# retain unique gene IDs (no isoforms)
tfs <- tfs %>%
  # retain info from major isoform
  distinct(gene_id, family) %>%
  rename(gene = gene_id) %>%
  # ensure gene IDs are all uppercase
  mutate(gene = toupper(gene))

# calculate number and proportion of each family
tfs <- tfs %>%
  group_by(family) %>%
  mutate(n_family = n_distinct(gene)) %>%
  ungroup() %>%
  mutate(prop_family = n_family/n_distinct(gene))

# Note that some TFs are assigned to more than one family (because of isoforms)
tfs %>% group_by(gene) %>% filter(n() > 1) %>% arrange(gene)

# write tidy data
write_csv(tfs, "data/external/transcription_factors/PlantTFDB_tidy.csv")



# Download regulatory network data ----------------------------------------

# Source: AGRIS

# Regulatory network
if(!file.exists("data/external/transcription_factors/AtRegNet.csv")){
  download.file("https://agris-knowledgebase.org/Downloads/AtRegNet.zip",
                destfile = "data/external/transcription_factors/AtRegNet.zip")
  unzip("data/external/transcription_factors/AtRegNet.zip",
        exdir = "data/external/transcription_factors/")
  unlink("data/external/transcription_factors/AtRegNet.zip")
}

# This file is mis-formatted, because they did not quote the reference column, which contains commas :(
# but the first few columns are imported OK, so we retain it as is
regnet <- read_csv("data/external/transcription_factors/AtRegNet.csv") %>%
  select(-Reference, -Note, -PubMedID)

# tidy up a bit
regnet <- regnet %>%
  # rename some variables
  rename(source = TFLocus, target = TargetLocus) %>%
  # ensure gene names are uppercase
  mutate(source = toupper(source), target = toupper(target))

write_csv(regnet, "data/external/transcription_factors/AGRIS_tidy.csv")



# Get genes annotated as regulators of transcription -----------------------

library(biomaRt)

# Specify the plants_mart
m <- useMart("plants_mart", host="plants.ensembl.org")

# Define dataset
m <- useMart("plants_mart", host="plants.ensembl.org", dataset="athaliana_eg_gene")

regulators <- getBM(attributes=c("ensembl_gene_id", "go_id"),
                    mart=m,
                    filters = "go",
                    values = c("GO:0006355")) %>%
  filter(go_id %in% c("GO:0006355")) %>%
  rename(gene = ensembl_gene_id)

write_csv(regulators, "data/external/transcription_factors/regulators_of_transcription.csv")
