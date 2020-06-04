library(tidyverse)
library(furrr)
library(patchwork)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))


# Read data ---------------------------------------------------------------

# get wilcox test results for each cluster
cluster_test <- read_csv("data/processed/gene_sets/cluster_markers_hardfilt.csv",
                         guess_max = Inf) %>%
  mutate(cluster = factor(cluster))

# make list of marker genes
cluster_markers <- cluster_test %>%
  # retain genes with FDR < 5%
  filter(FDR < 0.05) %>%
  distinct(id, cluster)

# inferred cluster annotation - empirically determined based on markers and curated genes
cluster_annot <- tribble(
  ~cluster, ~type,
  "1", "PPP",
  "2", "Unknown",
  "3", "Unknown",
  "4", "PPP",
  "5", "SE",
  "6", "Cycling",
  "7", "Unknown",
  "8", "Unknown",
  "9", "Cycling",
  "10", "Unknown",
  "11", "Outer Layers",
  "12", "CC"
)


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

# Retain only TFs for which we have expression data
tfs <- tfs %>%
  filter(gene %in% cluster_test$id)

# calculate number and proportion of each family
tfs <- tfs %>%
  group_by(family) %>%
  mutate(n_family = n_distinct(gene)) %>%
  ungroup() %>%
  mutate(prop_family = n_family/n_distinct(gene))

# Note that some TFs are assigned to more than one family (because of isoforms)
tfs %>% group_by(gene) %>% filter(n() > 1) %>% arrange(gene)


# Note: the data from AGRIS does not seem as up-to-date, so we use the one from plantTFDB
# # source: AGRIS
# # http://www.plantphysiol.org/content/140/3/818.abstract
#
# # Transcription factors
# download.file("https://agris-knowledgebase.org/Downloads/AtTFDB.zip",
#               destfile = "temp/AtTFDB.zip")


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

# Retain only genes in our expression data
regnet <- regnet %>%
  filter(target %in% cluster_test$id)



# TF family enrichment ----------------------------------------------------

# Check how many genes are TFs and/or DE
cluster_test %>%
  distinct(id) %>%
  mutate(is_de = id %in% cluster_markers$id,
         is_tf = id %in% tfs$gene) %>%
  with(table(is_de, is_tf)) #%>% fisher.test()

# Check TF enrichment per comparison and family
genes_annotated <- cluster_test %>%
  # indicate whether it's differentially expressed
  mutate(is_de = FDR < 0.05) %>%
  # add TF information
  mutate(is_tf = id %in% tfs$gene) %>%
  # turn these into factors (for use with fisher.test)
  mutate(is_de = factor(is_de, levels = c("TRUE", "FALSE")),
         is_tf = factor(is_tf, levels = c("TRUE", "FALSE")))

genes_annotated %>%
  group_by(cluster) %>%
  nest() %>%
  mutate(fisher = map(data, ~ broom::tidy(fisher.test(.$is_de, .$is_tf)))) %>%
  unnest(fisher) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "bonferroni")) %>%
  ggplot(aes(factor(cluster), estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  geom_text(aes(label = paste0("p = ", round(padj, 4))), nudge_x = 0.3) +
  geom_hline(yintercept = 1, colour = "grey", linetype = 2) +
  scale_y_log10() +
  coord_flip() +
  labs(x = "Comparison", y = "Odds ratio") +
  theme_classic()

# Get counts
genes_annotated %>%
  count(cluster, is_tf, is_de) %>%
  group_by(cluster) %>%
  mutate(pct = n/sum(n))


# Look for enrichment of specific families
genes_annotated %>%
  right_join(tfs, by = c("id" = "gene")) %>%
  group_by(cluster, family, n_family, prop_family) %>%
  summarise(de_n = sum(is_de == "TRUE")) %>%
  # run binomial test
  group_by(cluster, family) %>%
  mutate(binom_test = list(binom.test(x = de_n, n = n_family, p = prop_family, alternative = "greater"))) %>%
  ungroup() %>%
  mutate(binom_test = map(binom_test, broom::tidy)) %>%
  unnest() %>%
  mutate(padj = p.adjust(p.value, method = "bonferroni")) %>%
  filter(padj < 0.05)


# TF target enrichment  ---------------------------------------------------

# Check how many targets are DE
cluster_test %>%
  distinct(id) %>%
  mutate(is_de = id %in% cluster_markers$id,
         is_tf = id %in% tfs$gene,
         is_target = id %in% regnet$target) %>%
  with(table(is_de, is_target))

# Get list of targets for each TF
targets_per_tf <- split(regnet, f = regnet$source) %>%
  map(~ pull(., target))

# Run test for each TF
target_enrich <- future_map_dfr(targets_per_tf, function(targets){
  genes_annotated %>%
    mutate(is_target = factor(id %in% targets, levels = c("TRUE", "FALSE"))) %>%
    nest(data = (-cluster)) %>%
    mutate(fisher = map(data, ~ broom::tidy(fisher.test(.$is_de, .$is_target, alternative = "greater")))) %>%
    unnest(fisher) %>%
    select(-data)
}, .id = "tf")

# Correct multiple testing
target_enrich <- target_enrich %>%
  mutate(padj = p.adjust(p.value, method = "fdr"))

# Save table of enriched TFs
target_enrich %>%
  filter(padj < 0.05) %>%
  group_by(tf) %>%
  summarise(clusters = paste(sort(cluster), collapse = ", "),
            n_clusters = n_distinct(cluster)) %>%
  write_csv("data/processed/gene_sets/tf_target_enrichment.csv")


# Make plot of effect size
target_enrich %>%
  filter(padj < 0.05) %>%
  mutate(tf = fct_reorder(tf, estimate)) %>%
  group_by(tf) %>%
  mutate(n_comparison = n_distinct(cluster)) %>%
  ggplot(aes(tf, estimate, colour = cluster)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, colour = "grey", linetype = 2) +
  facet_grid(n_comparison ~ ., scales = "free", space = "free") +
  scale_y_log10(breaks = c(1, 2, 5, 15)) +
  coord_flip() +
  labs(x = "Transcription factor", y = "Odds ratio") +
  theme_classic() +
  scale_colour_brewer(palette = "Dark2")


target_enrich %>%
  filter(padj < 0.05) %>%
  count(tf) %>%
  count(n)
