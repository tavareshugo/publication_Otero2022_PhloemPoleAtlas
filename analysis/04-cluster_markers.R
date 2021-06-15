library(scater)
library(tidyverse)
library(UpSetR)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))

# source custom functions
source("analysis/functions/utils.R")
source("analysis/functions/plots.R")


# Read data ---------------------------------------------------------------

# sce object
ring <- readRDS("data/processed/SingleCellExperiment/ring_batches_strictfilt.rds")
all <- readRDS("data/processed/SingleCellExperiment/all_batches_strictfilt.rds")

# get wilcox test results for each cluster
cluster_test <- vroom::vroom("data/processed/gene_sets/ring_batches_strictfilt_cluster_markers_8in15.csv") %>%
  mutate(cluster = factor(cluster))

# inferred cluster annotation - empirically determined based on markers and curated genes
cluster_annot <- tribble(
  ~cluster, ~type,
  "1", "MSE?",
  "2", "PSE",
  "3", "CC",
  "4", "PPP",
  "5", "CC",
  "6", "PSE",
  "7", "PPP",
  "8", "Unknown",
  "9", "Unknown",
  "10", "Unknown",
  "11", "PPP",
  "12", "PSE (early)",
  "13", "Cycling",
  "14", "PPP",
  "15", "Outer Layers"
)

# make transcription factor list for convenience
tfs <- cluster_test %>%
  filter(is_tf) %>%
  pull(id) %>% unique()


# EDA ---------------------------------------------------------------------

# Count how many transcription factors in each cluster
cluster_test %>%
  filter(FDR < 0.05) %>%
  group_by(cluster) %>%
  summarise(n_tf = sum(is_tf),
            prop_tf = sum(is_tf)/n()) %>%
  mutate(cluster = fct_reorder(cluster, prop_tf)) %>%
  ggplot(aes(cluster, prop_tf)) +
  geom_col() +
  coord_flip()

# upset plot
temp <- cluster_test %>%
  filter(FDR < 0.05 & summary.AUC > 0.7) %>%
  distinct(cluster, id) %>%
  with(split(id, cluster))
upset(fromList(temp), order.by = "freq", nsets = 12)


# CC clusters ----------------------------------------------------

# focus on cluster 3 (presumed CC)
interesting_genes <- cluster_test %>%
  # focus on cluster 3
  filter(cluster == 3) %>%
  # false discovery < 5% and average AUC of 70%
  filter(FDR < 0.05 & summary.AUC > 0.8) %>%
  # get gene ids
  arrange(desc(summary.AUC)) %>%
  pull(id) %>% unique()

# more fine-tuned filtering
interesting_genes <- cluster_test %>%
  # FDR threshold
  filter(FDR < 0.05) %>%
  # for each gene
  group_by(id) %>%
  # retain genes that are significant in only some clusters
  filter(n_distinct(cluster) == 2 &
           cluster %in% c(3, 5)) %>%
  ungroup() %>%
  # filter AUC threshold for better signal
  filter(summary.AUC > 0.7) %>%
  # or could filter to be transcription factor
  # filter(is_tf) %>%
  # arrange in descending order of AUC and get the gene IDs
  arrange(desc(summary.AUC)) %>%
  pull(id) %>% unique()

# violin plots of expression
plotExpression(ring,
               features = unique(c(interesting_genes[1:8], "AT1G22710", "AT5G02600")),
               x = "cluster_mnn_logvst", "logcounts")

# visualise with the Denyer dataset
all %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = interesting_genes, melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression")


# focus on cluster 3 but ensuring genes distinguish from clusters
# 14 & 4 (presumed PPP) and 6 (presumed SE)
interesting_genes <- cluster_test %>%
  filter(cluster == 3 & FDR < 0.05) %>%
  filter(AUC.14 > 0.7 & AUC.4 > 0.7 & AUC.6 > 0.7) %>%
  pull(id)

plotExpression(object = ring,
               features = interesting_genes,
               x = "cluster_mnn_logvst", "logcounts")


# Early clusters ----------------------------------------------------------

# use these to identify genes highly correlated with them
curated_early <- "AT1G29160"

# identify cyclin genes (from their names)
cyclins <- rownames(ring)[grepl("CYC[A,B,D]",
                                rowData(ring)$alternative_name,
                                ignore.case = TRUE)]

# cluster 13, 8, 9
interesting_genes <- cluster_test %>%
  # FDR threshold
  filter(FDR < 0.05) %>%
  # for each gene
  group_by(id) %>%
  # retain genes that are significant in only certain clusters
  filter(n_distinct(cluster) == 3 &
           all(cluster %in% c(13, 8, 9))) %>%
  ungroup() %>%
  # filter AUC threshold for better signal
  filter(summary.AUC > 0.8) %>%
  #filter(is_tf) %>%
  # arrange in descending order of AUC and get the gene IDs
  arrange(desc(summary.AUC)) %>%
  pull(id) %>% unique()

# none of these are cyclins
sum(interesting_genes %in% cyclins)

# plot
plotExpression(object = ring,
               features = c(interesting_genes, curated_early),
               x = "cluster_mnn_logvst", "logcounts")

# visualise with the Denyer dataset
all %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = "AT1G29160", melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression")


# PPP clusters -----------------------------------------------------------

##### !!! CLUSTERS ARE OUTDATED FROM HERE ON #####

# focus on clusters 1 & 4
interesting_genes <- cluster_test %>%
  # FDR threshold
  filter(FDR < 0.05) %>%
  # for each gene
  group_by(id) %>%
  # retain genes that are significant in only certain clusters
  filter(n_distinct(cluster) == 2 &
           all(cluster %in% c(1, 4))) %>%
  ungroup() %>%
  # filter AUC threshold for better signal
  filter(summary.AUC > 0.8) %>%
  # arrange in descending order of AUC and get the gene IDs
  arrange(desc(summary.AUC)) %>%
  pull(id) %>% unique()


plotExpression(object = ring,
               features = interesting_genes,
               x = "cluster_mnn_logvst", "logcounts")

# visualise with the Denyer dataset
all %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = interesting_genes, melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression")

# other genes that have been cloned before
plotExpression(object = ring,
               features = c("AT5G47450", "AT1G70990", "AT2G23560"),
               x = "cluster_mnn_logvst", "logcounts")


# Distinguish CC + PPP clusters from the rest ------------------------------

curated_ring <- c("AT3G16330", "AT3G21770", "AT2G02230", "AT1G52140")

# All clusters annotated as CC and PPP
interesting_genes <- cluster_test %>%
  # FDR threshold
  filter(FDR < 0.05) %>%
  # for each gene
  group_by(id) %>%
  # retain genes that are significant in only certain clusters
  filter(n_distinct(cluster) == 3 &
           all(cluster %in% c(1, 4, 12))) %>%
  ungroup() %>%
  # filter AUC threshold for better signal
  # filter(summary.AUC > 0.8) %>%
  filter(is_tf) %>%
  # arrange in descending order of AUC and get the gene IDs
  arrange(desc(summary.AUC)) %>%
  pull(id) %>% unique()


plotExpression(object = ring,
               features = c(interesting_genes[5:9], curated_ring),
               x = "cluster_mnn_logvst")

# visualise with the Denyer dataset
all %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = curated_ring, melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression")

ring %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = "AT3G01470", melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression")


ring_cor <- apply(assay(ring, "logvst")[curated_ring, ], 1, function(i){
  apply(assay(ring, "logvst"), 1, function(j){ cor(i, j)})
})



# Correlation approach ----------------------------------------------------

curated_ring <- c("AT3G16330", "AT3G21770", "AT2G02230", "AT1G52140")
curated_early <- "AT1G29160"

cormat <- cor(t(as.matrix(assay(ring, "logvst"))),
              t(as.matrix(assay(ring, "logvst")[c(curated_early, curated_ring), ])))

# write result
cormat %>%
  as_tibble(rownames = ".feature") %>%
  pivot_longer(starts_with("AT")) %>%
  group_by(name) %>%
  mutate(rank = rank(-value)) %>%
  arrange(rank) %>%
  ungroup() %>%
  select(id = .feature, reference_gene = name, correlation = value, rank) %>%
  write_csv("data/processed/gene_sets/correlation_with_reference_genes.csv")

# early ones
candidates <- cormat %>%
  as_tibble(rownames = ".feature") %>%
  pivot_longer(starts_with("AT")) %>%
  group_by(name) %>%
  mutate(rank = rank(-value)) %>%
  ungroup() %>%
  filter(name %in% curated_early) %>%
  filter(rank <= 8) %>%
  arrange(rank) %>%
  pull(.feature) %>%
  unique()

plotExpression(object = ring,
               features = c(candidates),
               x = "cluster_mnn_logvst")

all %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = c(candidates), melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  group_by(id) %>%
  mutate(expr_weighted = expr_weighted/max(expr_weighted, na.rm = TRUE)) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  # geom_label(stat = "centroid",
  #            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
  #            alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression")


# ring ones
candidates <- cormat %>%
  as_tibble(rownames = ".feature") %>%
  pivot_longer(starts_with("AT")) %>%
  group_by(name) %>%
  mutate(rank = rank(-value)) %>%
  ungroup() %>%
  filter(name %in% curated_ring[1] & !(.feature %in% curated_ring[-1])) %>%
  filter(rank <= 8) %>%
  arrange(rank) %>%
  pull(.feature) %>%
  unique()

plotExpression(object = ring,
               features = c(candidates),
               x = "cluster_mnn_logvst")

all %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = c(candidates), melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  group_by(id) %>%
  mutate(expr_weighted = expr_weighted/max(expr_weighted, na.rm = TRUE)) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  # geom_label(stat = "centroid",
  #            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
  #            alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression")


