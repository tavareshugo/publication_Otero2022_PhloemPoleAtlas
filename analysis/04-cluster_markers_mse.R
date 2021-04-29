library(scran)
library(scater)
library(tidyverse)
library(patchwork)

theme_set(theme_minimal() + theme(text = element_text(size = 16)))
source("./analysis/functions/utils.R")

# Read data ---------------------------------------------------------------

ring <- readRDS("data/processed/SingleCellExperiment/ring_batches_strictfilt.rds")

plotReducedDim(ring, "UMAP30_MNN_logvst",
               colour_by = "cluster_mnn_logvst", text_by = "cluster_mnn_logvst")


# Find markers ------------------------------------------------------------

markers <- vroom::vroom("data/processed/gene_sets/ring_batches_strictfilt_cluster_markers_12in15.csv")
markers$cluster = factor(as.numeric(markers$cluster))

# add TF information
tfs <- read_csv("data/external/transcription_factors/PlantTFDB_tidy.csv")
markers$is_tf <- markers$id %in% tfs$gene

markers %>%
  filter(FDR < 0.01 & summary.AUC > 0.7) %>%
  distinct(cluster, id) %>%
  write_csv("~/temp/cluster_markers.csv")


# Identify candidate MSE markers - cluster 10 ------------------------------

# get markers for cluster 10
interesting_genes <- markers %>%
  filter(cluster == 10 & FDR < 0.01 & summary.AUC > 0.7) %>%
  pull(id)

# plot them all in a PDF document
pdf("~/temp/cluster10_markers.pdf")
for(i in interesting_genes){
  p1 <- ring %>%
    getReducedDim("UMAP30_MNN_logvst", genes = i, zeros_to_na = FALSE) %>%
    ggplot(aes(cluster_mnn_logvst, expr)) +
    geom_violin(scale = "width", fill = "lightgrey") +
    theme_classic() +
    labs(x = "Cluster", y = "logcounts")

  p2 <- ring %>%
    getReducedDim("UMAP30_MNN_logvst",
                  genes = i, melted = TRUE,
                  exprs_values = "logcounts") %>%
    arrange(!is.na(expr)) %>%
    group_by(cluster_mnn_logvst, id) %>%
    mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
    ggplot(aes(V1, V2)) +
    geom_point(aes(colour = expr), size = .5) +
    geom_label(stat = "centroid",
               aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
               alpha = 0.8, label.padding = unit(0.1, "lines")) +
    scale_colour_viridis_c(na.value = "lightgrey") +
    labs(x = "UMAP1", y = "UMAP2",
         colour = "Cluster-weighted\nNormalised\nExpression") +
    theme_void() + coord_equal() + theme(legend.position = "none")
  print((p1 | p2) + plot_annotation(title = i))
}
dev.off()

# sAPL - is it more expressed in cluster 10?
# not from this
markers %>%
  filter(id == "AT3G12730" & cluster == 10)

# but the percentage of cells where it is detected is reasonable
ring %>%
  getReducedDim("UMAP30_MNN_logvst", genes = "AT3G12730") %>%
  select(cluster_mnn_logvst, expr) %>%
  group_by(cluster_mnn_logvst) %>%
  summarise(pct = sum(!is.na(expr))/n()*100) %>%
  ggplot(aes(cluster_mnn_logvst, pct)) +
  geom_col() +
  labs(x = "Cluster", y = "% cells expressing AT3G12730")


# Identify candidate MSE markers - cluster 1 ------------------------------

# get markers for cluster 1
# relaxing threshold a bit to get a few more candidates
interesting_genes <- markers %>%
  filter(cluster == 1 & FDR < 0.01 & summary.AUC > 0.6) %>%
  pull(id)

# plot them all in a PDF document
pdf("~/temp/cluster1_markers.pdf")
for(i in interesting_genes){
  p1 <- ring %>%
    getReducedDim("UMAP30_MNN_logvst", genes = i, zeros_to_na = FALSE) %>%
    ggplot(aes(cluster_mnn_logvst, expr)) +
    geom_violin(scale = "width", fill = "lightgrey") +
    theme_classic() +
    labs(x = "Cluster", y = "logcounts")

  p2 <- ring %>%
    getReducedDim("UMAP30_MNN_logvst",
                  genes = i, melted = TRUE,
                  exprs_values = "logcounts") %>%
    arrange(!is.na(expr)) %>%
    group_by(cluster_mnn_logvst, id) %>%
    mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
    ggplot(aes(V1, V2)) +
    geom_point(aes(colour = expr), size = .5) +
    geom_label(stat = "centroid",
               aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
               alpha = 0.8, label.padding = unit(0.1, "lines")) +
    scale_colour_viridis_c(na.value = "lightgrey") +
    labs(x = "UMAP1", y = "UMAP2",
         colour = "Cluster-weighted\nNormalised\nExpression") +
    theme_void() + coord_equal() + theme(legend.position = "none")
  print((p1 | p2) + plot_annotation(title = i))
}
dev.off()



# Early clusters ----------------------------------------------------------

# get markers for cluster 8
interesting_genes <- markers %>%
  filter(cluster == 8 & FDR < 0.01 & summary.AUC > 0.7) %>%
  pull(id)

# plot them all in a PDF document
pdf("~/temp/cluster8_markers.pdf")
for(i in interesting_genes){
  p1 <- ring %>%
    getReducedDim("UMAP30_MNN_logvst", genes = i, zeros_to_na = FALSE) %>%
    ggplot(aes(cluster_mnn_logvst, expr)) +
    geom_violin(scale = "width", fill = "lightgrey") +
    theme_classic() +
    labs(x = "Cluster", y = "logcounts")

  p2 <- ring %>%
    getReducedDim("UMAP30_MNN_logvst",
                  genes = i, melted = TRUE,
                  exprs_values = "logcounts") %>%
    arrange(!is.na(expr)) %>%
    group_by(cluster_mnn_logvst, id) %>%
    mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
    ggplot(aes(V1, V2)) +
    geom_point(aes(colour = expr), size = .5) +
    geom_label(stat = "centroid",
               aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
               alpha = 0.8, label.padding = unit(0.1, "lines")) +
    scale_colour_viridis_c(na.value = "lightgrey") +
    labs(x = "UMAP1", y = "UMAP2",
         colour = "Cluster-weighted\nNormalised\nExpression") +
    theme_void() + coord_equal() + theme(legend.position = "none")
  print((p1 | p2) + plot_annotation(title = i))
}
dev.off()


# no markers for cluster 9
markers %>%
  filter(cluster == 9 & FDR < 0.01 & summary.AUC > 0.7) %>%
  pull(id)

# pct cells in each cluster per sample
ring %>%
  colData() %>%
  as_tibble() %>%
  count(Sample, cluster_mnn_logvst) %>%
  group_by(Sample) %>%
  mutate(pct = n/sum(n) * 100) %>%
  ungroup() %>%
  ggplot(aes(Sample, pct)) +
  geom_col(aes(fill = cluster_mnn_logvst)) +
  ggthemes::scale_fill_tableau("Tableau 20") +
  labs(x = "Sample", y = "% cells", fill = "Cluster")

