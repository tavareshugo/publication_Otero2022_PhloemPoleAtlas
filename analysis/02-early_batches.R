library(SingleCellExperiment)
library(tidyverse)
library(ggpointdensity)
library(patchwork)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")



# Read Data ---------------------------------------------------------------

# hard filtering
ring_hard <- readRDS("data/processed/SingleCellExperiment/ring_batches_hardfilt.rds")
early_hard <- readRDS("data/processed/SingleCellExperiment/early_batches_hardfilt.rds")

# curated list of genes
curated_genes <- list(`Outer Layers` = tribble(~id, ~name,
                                               "AT1G79580", "SMB",
                                               "AT1G79840", "GL2",
                                               "AT5G14750", "WER",
                                               "AT3G54220", "SCR"),
                      `SE (early)` = tribble(~id, ~name,
                                             "AT1G05470", "CVP2",
                                             "AT1G54330", "NAC020",
                                             "AT2G37590", "PEAR1",
                                             "AT5G02460", "PEAR2"),
                      `SE (late)` = tribble(~id, ~name,
                                            "AT1G06490", "CALS7",
                                            "AT5G17260", "NAC086"),
                      `PPP` = tribble(~id, ~name,
                                      "AT2G22850", "S17",
                                      "AT3G14570", "CALS8"),
                      `CC` = tribble(~id, ~name,
                                     "AT2G38640", "AT2G38640",
                                     "AT3G12730", "SAPL",
                                     "AT1G22710", "SUC2",
                                     "AT5G02600", "NAKR1",
                                     "AT5G57350", "AHA3"),
                      `early` = tribble(~id, ~name,
                                        "AT1G29160", "DOF1.5",
                                        "AT2G36400", "GRF3",
                                        "AT5G52870", "MAKR5")) %>%
  bind_rows(.id = "tissue")



# cyclins
getReducedDim(early_hard, "UMAP30_MNN_logvst",
              genes = cyclins$ID, melted = TRUE,
              exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(cyclins, by = c("id" = "ID")) %>%
  group_by(id) %>%
  mutate(expr_center = expr - mean(expr)) %>%
  group_by(cluster_mnn_logvst, id, alternative_name, cyclin) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  ggplot(aes(factor(cluster_mnn_logvst), alternative_name)) +
  geom_point(aes(colour = mean_expr, size = pct_expressing)) +
  facet_grid(cyclin ~ ., scales = "free", space = "free") +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "#f7f7f7") +
  scale_size_continuous(limits = c(0, 100)) +
  labs(x = "Cluster", y = "Cyclin", size = "% Expressing", colour = "Mean-centered\nExpression")

early_hard %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = cyclins$ID, melted = TRUE,
                exprs_values = "logvst") %>%
  left_join(cyclins, by = c("id" = "ID")) %>%
  arrange(!is.na(expr), expr) %>%
  filter(alternative_name == "CYCB1-1") %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ alternative_name) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2", colour = "Cluster-weighted\nNormalised\nExpression",
       title = "Early batches (MAKR5, MAKR5diff, PEARdel)")


early_hard %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = c("AT2G34140"), melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr), expr) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(. ~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2", colour = "Cluster-weighted\nNormalised\nExpression",
       title = "Early batches (MAKR5, MAKR5diff, PEARdel)")

early_hard %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = c("AT5G65470"), melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr), expr) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = Sample)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ Sample)


# Painted UMAPs -------------------------------------------------------------

tissue_of_interest <- "early"

p1 <- ring_hard %>%
  getReducedDim("UMAP30_MNN_logvst",
              genes = curated_genes$id, melted = TRUE,
              exprs_values = "logvst") %>%
  left_join(curated_genes, by = "id") %>%
  arrange(!is.na(expr), expr) %>%
  filter(tissue == tissue_of_interest) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ tissue + name) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2", colour = "Cluster-weighted\nNormalised\nExpression",
       title = "All batches")

p2 <- early_hard %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = curated_genes$id, melted = TRUE,
                exprs_values = "logvst") %>%
  left_join(curated_genes, by = "id") %>%
  arrange(!is.na(expr), expr) %>%
  filter(tissue == tissue_of_interest) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ tissue + name) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2", colour = "Cluster-weighted\nNormalised\nExpression",
       title = "Early batches (MAKR5, MAKR5diff, PEARdel)")

p1 / p2


# boxplot of log counts
p1.1 <- getReducedDim(ring_hard, "UMAP30_MNN_logvst",
              genes = curated_genes$id, melted = TRUE,
              exprs_values = "logcounts", zeros_to_na = FALSE) %>%
  left_join(curated_genes, by = "id") %>%
  filter(tissue == tissue_of_interest) %>%
  ggplot(aes(factor(cluster_mnn_logvst), expr)) +
  geom_violin(scale = "width") +
  facet_wrap(~ tissue + name)

p2.1 <- getReducedDim(early_hard, "UMAP30_MNN_logvst",
              genes = curated_genes$id, melted = TRUE,
              exprs_values = "logcounts", zeros_to_na = FALSE) %>%
  left_join(curated_genes, by = "id") %>%
  filter(tissue == tissue_of_interest) %>%
  ggplot(aes(factor(cluster_mnn_logvst), expr)) +
  geom_violin(scale = "width") +
  facet_wrap(~ tissue + name)

p1 / p1.1
p2 / p2.1



# Cell types --------------------------------------------------------------

# UMAP with coloured clusters
p1 <- ring_hard %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = cluster_mnn_logvst)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  ggthemes::scale_colour_tableau("Tableau 20") +
  labs(x = "UMAP1", y = "UMAP2",
       title = "All batches")

p1.1 <- getReducedDim(ring_hard, "UMAP30_MNN_logvst",
                      genes = curated_genes$id, melted = TRUE,
                      exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(curated_genes, by = "id") %>%
  group_by(id) %>%
  mutate(expr_center = (expr - mean(expr))/sd(expr)) %>%
  group_by(cluster_mnn_logvst, id, name, tissue) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  ggplot(aes(factor(cluster_mnn_logvst), name)) +
  geom_point(aes(colour = mean_expr, size = pct_expressing)) +
  facet_grid(tissue ~ ., scales = "free", space = "free") +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "#f7f7f7") +
  scale_size_continuous(limits = c(0, 100)) +
  guides(shape = "none") +
  labs(x = "Cluster", y = "", size = "% Expressing",
       colour = "Z-score\nExpression")

p2 <- early_hard %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = cluster_mnn_logvst)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  ggthemes::scale_colour_tableau("Tableau 20") +
  labs(x = "UMAP1", y = "UMAP2",
       title = "Early batches (MAKR5, MAKR5diff, PEARdel)")

p2.1 <- getReducedDim(early_hard, "UMAP30_MNN_logvst",
                      genes = curated_genes$id, melted = TRUE,
                      exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(curated_genes, by = "id") %>%
  group_by(id) %>%
  mutate(expr_center = (expr - mean(expr))/sd(expr)) %>%
  group_by(cluster_mnn_logvst, id, name, tissue) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  ggplot(aes(factor(cluster_mnn_logvst), name)) +
  geom_point(aes(colour = mean_expr, size = pct_expressing)) +
  facet_grid(tissue ~ ., scales = "free", space = "free") +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "#f7f7f7") +
  scale_size_continuous(limits = c(0, 100)) +
  guides(shape = "none") +
  labs(x = "Cluster", y = "", size = "% Expressing",
       colour = "Z-score\nExpression")

p1 / p1.1
p2 / p2.1


# cross-reference clusters
clusters <- colData(early_hard)[, c("cluster_mnn_logvst"), drop = FALSE]
colnames(clusters) <- c("early_clusters")
clusters$ring_clusters <- ring_hard[, rownames(clusters)]$cluster_mnn_logvst

clusters %>%
  ggplot(aes(early_clusters, ring_clusters)) +
  geom_count()

