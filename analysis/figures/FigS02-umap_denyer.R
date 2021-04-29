library(SingleCellExperiment)
library(tidyverse)
library(ggpointdensity)
library(ggridges)
library(patchwork)

# set ggplot2 theme
theme_set(theme_classic() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")



# Read Data ---------------------------------------------------------------

# hard filtering
all_hard <- readRDS("data/processed/SingleCellExperiment/all_batches_hardfilt.rds")

# curated list of genes
curated_genes <- list(`Outer Layers` = tribble(~id, ~name,
                                               "AT1G79580", "SMB",
                                               "AT1G79840", "GL2",
                                               "AT5G14750", "WER"),
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
                                     "AT5G57350", "AHA3")) %>%
  bind_rows(.id = "tissue")

# plot data
pdata <- getReducedDim(all_hard, "UMAP30_MNN_logvst",
                       genes = curated_genes$id, melted = TRUE,
                       exprs_values = "logvst") %>%
  left_join(curated_genes, by = "id") %>%
  arrange(!is.na(expr), expr)


# Figure ------------------------------------------------------------------

# UMAP with clusters
p1 <- ggplot(getReducedDim(all_hard, "UMAP30_MNN_logvst"),
       aes(V1, V2)) +
  geom_point(aes(colour = cluster_mnn_logvst),
             show.legend = FALSE) +
  geom_point(stat = "centroid", aes(group = cluster_mnn_logvst),
             shape = 21, size = 5, fill = "white", alpha = 0.7) +
  geom_text(stat = "centroid",
            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst)) +
  ggthemes::scale_colour_tableau(palette = "Tableau 20") +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Cluster")

# Data in two batches
p2 <- all_hard %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(dataset = ifelse(grepl("denyer", Sample), "Denyer et al", "Ring-enriched")) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(all_hard, "UMAP30_MNN_logvst"),
             colour = "lightgrey") +
  geom_pointdensity(show.legend = FALSE) +
  facet_grid(dataset ~ .) +
  scale_colour_viridis_c(option = "magma") +
  labs(x = "UMAP1", y = "UMAP2", colour = "Point\nDensity") +
  theme(axis.title = element_blank())

# Expression of outer layer genes in each batch
p3 <- pdata %>%
  mutate(dataset = ifelse(grepl("denyer", Sample), "Denyer et al", "Ring-enriched")) %>%
  filter(tissue == "Outer Layers") %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(all_hard, "UMAP30_MNN_logvst"),
             colour = "lightgrey") +
  geom_point(aes(colour = expr_weighted)) +
  # geom_label(stat = "centroid",
  #            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
  #            alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_grid(dataset ~ name) +
  scale_colour_viridis_c(na.value = NA) +
  labs(x = "UMAP1", y = "UMAP2", colour = "Cluster-weighted\nNormalised\nExpression")

(((p1 | p2) + plot_layout(widths = c(2, 1))) / p3) +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(axis.text = element_blank(), axis.ticks = element_blank())


# Another possibility

# Data in two batches
p1 <- all_hard %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(dataset = ifelse(grepl("denyer", Sample), "Denyer et al", "Ring-enriched")) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(all_hard, "UMAP30_MNN_logvst"),
             colour = "lightgrey") +
  geom_pointdensity(show.legend = FALSE) +
  facet_grid( ~ dataset) +
  scale_colour_viridis_c(option = "magma") +
  labs(x = "UMAP1", y = "UMAP2", colour = "Point\nDensity")

# Expression of outer layer genes in each batch
p2 <- pdata %>%
  mutate(dataset = ifelse(grepl("denyer", Sample), "Denyer et al", "Ring-enriched")) %>%
  filter(tissue == "Outer Layers") %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(all_hard, "UMAP30_MNN_logvst"),
             colour = "lightgrey") +
  geom_point(aes(colour = expr_weighted)) +
  facet_grid(name ~ dataset) +
  scale_colour_viridis_c(na.value = NA) +
  labs(x = "UMAP1", y = "UMAP2", colour = "Cluster-weighted\nNormalised\nExpression")

p1 / p2 +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 2)) &
  theme(axis.text = element_blank(), axis.ticks = element_blank())
