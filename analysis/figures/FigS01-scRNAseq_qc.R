library(SingleCellExperiment)
library(tidyverse)
library(ggpointdensity)
library(ggridges)
library(patchwork)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")



# Read Data ---------------------------------------------------------------

ring <- readRDS("data/processed/SingleCellExperiment/ring_batches_strictfilt.rds")



# Summary stats -----------------------------------------------------------

# quantiles of number of detected genes
detected_genes <- round(quantile(colSums(counts(ring) > 0), c(0.5, 0.1, 0.9)))

# quantiles of total UMIs
total_umis <- round(quantile(colSums(counts(ring)), c(0.5, 0.1, 0.9)))

# put a title together from these stats
title <- glue::glue("{ncol(ring)} Cells with > 2000 genes
{nrow(ring)} Genes detected in > 100 cells")
subtitle <- glue::glue("Genes per cell: {detected_genes[1]} ({detected_genes[2]}-{detected_genes[3]})
UMIs per cell: {total_umis[1]} ({total_umis[2]}-{total_umis[3]})")


# Figure ------------------------------------------------------------------

p1 <- ggplot(tibble(x = 1:4, y = 1:4), aes(x, y)) +
  annotate("text", x = 1, y = 4,
           label = glue::glue("{ncol(ring)} Cells with > 2000 genes"),
           hjust = 0) +
  annotate("text", x = 1, y = 3,
           label = glue::glue("{nrow(ring)} Genes detected in > 100 cells"),
           hjust = 0) +
  annotate("text", x = 1, y = 2,
           label = glue::glue("{detected_genes[1]} ({detected_genes[2]}-{detected_genes[3]}) median (10%-90%) genes/cell"),
           hjust = 0) +
  annotate("text", x = 1, y = 1,
           label = glue::glue("{total_umis[1]} ({total_umis[2]}-{total_umis[3]}) median (10%-90%) UMIs/cell"),
           hjust = 0) +
  xlim(1, 5) + ylim(1, 5) +
  theme_void()

p2 <- ring %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = cluster_mnn_logvst), size = 0.5) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             size = 2,
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  ggthemes::scale_colour_tableau("Tableau 20") +
  labs(x = "UMAP1", y = "UMAP2") +
  guides(colour = "none") +
  coord_equal() +
  theme_void()

pdf("documents/pdf for figures/umap_with_clusters.pdf", width = 7.5, height = 3)
(p1 / plot_spacer() + plot_layout(heights = c(1, 2))) | p2
dev.off()

# Figure - Sup? ------------------------------------------------------------

p1 <- ring %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = detected)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  scale_colour_viridis_c(trans = "log10", limits = c(2e3, 14e3)) +
  labs(x = "UMAP1", y = "UMAP2", colour = "Detected\nGenes",
       subtitle = title) +
  coord_equal() +
  theme_void()

p2 <- ggplot(colData(ring),
       aes(detected, cluster_mnn_logvst,
           height = ..density..)) +
  geom_density_ridges(fill = "lightgrey", alpha = 0.5, stat = "density", trim = TRUE,
                      size = .3) +
  geom_vline(xintercept = 2000, linetype = 2) +
  scale_x_continuous(trans = "log10") +
  annotation_logticks(sides = "b") +
  labs(x = "No. detected genes", y = "Cluster",
       subtitle = subtitle) +
  theme_ridges()

p3 <- ggplot(colData(ring),
             aes(total, cluster_mnn_logvst,
                 height = ..density..)) +
  geom_density_ridges(fill = "lightgrey", alpha = 0.5, stat = "density", trim = TRUE,
                      size = .3) +
  geom_vline(xintercept = 3000, linetype = 2) +
  scale_x_continuous(trans = "log10") +
  annotation_logticks(sides = "b") +
  labs(x = "Total UMIs", y = "Cluster") +
  theme_ridges()

(
  ((plot_spacer() + p1 + plot_spacer()) + plot_layout(widths = c(1, 3, 1))) /
    (p2 | p3)
) + plot_annotation(tag_levels = "A")

