library(SingleCellExperiment)
library(tidyverse)
library(ggpointdensity)
library(patchwork)

# set ggplot2 theme
theme_set(theme_classic() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")



# Read Data ---------------------------------------------------------------

# hard filtering
ring <- readRDS("data/processed/SingleCellExperiment/ring_batches_strictfilt.rds")


# Percent cells expressing genes ------------------------------------------

# genes of interest
gois <- c("AT1G29160", "AT5G59090", "AT5G47920", "AT1G79430", "AT3G12730")

pct_expressing <- ring %>%
  getReducedDim("UMAP30_MNN_logvst", genes = gois) %>%
  select(cluster_mnn_logvst, gene = id, expr) %>%
  group_by(cluster_mnn_logvst, gene) %>%
  summarise(pct = sum(!is.na(expr))/n()*100) %>%
  ungroup() %>%
  mutate(cluster_mnn_logvst = fct_inseq(cluster_mnn_logvst))

pct_expressing %>%
  ggplot(aes(cluster_mnn_logvst, pct)) +
  geom_col() +
  facet_wrap(~ gene) +
  theme_minimal()

pct_expressing %>%
  mutate(cluster_mnn_logvst = tidytext::reorder_within(cluster_mnn_logvst,
                                                       pct, gene)) %>%
  ggplot(aes(pct, cluster_mnn_logvst)) +
  geom_segment(aes(x = 0, xend = pct,
                   y = cluster_mnn_logvst, yend = cluster_mnn_logvst),
               colour = "grey") +
  geom_point() +
  facet_wrap(~ gene, ncol = 1, scales = "free") +
  theme_minimal() +
  tidytext::scale_y_reordered()

pct_expressing %>%
  ggplot(aes(cluster_mnn_logvst, pct)) +
  geom_line(aes(colour = gene, group = gene)) +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal() +
  ylim(0, 100)

pct_expressing %>%
  ggplot(aes(cluster_mnn_logvst, gene)) +
  geom_point(aes(size = pct, alpha = pct)) +
  # geom_text(aes(label = round(ifelse(pct > 10, pct, NA))),
  #           nudge_y = 0.25) +
  scale_size_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_alpha_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_minimal() +
  labs(x = "Cluster", y = "Gene",
       size = "% Cells",
       alpha = "% Cells")

pct_expressing %>%
  ggplot(aes(cluster_mnn_logvst, gene)) +
  geom_tile(aes(fill = pct)) +
  # geom_text(aes(label = round(ifelse(pct > 10, pct, NA))),
  #           nudge_y = 0.25) +
  scale_fill_gradient(low = "white", high = "black",
                      limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(x = "Cluster", y = "Gene", fill = "% cells")

pct_expressing %>%
  ggplot(aes(gene, pct)) +
  geom_col(aes(fill = gene)) +
  facet_wrap(~ cluster_mnn_logvst, nrow = 1) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "% cells gene detected", fill = "")

# faceted one is the one
pdf("documents/pdf for figures/pct_cells_ring_genes_detected.pdf",
    width = 7, height = 4)
pct_expressing %>%
  ggplot(aes(cluster_mnn_logvst, pct)) +
  geom_col() +
  facet_wrap( ~ gene, ncol = 1, strip.position = "right") +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  guides(fill = "none") +
  labs(x = "Cluster", y = "% cells gene detected", fill = "")
dev.off()
