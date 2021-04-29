library(SingleCellExperiment)
library(tidyverse)
library(readxl)
library(patchwork)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))

source("analysis/functions/utils.R")

# Read SCE data -----------------------------------------------------------

# slingshot
sling <- readRDS("data/processed/trajectories/ring_batches_strictfilt_slingshot.rds")



# APL and NEN85 -----------------------------------------------------------

genes_of_interest <- tribble(
  ~name, ~gene,
  "APL", "AT1G79430",
  "NAC86", "AT5G17260",
  "NEN4", "AT4G39810")

p1 <- sling %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = genes_of_interest$gene) %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  select(cell, cluster_mnn_logvst, matches("slingPseudotime"), id, expr) %>%
  pivot_longer(cols = matches("slingPseudotime"),
               names_to = "traj", values_to = "pseudotime") %>%
  drop_na() %>%
  mutate(traj = str_replace(traj, "slingPseudotime_", "Traj")) %>%
  filter(traj %in% paste0("Traj", 1:3)) %>%
  # group_by(traj) %>%
  # mutate(pseudotime = min_max(pseudotime, TRUE)) %>%
  # ungroup() %>%
  ggplot(aes(pseudotime, expr)) +
  geom_point(aes(colour = cluster_mnn_logvst), size = 0.5) +
  geom_smooth(se = FALSE, colour = "black",
              method = mgcv::gam, formula = y ~ s(x, bs = "cr", k = 7)) +
  facet_grid(traj ~ id) +
  ggthemes::scale_colour_tableau("Tableau 20") +
  guides(colour = "none") +
  labs(colour = "Cluster", y = "logcounts") +
  theme_minimal()

p2 <- sling %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = genes_of_interest$gene) %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  select(cell, cluster_mnn_logvst, matches("slingPseudotime"), id, expr) %>%
  pivot_longer(cols = matches("slingPseudotime"),
               names_to = "traj", values_to = "pseudotime") %>%
  drop_na() %>%
  mutate(traj = str_replace(traj, "slingPseudotime_", "Traj")) %>%
  filter(traj %in% paste0("Traj", 1:3)) %>%
  ggplot(aes(pseudotime, expr)) +
  geom_smooth(aes(colour = traj),
              se = FALSE,
              method = mgcv::gam, formula = y ~ s(x, bs = "cr", k = 7)) +
  facet_grid( ~ id, scales = "free", space = "free") +
  scale_colour_brewer(palette = "Dark2") +
  labs(colour = "Trajectory", y = "logcounts") +
  theme_minimal()

pdf("documents/pdf for figures/trajectories_apl.pdf",
    width = 4.68, height = 3)
p1 + theme(panel.grid.minor = element_blank())
dev.off()

pdf("documents/pdf for figures/trajectories_apl_compact.pdf",
    width = 4.68, height = 1.5)
p2 + theme(panel.grid.minor = element_blank())
dev.off()



# Soft data ---------------------------------------------------------------

ring <- readRDS("data/processed/SingleCellExperiment/ring_batches_hardfilt.rds")

temp <- logcounts(ring[genes_of_interest$gene, ]) %>%
  as.matrix() %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "expr")

p1 <- sling %>%
  colData() %>%
  as_tibble(rownames = "cell") %>%
  select(cell, cluster_mnn_logvst, matches("slingPseudo")) %>%
  pivot_longer(matches("sling"), names_to = "traj", values_to = "pseudotime") %>%
  drop_na() %>%
  mutate(traj = str_replace(traj, "slingPseudotime_", "Traj")) %>%
  filter(traj %in% paste0("Traj", 1:3)) %>%
  inner_join(temp, by = "cell") %>%
  inner_join(genes_of_interest, by = "gene") %>%
  # scale pseudotime between 0-1
  group_by(traj, name) %>%
  mutate(pseudotime = (pseudotime - min(pseudotime))/(max(pseudotime) - min(pseudotime))) %>%
  ggplot(aes(pseudotime, expr)) +
  geom_point(aes(colour = cluster_mnn_logvst), size = 0.5) +
  geom_smooth(se = FALSE, colour = "black",
              method = mgcv::gam, formula = y ~ s(x, bs = "cr", k = 7)) +
  facet_grid(traj ~ name) +
  ggthemes::scale_colour_tableau("Tableau 20") +
  guides(colour = "none") +
  labs(colour = "Cluster", y = "logcounts", x = "Pseudotime (scaled)") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

p2 <- sling %>%
  colData() %>%
  as_tibble(rownames = "cell") %>%
  select(cell, cluster_mnn_logvst, matches("slingPseudo")) %>%
  pivot_longer(matches("sling"), names_to = "traj", values_to = "pseudotime") %>%
  drop_na() %>%
  mutate(traj = str_replace(traj, "slingPseudotime_", "Traj")) %>%
  filter(traj %in% paste0("Traj", 1:3)) %>%
  inner_join(temp, by = "cell") %>%
  inner_join(genes_of_interest, by = "gene") %>%
  # scale pseudotime between 0-1
  group_by(traj, name) %>%
  mutate(pseudotime = (pseudotime - min(pseudotime))/(max(pseudotime) - min(pseudotime))) %>%
  ggplot(aes(pseudotime, expr)) +
  geom_smooth(aes(colour = traj),
              se = FALSE,
              method = mgcv::gam, formula = y ~ s(x, bs = "cr", k = 7)) +
  facet_grid( ~ name, scales = "free", space = "free") +
  scale_colour_brewer(palette = "Dark2") +
  labs(colour = "Trajectory", y = "logcounts", x = "Pseudotime (scaled)") +
  theme_minimal() +
  theme(axis.text.x = element_blank())


pdf("documents/pdf for figures/trajectories_apl.pdf",
    width = 4.68, height = 3)
p1 + theme(panel.grid.minor.y = element_blank(),
           panel.grid.minor.x = element_blank())
dev.off()

pdf("documents/pdf for figures/trajectories_apl_compact.pdf",
    width = 4.68, height = 1.5)
p2 + theme(panel.grid.minor.y = element_blank(),
           panel.grid.minor.x = element_blank())
dev.off()

