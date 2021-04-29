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
ring <- readRDS("data/processed/SingleCellExperiment/ring_batches_strictfilt.rds")

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
                                     "AT5G57350", "AHA3"),
                      `early` = tribble(~id, ~name,
                                        "AT1G29160", "DOF1.5",
                                        "AT2G36400", "GRF3",
                                        "AT5G52870", "MAKR5")) %>%
  bind_rows(.id = "tissue")

# plot data
pdata <- getReducedDim(ring, "UMAP30_MNN_logvst",
                       genes = curated_genes$id, melted = TRUE,
                       exprs_values = "logvst") %>%
  left_join(curated_genes, by = "id") %>%
  arrange(!is.na(expr))


# Figure ------------------------------------------------------------------

# expression - weighted by fraction of cells expressing gene in the cluster
# this weighting is to somehow combine the two visualisations below into one
pdata %>%
  filter(tissue == "early") %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ tissue + name) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2", colour = "Cluster-weighted\nNormalised\nExpression")

# boxplot of logcounts
getReducedDim(ring, "UMAP30_MNN_logvst",
              genes = curated_genes$id, melted = TRUE,
              exprs_values = "logcounts", zeros_to_na = FALSE) %>%
  left_join(curated_genes, by = "id") %>%
  filter(tissue == "early") %>%
  ggplot(aes(factor(cluster_mnn_logvst), expr)) +
  geom_violin(scale = "width") +
  facet_wrap(~ tissue + name)

# expression - as is
pdata %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ tissue + name) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2", colour = "Normalised\nExpression")


# density of points
pdata %>%
  ggplot(aes(V1, V2)) +
  geom_point(colour = "lightgrey") +
  geom_pointdensity(data = pdata %>% drop_na(expr)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ tissue + name) +
  scale_colour_viridis_c(option = "magma") +
  labs(x = "UMAP1", y = "UMAP2", colour = "Point\nDensity")
