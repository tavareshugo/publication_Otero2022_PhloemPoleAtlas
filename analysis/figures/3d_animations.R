library(scater)
library(tidyverse)
library(rgl)

# source util functions
source("analysis/functions/utils.R")

ring_hard <- readRDS("data/processed/SingleCellExperiment/ring_batches_strictfilt.rds")


# 3D UMAP -----------------------------------------------------------------

# UMAP with 3 dimensions
reducedDim(ring_hard, "temp") <-
  calculateUMAP(ring_hard,
                ncomponents = 3,
                dimred = "MNN_logvst",
                n_neighbors = 30)

# get the data to a tidy table
temp <- getReducedDim(ring_hard, "temp")

# add colours to a column
temp$col <- ggthemes::tableau_color_pal("Tableau 20")(length(unique(temp$cluster_mnn_logvst)))[as.numeric(as.factor(temp$cluster_mnn_logvst))]

rm(ring_hard); gc()

# Static chart
plot3d(temp$V1, temp$V2, -temp$V3,
       col = temp$col,
       type = "p", radius = .2, axes = FALSE, xlab = "", ylab = "", zlab = "")

# We can indicate the axis and the rotation velocity
play3d( spin3d( axis = c(0, 0, 1), rpm = 8), duration = 1/8*60 )

movie3d(
  spin3d(axis = c(0, 0, 1), rpm = 8),
  duration = 1/8*60, # 1/rpm * 60s
  fps = 50,
  movie="phloem_umap",
  dir = "~/Desktop",
  type = "gif",
  clean = TRUE
)



# Rayshader ---------------------------------------------------------------

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
pdata <- getReducedDim(ring_hard, "UMAP30_MNN_logvst",
                       genes = curated_genes$id, melted = TRUE,
                       exprs_values = "logvst") %>%
  left_join(curated_genes, by = "id") %>%
  arrange(!is.na(expr), expr)

p <- pdata %>%
  filter(name %in% c("SMB", "CALS7", "S17", "AHA3")) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  group_by(name) %>%
  mutate(expr_weighted = expr_weighted/max(expr_weighted, na.rm = TRUE)) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(ring_hard, "UMAP30_MNN_logvst"),
             colour = "lightgrey") +
  geom_point(aes(colour = expr_weighted)) +
  facet_wrap( ~ name) +
  scale_colour_viridis_c(na.value = NA) +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression") +
  theme(legend.position = "none")
p

library(rayshader)
plot_gg(p)
plot_gg(p, fov = 40)
