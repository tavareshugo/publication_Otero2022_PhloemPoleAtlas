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


# Figure ------------------------------------------------------------------

# make prettier sample labels
sample_labels <- colData(ring) %>%
  as_tibble() %>%
  count(Sample) %>%
  mutate(sample_label = case_when(Sample == "MAKR5diff" ~ "MAKR5 differentiated",
                            Sample == "PEARdel" ~ "PEAR1(delta)",
                            Sample == "sAPL" ~ "SAPL",
                            TRUE ~ Sample)) %>%
  mutate(sample_label = paste0(sample_label, "\nn = ", n, "")) %>%
  select(Sample, sample_label)

colData(ring) <- merge(colData(ring), sample_labels)


# plot
p1 <- ring %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(ring, "UMAP30_MNN_logvst") %>% select(-sample_label),
             colour = "lightgrey", size = 0.2) +
  geom_pointdensity(show.legend = FALSE, size = 0.2) +
  # geom_label(data = getReducedDim(ring, "UMAP30_MNN_logvst") %>% select(-sample_label),
  #            stat = "centroid",
  #            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
  #            alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ sample_label) +
  scale_colour_viridis_c(option = "magma") +
  labs(x = "UMAP1", y = "UMAP2", colour = "Point\nDensity") +
  theme_void() + coord_equal()

pdf("documents/pdf for figures/umap_by_sample.pdf", width = 7.5, height = 4)
p1
dev.off()
