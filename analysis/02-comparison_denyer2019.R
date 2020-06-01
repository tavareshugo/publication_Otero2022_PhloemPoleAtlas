library(data.table)
library(scater)
library(scran)
library(ggplot2)
library(patchwork)

# set seed for reproducible results
set.seed(1001)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")

# Marker genes whose promoters were used for cell sorting
markers <- data.table(name = c("APL", "MAKR5", "PEARdel", "S17", "sAPL"),
                      id = c("AT1G79430", "AT5G52870", "AT2G37590", "AT2G22850", "AT3G12730"))


# Read data ---------------------------------------------------------------

# ring data both soft and hard filtered
ring_soft <- readRDS("data/processed/SingleCellExperiment/ring_batches_softfilt.rds")
ring_hard <- readRDS("data/processed/SingleCellExperiment/ring_batches_hardfilt.rds")

# all batches, both soft and hard filtered
all_soft <- readRDS("data/processed/SingleCellExperiment/all_batches_softfilt.rds")
all_soft$dataset <- ifelse(grepl("denyer", all_soft$Sample), "Denyer et al 2019", "ring")

all_hard <- readRDS("data/processed/SingleCellExperiment/all_batches_hardfilt.rds")
all_hard$dataset <- ifelse(grepl("denyer", all_hard$Sample), "Denyer et al 2019", "ring")



# Comparison with Denyer - hard filter ------------------------------------

# plot things together
p1 <- ggplot(getReducedDim(all_hard, "UMAP-MNN_30"), aes(V1, V2, group = cluster_mnn)) +
  geom_point(aes(colour = cluster_mnn), show.legend = FALSE) +
  geom_label(stat = "centroid", aes(label = cluster_mnn), size = 3,
             label.padding = unit(0.1, "lines")) +
  ggthemes::scale_colour_tableau(palette = "Tableau 20") +
  labs(x = "UMAP 1", y = "UMAP 2")

p2 <- ggplot(getReducedDim(all_hard, "UMAP-MNN_30"), aes(V1, V2)) +
  geom_point(aes(colour = dataset)) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_colour_brewer(palette = "Dark2")

p3 <- ggplot(getReducedDim(all_hard, "UMAP-MNN_30"), aes(V1, V2)) +
  ggpointdensity::geom_pointdensity(show.legend = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_colour_viridis_c(option = "magma") +
  facet_grid(~ dataset)

(p1 | p2) / p3 + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")


# plot genes mentioned in Denyer2019 Fig 2D
temp <- getReducedDim(all_hard, "UMAP-MNN_30", genes = c("AT3G58190", "AT1G79430", "AT5G47450"),
              melted = TRUE)
ggplot(temp[order(expr, na.last = FALSE)], aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  scale_colour_viridis_c() +
  facet_grid(dataset ~ id) +
  labs(x = "UMAP 1", y = "UMAP 2")


# Comparison with Denyer - soft filter -------------------------------------

# plot things together
p1 <- ggplot(getReducedDim(all_soft, "UMAP-MNN_30"), aes(V1, V2, group = cluster_mnn)) +
  geom_point(aes(colour = cluster_mnn), show.legend = FALSE) +
  geom_label(stat = "centroid", aes(label = cluster_mnn), size = 3,
             label.padding = unit(0.1, "lines")) +
  ggthemes::scale_colour_tableau(palette = "Tableau 20") +
  labs(x = "UMAP 1", y = "UMAP 2")

p2 <- ggplot(getReducedDim(all_soft, "UMAP-MNN_30"), aes(V1, V2)) +
  geom_point(aes(colour = dataset)) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_colour_brewer(palette = "Dark2")

p3 <- ggplot(getReducedDim(all_soft, "UMAP-MNN_30"), aes(V1, V2)) +
  ggpointdensity::geom_pointdensity(show.legend = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_colour_viridis_c(option = "magma") +
  facet_grid(~ dataset)

(p1 | p2) / p3 + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")


# plot genes mentioned in Denyer2019 Fig 2D
temp <- getReducedDim(all_soft, "UMAP-MNN_30", genes = c("AT3G58190", "AT1G79430", "AT5G47450"),
                      melted = TRUE)
temp$expr <- ifelse(temp$expr == 0, NA, temp$expr)
ggplot(temp[order(expr, na.last = FALSE)], aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  scale_colour_viridis_c() +
  facet_grid(dataset ~ id) +
  labs(x = "UMAP 1", y = "UMAP 2")


# Contaminating layers ----------------------------------------------------
# plot genes identifying collumela and epidermis

# hard filtered data
temp <- getReducedDim(all_hard, "UMAP-MNN_30",
                      genes = c("AT1G79580", "AT3G29810", "AT5G58580", "AT5G14750"),
                      melted = TRUE)
temp$expr <- ifelse(temp$expr == 0, NA, temp$expr)
ggplot(temp[order(expr, na.last = FALSE)], aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  scale_colour_viridis_c() +
  facet_grid(dataset ~ id) +
  labs(x = "UMAP 1", y = "UMAP 2")


# soft filtered data
temp <- getReducedDim(all_soft, "UMAP-MNN_30",
                      genes = c("AT1G79580", "AT3G29810", "AT5G58580", "AT5G14750"),
                      melted = TRUE)
temp$expr <- ifelse(temp$expr == 0, NA, temp$expr)
ggplot(temp[order(expr, na.last = FALSE)], aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  scale_colour_viridis_c() +
  facet_grid(dataset ~ id) +
  labs(x = "UMAP 1", y = "UMAP 2")


# Collumela
temp <- getReducedDim(all_hard, "UMAP-MNN_30", genes = c("AT1G79580", "AT1G13620", "AT2G04025"),
                      melted = TRUE)
ggplot(temp[order(expr, na.last = FALSE)], aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  scale_colour_viridis_c() +
  facet_grid(dataset ~ id) +
  labs(x = "UMAP 1", y = "UMAP 2")

# Epidermis - actually AT2G46410 is also expressed in "stele"
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC150583/
temp <- getReducedDim(all_hard, "UMAP-MNN_30", genes = c("AT5G14750", "AT1G79840", "AT2G46410"),
                      melted = TRUE)
ggplot(temp[order(expr, na.last = FALSE)], aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  scale_colour_viridis_c() +
  facet_grid(dataset ~ id) +
  labs(x = "UMAP 1", y = "UMAP 2")

# Genes used for sorting
temp <- getReducedDim(all_hard, "UMAP-MNN_30",
                      genes = unique(markers$id), melted = TRUE)
temp <- merge(temp, markers, by = "id")
ggplot(temp[order(expr, na.last = FALSE)], aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  scale_colour_viridis_c() +
  facet_grid(dataset ~ name) +
  labs(x = "UMAP 1", y = "UMAP 2")


