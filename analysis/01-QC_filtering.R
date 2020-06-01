library(data.table)
library(scater)
library(scran)
library(ggplot2)
library(ggpointdensity)
library(ggridges)
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


# read data ---------------------------------------------------------------

# soft filtering (only on mitochondrial genes)
ring_soft <- readRDS("data/processed/SingleCellExperiment/ring_batches_softfilt.rds")

# hard filtering (also total UMIs > 1000)
ring_hard <- readRDS("data/processed/SingleCellExperiment/ring_batches_hardfilt.rds")


# UMI statistics and filters ----------------------------------------------

# UMI count distribution
p1 <- ggplot(colData(ring_soft), aes(total, Sample)) +
  geom_density_ridges(fill = "lightgrey", alpha = 0.5) +
  geom_vline(xintercept = 1000, linetype = "dashed") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  labs(x = "Total UMIs per cell", y = "")

# Detected genes distribution
p2 <- ggplot(colData(ring_soft), aes(detected, Sample)) +
  geom_density_ridges(fill = "lightgrey", alpha = 0.5) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  labs(x = "No. detected genes", y = "")

p3 <- ggplot(colData(ring_soft), aes(total, detected)) +
  geom_pointdensity() +
  geom_vline(xintercept = 1000, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks(sides = "bl") +
  scale_colour_viridis_c(guide = "none", option = "magma") +
  facet_wrap(~ Sample)

(p1 + p2) / p3 + plot_annotation(tag_levels = "A")

# Checking how many cells pass different thresholds
coldat <- as.data.table(colData(ring_soft))
coldat[, .(no_filt = .N,
           mito_filt = sum(subsets_mitochondria_percent < 10),
           total_filt = sum(subsets_mitochondria_percent < 10 & total > 1e3)),
       by = Sample]


# Total UMIs - to filter or not to filter? -------------------------------

# Correlation of PC scores and total UMI
pc_cor <- data.frame(PC = factor(1:12),
                     cor = sapply(1:12, function(i) cor(ring_soft$total,
                                                        reducedDim(ring_soft, "PCA")[, i],
                                                        method = "spearman")),
                     var = attr(reducedDim(ring_soft, "PCA"), "percentVar")[1:12])

p1 <- ggplot(pc_cor, aes(PC, cor)) +
  geom_label(aes(label = round(var))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "PC", y = "Spearman Correlation") +
  theme_classic()

# Visualise this correlation
p2 <- ggplot(getReducedDim(ring_soft, "PCA"), aes(total, PC2)) +
  geom_point() +
  scale_x_log10() +
  labs(x = "Total UMIs",
       y = paste0("PC2 (", round(attr(reducedDim(ring_soft, "PCA"), "percentVar")[2]), "%)"),
       subtitle = paste0("Spearman correlation = ",
                         round(
                           cor(ring_soft$total,
                               getReducedDim(ring_soft, "PCA")$PC2, method = "spearman"),
                           2)))

p3 <- ggplot(getReducedDim(ring_soft, "PCA"), aes(PC1, PC2, colour = total)) +
  geom_point() +
  scale_colour_viridis_c(trans = "log10") +
  labs(x = paste0("PC1 (", round(attr(reducedDim(ring_soft, "PCA"), "percentVar")[1]), "%)"),
       y = paste0("PC2 (", round(attr(reducedDim(ring_soft, "PCA"), "percentVar")[2]), "%)")) +
  coord_fixed(ratio = 0.8)

p1 + p2 + p3 + plot_layout(ncol = 2)


# clustering without PC2
reducedDim(ring_soft, "temp") <- reducedDim(ring_soft, "PCA")[, -2]
graph_pca <- buildSNNGraph(ring_soft, k = 100, use.dimred = "temp",
                           type = "jaccard")
ring_soft$cluster_pca_minus2 <- factor(igraph::cluster_louvain(graph_pca)$membership)
rm(graph_pca)

# UMAP without using PC2
reducedDim(ring_soft, "UMAPtemp") <- calculateUMAP(ring_soft,
                                               dimred = "temp",
                                               n_neighbors = 30)

p1 <- ggplot(getReducedDim(ring_soft, "UMAPall_30"), aes(V1, V2)) +
  geom_point(aes(colour = total)) +
  geom_label(stat = "centroid",
             aes(group = cluster_pca, label = cluster_pca),
             alpha = 0.5) +
  scale_colour_viridis_c(trans = "log10") +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "UMAP including PC2") +
  coord_fixed()

p2 <- ggplot(getReducedDim(ring_soft, "UMAPtemp"), aes(V1, V2)) +
  geom_point(aes(colour = total)) +
  geom_label(stat = "centroid",
             aes(group = cluster_pca_minus2, label = cluster_pca_minus2),
             alpha = 0.5) +
  scale_colour_viridis_c(trans = "log10") +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "UMAP excluding PC2") +
  coord_fixed()

p3 <- ggplot(colData(ring_soft), aes(cluster_pca, cluster_pca_minus2)) +
  geom_count() +
  coord_fixed() +
  labs(x = "Cluster (all PCs)", y = "Cluster (excluding PC2)")

(p1 | p2) / (p3 | plot_spacer()) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")


# Batch effects - hard filter ---------------------------------------------

# UMAP on original and corrected data
p1 <- ggplot(getReducedDim(ring_hard, type = "UMAPall_30"),
             aes(V1, V2)) +
  geom_point(aes(colour = Sample), alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "No batch correction",
       caption = "UMAP ran on all PCs")

pcvar <- data.table(pc = 1:50,
                    var = attr(reducedDim(ring_hard, "PCA"), "percentVar"))
p2 <- ggplot(pcvar, aes(pc, var)) +
  geom_col() +
  geom_point(aes(y = cumsum(var))) +
  geom_line(aes(y = cumsum(var))) +
  geom_vline(xintercept = 10, colour = "dodgerblue", linetype = 2) +
  labs(x = "PC", y = "% Variance") +
  scale_y_continuous(limits = c(0, 80))

p3 <- ggplot(getReducedDim(ring_hard, type = "UMAP_30"),
       aes(V1, V2)) +
  geom_point(aes(colour = Sample), alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "No batch correction",
       caption = "UMAP ran on first 10 PCs")

p4 <- ggplot(getReducedDim(ring_hard, type = "UMAP-MNN_30"),
       aes(V1, V2)) +
  geom_point(aes(colour = Sample), alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "Batch correction",
       caption = "UMAP ran on MNN-corrected data")

(p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")


# comparing clustering using 50 PCs, 10 PCs and MNN-corrected data
graph_50pcs <- buildSNNGraph(ring_hard, k = 20, d = 50, type = "jaccard")
ring_hard$cluster_50pcs <- factor(igraph::cluster_louvain(graph_50pcs)$membership)

graph_10pcs <- buildSNNGraph(ring_hard, k = 20, d = 10, type = "jaccard")
ring_hard$cluster_10pcs <- factor(igraph::cluster_louvain(graph_10pcs)$membership)

graph_mnn <- buildSNNGraph(ring_hard, k = 20, use.dimred="MNN_corrected", type = "jaccard")
ring_hard$cluster_mnn <- factor(igraph::cluster_louvain(graph_mnn)$membership)


# visualise distribution of cells in clusters and batches
p1 <- plotReducedDim(ring_hard, "UMAPall_30",
               colour_by = "cluster_50pcs", text_by = "cluster_50pcs") +
  theme(legend.position = "none") +
  labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP on 50 PCs") +
  scale_fill_viridis_d()

temp <- as.data.table(colData(ring_hard))
temp <- temp[, .(n = .N), by = c("cluster_50pcs", "Sample")]
temp <- temp[, pct := n/sum(n), by = "Sample"]
p2 <- ggplot(temp, aes(Sample, cluster_50pcs)) +
  geom_point(aes(size = pct))

p3 <- plotReducedDim(ring_hard, "UMAP_30",
                     colour_by = "cluster_10pcs", text_by = "cluster_10pcs") +
  theme(legend.position = "none") +
  labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP on 10 PCs") +
  scale_fill_viridis_d()

temp <- as.data.table(colData(ring_hard))
temp <- temp[, .(n = .N), by = c("cluster_10pcs", "Sample")]
temp <- temp[, pct := n/sum(n), by = "Sample"]
p4 <- ggplot(temp, aes(Sample, cluster_10pcs)) +
  geom_point(aes(size = pct))

p5 <- plotReducedDim(ring_hard, "UMAP-MNN_30",
                     colour_by = "cluster_mnn", text_by = "cluster_mnn") +
  theme(legend.position = "none") +
  labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP on MNN") +
  scale_fill_viridis_d()

temp <- as.data.table(colData(ring_hard))
temp <- temp[, .(n = .N), by = c("cluster_mnn", "Sample")]
temp <- temp[, pct := n/sum(n), by = "Sample"]
p6 <- ggplot(temp, aes(Sample, cluster_mnn)) +
  geom_point(aes(size = pct))

(p1 + p2) / (p3 + p4) / (p5 + p6) &
  theme(legend.position = "none")


# Batch effects - soft filter ---------------------------------------------

# UMAP on original and corrected data
p1 <- ggplot(getReducedDim(ring_soft, type = "UMAPall_30"),
             aes(V1, V2)) +
  geom_point(aes(colour = Sample), alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "No batch correction",
       caption = "UMAP ran on all PCs")

pcvar <- data.table(pc = 1:50,
                    var = attr(reducedDim(ring_soft, "PCA"), "percentVar"))
p2 <- ggplot(pcvar, aes(pc, var)) +
  geom_col() +
  geom_point(aes(y = cumsum(var))) +
  geom_line(aes(y = cumsum(var))) +
  geom_vline(xintercept = 10, colour = "dodgerblue", linetype = 2) +
  labs(x = "PC", y = "% Variance") +
  scale_y_continuous(limits = c(0, 80))

p3 <- ggplot(getReducedDim(ring_soft, type = "UMAP_30"),
             aes(V1, V2)) +
  geom_point(aes(colour = Sample), alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "No batch correction",
       caption = "UMAP ran on first 10 PCs")

p4 <- ggplot(getReducedDim(ring_soft, type = "UMAP-MNN_30"),
             aes(V1, V2)) +
  geom_point(aes(colour = Sample), alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "Batch correction",
       caption = "UMAP ran on MNN-corrected data")

(p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")


# highlight those cells with low UMIs - batch correction has large effect on there!
p1 <- ggplot(getReducedDim(ring_soft, type = "UMAPall_30"),
             aes(V1, V2)) +
  geom_point(aes(colour = total)) +
  scale_colour_viridis_c(trans = "log10") +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "No batch correction",
       caption = "UMAP ran 50 PCs")

p2 <- ggplot(getReducedDim(ring_soft, type = "UMAP-MNN_30"),
             aes(V1, V2)) +
  geom_point(aes(colour = total)) +
  scale_colour_viridis_c(trans = "log10") +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "Batch correction",
       caption = "UMAP ran on MNN-corrected data")

p3 <- ggplot(getReducedDim(ring_soft, type = "UMAPall_30"),
       aes(V1, V2)) +
  geom_point(aes(colour = total < 1000)) +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "No batch correction",
       caption = "UMAP ran 50 PCs")

p4 <- ggplot(getReducedDim(ring_soft, type = "UMAP-MNN_30"),
       aes(V1, V2)) +
  geom_point(aes(colour = total < 1000)) +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "Batch correction",
       caption = "UMAP ran on MNN-corrected data")

(p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")



# Marker genes ------------------------------------------------------------

# visualise expression of marker genes in different datasets
dat <- getReducedDim(ring_hard, "UMAP-MNN_30",
                     genes = unique(markers$id), melted = TRUE)
dat <- merge(dat, markers, by = "id")
dat$expr <- ifelse(dat$expr == 0, NA, dat$expr)
ggplot(dat[order(expr, na.last = FALSE), ], aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  scale_colour_viridis_c() +
  facet_wrap( ~ name) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  coord_fixed()

# Same but with soft-filtered data
dat <- getReducedDim(ring_soft, "UMAP-MNN_30",
                     genes = unique(markers$id), melted = TRUE)
dat <- merge(dat, markers, by = "id")
dat$expr <- ifelse(dat$expr == 0, NA, dat$expr)
ggplot(dat[order(expr, na.last = FALSE), ], aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  scale_colour_viridis_c() +
  facet_wrap( ~ name) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  coord_fixed()

# Look at correlation of each gene and total UMIs
dat <- getReducedDim(ring_soft, "UMAP-MNN_30",
                     genes = unique(markers$id), melted = TRUE)
dat <- merge(dat, markers, by = "id")
ggplot(dat,
       aes(total, expr)) +
  geom_pointdensity(show.legend = FALSE, size = 0.3) +
  facet_wrap(~ name) +
  scale_x_log10() +
  scale_colour_viridis_c(option = "magma")


# Some random sample of genes
set.seed(1588669631)
dat <- getReducedDim(ring_soft, "UMAP-MNN_30",
                     genes = sample(rownames(ring_soft), 9), melted = TRUE)
ggplot(dat,
       aes(total, expr)) +
  geom_pointdensity(show.legend = FALSE, size = 0.3) +
  facet_wrap(~ id) +
  scale_x_log10() +
  scale_colour_viridis_c(option = "magma") +
  labs(x = "Total UMIs", y = "logcounts", title = "Random Sample of Genes")



# Figure ------------------------------------------------------------------
# a figure summarising all the QC issues and justification for removing low-UMI clusters

# Total UMI distributions
ggplot(colData(ring_soft),
       aes(total, Sample)) +
  ggridges::geom_density_ridges(alpha = 0.5) +
  scale_x_log10(labels = scales::label_number_si()) +
  annotation_logticks(sides = "b") +
  labs(x = "Total UMIs", y = "")

# detected genes
ggplot(colData(ring_soft),
       aes(detected, Sample)) +
  ggridges::geom_density_ridges(alpha = 0.5) +
  scale_x_log10(labels = scales::label_number_si()) +
  annotation_logticks(sides = "b") +
  labs(x = "No. detected genes", y = "")

# UMAP coloured by total counts
ggplot(getReducedDim(ring_soft, "UMAP-MNN_30"),
       aes(V1, V2)) +
  geom_point(aes(colour = total)) +
  geom_label(stat = "centroid", aes(group = cluster_mnn, label = cluster_mnn),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  scale_colour_viridis_c(trans = "log10") +
  labs(x = "UMAP1", y = "UMAP2", colour = "Total UMIs")

# PCA showing high correlation with PC2
ggplot(getReducedDim(ring_soft, "PCA"),
       aes(PC1, PC2)) +
  geom_point(aes(colour = total)) +
  scale_colour_viridis_c(trans = "log10",
                         labels = scales::label_number_si()) +
  labs(colour = "Total UMIs",
       x = paste0("PC1 (", round(attr(reducedDim(ring_soft, "PCA"), "percentVar")[1]), "%)"),
       y = paste0("PC1 (", round(attr(reducedDim(ring_soft, "PCA"), "percentVar")[2]), "%)"))
