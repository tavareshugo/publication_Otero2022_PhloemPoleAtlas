# library(data.table)
# library(scater)
library(SingleCellExperiment)
library(slingshot)
library(tidyverse)
library(patchwork)

# set seed for reproducible results
set.seed(1001)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")


# Read data ---------------------------------------------------------------

ring_hard <- readRDS("data/processed/SingleCellExperiment/ring_batches_hardfilt.rds")


# Trajectory fitting with slingshot -------------------------------------

# using only first 10 dimensions of MNN (colinearity at higher dimensions?)
# I checked running slingshot with 10 or 40 MNN axis and it's the same essentially
reducedDim(ring_hard, "temp") <- reducedDim(ring_hard, "MNN_logvst")[, 1:10]
sling1 <- slingshot(ring_hard,
                       clusterLabels = ring_hard$cluster_mnn_logvst,
                       reducedDim = "temp")

ring_hard <- slingshot(ring_hard,
                       clusterLabels = ring_hard$cluster_mnn_logvst,
                       reducedDim = "DIFFMAP_MNN_logvst")



# Visualise ---------------------------------------------------------------

temp <- melt(getReducedDim(ring_hard, "UMAP-MNN_30"),
             measure.vars = patterns("^slingPseudotime"),
             variable.name = "trajectory", value.name = "pseudotime")
temp$trajectory <- gsub("slingPseudotime_", "Trajectory ", temp$trajectory)

ring_hard %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  pivot_longer(matches("^slingPseudotime"),
               values_to = "pseudotime", names_to = "trajectory") %>%
  mutate(trajectory = gsub("slingPseudotime_", "Trajectory ", trajectory)) %>%
  arrange(trajectory, !is.na(pseudotime), pseudotime) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = pseudotime), alpha = 0.1) +
  geom_point(stat = "centroid", aes(group = cluster_mnn_logvst), shape = 21, size = 5,
             fill = "white", alpha = 0.7) +
  geom_text(stat = "centroid",
            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst)) +
  facet_wrap(~ trajectory, ncol = 2) +
  scale_colour_viridis_c(option = "inferno", na.value = "lightgrey") +
  coord_fixed() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())


# this is not working quite like I wanted it to...
ggplot(temp,
       aes(x = tidytext::reorder_within(cluster_mnn, pseudotime, trajectory),
           y = pseudotime)) +
  geom_boxplot() +
  facet_wrap(~ trajectory, scales = "free_x") +
  tidytext::scale_x_reordered()



# Expression along pseudotime ---------------------------------------------

# read known genes
gene_sets <- readxl::read_xlsx("data/raw/gene_expression_patterns_root_SO.xlsx")
gene_sets <- as.data.table(gene_sets)

gene_sets <- melt(gene_sets, id.vars = c("gene_name", "gene_name2", "gene_ID", "NOTES"),
                  variable.name = "tissue", value.name = "expressed")
gene_sets <- gene_sets[expressed == "YES"]
gene_sets[, alternative_name := ifelse(is.na(gene_name), gene_ID,
                                       ifelse(is.na(gene_name2), gene_name,
                                              paste(gene_name, gene_name2, sep = "/")))]

# temp <- getReducedDim(ring_hard, "UMAP-MNN_30", genes = gene_sets$gene_ID, melted = FALSE)
pseudotime <- ring_hard$slingPseudotime_1[which(!is.na(ring_hard$slingPseudotime_1))]
temp <- logcounts(ring_hard)[, which(!is.na(ring_hard$slingPseudotime_1))]
temp <- temp[which(rownames(temp) %in% gene_sets$gene_ID), ]
temp <- temp[, order(pseudotime)]

pheatmap::pheatmap(temp, cluster_rows = TRUE, cluster_cols = FALSE,
                   show_colnames = FALSE)



# Focus on cell populations of interest -----------------------------------
library(scran); library(batchelor)

# Exclude 7 & 11 (outer layers), and 2 & 6 & 9 (low counts)
# cluster 1 - PPP
# cluster 3 + 10 - CC
# cluster 5 - pSE
# cluster 4 (& 8?) - cell cycle
ring_subset <- ring_hard[, which(ring_hard$cluster_mnn %in% c(1, 3, 10, 5, 4, 8))]

# # remove old dimreds
# for(i in reducedDimNames(ring_subset)){
#   reducedDim(ring_subset, i) <- NULL
#   rm(i)
# }
#
# # get new set of variable genes
# metadata(ring_subset)[["genevar"]] <- modelGeneVar(ring_subset, block = colData(ring_subset)$Sample,
#                                                subset.row = metadata(ring_subset)$nucgenes)
# metadata(ring_subset)[["hvgs"]] <- getTopHVGs(metadata(ring_subset)$genevar,
#                                           n = 2000)
#
# # create object with MNN correction
# ring_subset_mnn <- fastMNN(ring_subset,
#                            batch = colData(ring_subset)$Sample,
#                            subset.row = metadata(ring_subset)$nucgenes)
# reducedDim(ring_subset, "MNN_corrected") <- reducedDim(ring_subset_mnn, "corrected")
# rm(ring_subset_mnn)

# Clustering
graph_subset <-buildSNNGraph(ring_subset, k = 100,
                             use.dimred = "MNN_corrected",
                             type = "jaccard")
ring_subset$cluster_subset <- factor(igraph::cluster_louvain(graph_subset)$membership)
rm(graph_subset)

# Add new cluster labels
temp <- as.data.table(colData(ring_subset))[
  , .N, by = .(cluster_subset, cluster_mnn)
  ][, .SD[N == max(N)], by = .(cluster_subset)
    ][, cluster := paste(cluster_mnn, cluster_subset, sep = ".")][
        , `:=`(cluster_mnn = NULL, N = NULL)]
x <- temp$cluster
names(x) <- temp$cluster_subset
ring_subset$cluster_subset <- x[as.character(ring_subset$cluster_subset)]

# PCA
reducedDim(ring_subset, "PCA") <- calculatePCA(ring_subset,
                                               subset_row = metadata(ring_subset)$hvgs)

# UMAP
reducedDim(ring_subset, "UMAP-MNN_30") <- calculateUMAP(ring_subset,
                                                        dimred = "MNN_corrected",
                                                        n_neighbors = 30)

# Diffusion map
reducedDim(ring_subset, "DIFFMAP-MNN") <- calculateDiffusionMap(ring_subset,
                                                                   ncomponents = 10,
                                                                dimred = "MNN_corrected")

# visualise
ggplot(getReducedDim(ring_subset, "UMAP-MNN_30"),
       aes(V1, V2)) +
  geom_point(aes(colour = cluster_subset)) +
  geom_label(stat = "centroid",
            aes(group = cluster_subset, label = cluster_subset,
                fill = cluster_subset),
            alpha = 0.8) +
  ggthemes::scale_colour_tableau(palette = "Tableau 20") +
  ggthemes::scale_fill_tableau(palette = "Tableau 20")

ggplot(getReducedDim(ring_subset, "DIFFMAP-MNN"),
       aes(DC3, DC4)) +
  geom_point(aes(colour = cluster_mnn)) +
  geom_point(stat = "centroid", aes(group = cluster_mnn),
             shape = 21, size = 5, fill = "white", alpha = 0.7) +
  geom_text(stat = "centroid",
            aes(group = cluster_mnn, label = cluster_mnn)) +
  ggthemes::scale_colour_tableau(palette = "Tableau 20")


# running on this subset of clusters
reducedDim(ring_subset, "temp") <- reducedDim(ring_subset, "MNN_corrected")[, 1:10]
ring_subset <- slingshot(ring_subset,
                         clusterLabels = ring_subset$cluster_mnn,
                         reducedDim = "temp")

temp <- melt(getReducedDim(ring_subset, "UMAP-MNN_30"),
             measure.vars = patterns("^slingPseudotime"),
             variable.name = "trajectory", value.name = "pseudotime")
temp$trajectory <- gsub("slingPseudotime_", "Trajectory ", temp$trajectory)

ggplot(temp[order(trajectory, pseudotime, na.last = FALSE)],
       aes(V1, V2)) +
  geom_point(aes(colour = pseudotime), alpha = 0.1) +
  geom_label(stat = "centroid",
            aes(group = cluster_mnn, label = cluster_mnn)) +
  facet_wrap(~ trajectory) +
  scale_colour_viridis_c(option = "inferno") +
  coord_fixed() +
  labs(x = "UMAP1", y = "UMAP2",
       subtitle = "Remove outer layers and low-count clusters")


# Enforcing start and end clusters
ring_subset <- slingshot(ring_subset,
                         clusterLabels = ring_subset$cluster_mnn,
                         start.clus = c(4),
                         end.clus = c(5, 1, 3),
                         reducedDim = "temp")

temp <- melt(getReducedDim(ring_subset, "UMAP-MNN_30"),
             measure.vars = patterns("^slingPseudotime"),
             variable.name = "trajectory", value.name = "pseudotime")
temp$trajectory <- gsub("slingPseudotime_", "Trajectory ", temp$trajectory)

ggplot(temp[order(trajectory, pseudotime, na.last = FALSE)],
       aes(V1, V2)) +
  geom_point(aes(colour = pseudotime), alpha = 0.1) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn, label = cluster_mnn)) +
  facet_wrap(~ trajectory) +
  scale_colour_viridis_c(option = "inferno") +
  coord_fixed() +
  labs(x = "UMAP1", y = "UMAP2",
       subtitle = "Enforce start and end clusters")



# Trajectory on individual clusters ---------------------------------------

# Clusters inferred to be PPP
ppp <- runSlingshot(ring_hard[, which(ring_hard$cluster_mnn %in% c(3, 10))],
                    use.dimred = "DIFFMAP-MNN")

# Clusters inferred to be CC
cc <- runSlingshot(ring_hard[, which(ring_hard$cluster_mnn %in% c(1))],
                    use.dimred = "DIFFMAP-MNN")

# Clusters inferred to be MSE/PSE
se <- runSlingshot(ring_hard[, which(ring_hard$cluster_mnn %in% c(5))],
                   use.dimred = "DIFFMAP-MNN")

# Visualise
temp <- melt(getReducedDim(ppp, "UMAP-MNN_30"),
             measure.vars = patterns("^slingPseudotime"),
             variable.name = "trajectory", value.name = "pseudotime")
temp$trajectory <- gsub("slingPseudotime_", "Trajectory ", temp$trajectory)

ggplot(temp[order(trajectory, pseudotime, na.last = FALSE)],
       aes(V1, V2)) +
  geom_point(aes(colour = pseudotime), alpha = 0.1) +
  geom_label(stat = "centroid",
             aes(group = cluster, label = cluster)) +
  facet_wrap(~ trajectory) +
  scale_colour_viridis_c(option = "inferno") +
  coord_fixed() +
  labs(x = "UMAP1", y = "UMAP2")


temp <- melt(getReducedDim(cc, "UMAP-MNN_30"),
             measure.vars = patterns("^slingPseudotime"),
             variable.name = "trajectory", value.name = "pseudotime")
temp$trajectory <- gsub("slingPseudotime_", "Trajectory ", temp$trajectory)

ggplot(temp[order(trajectory, pseudotime, na.last = FALSE)],
       aes(V1, V2)) +
  geom_point(aes(colour = pseudotime), alpha = 0.1) +
  geom_label(stat = "centroid",
             aes(group = cluster, label = cluster)) +
  facet_wrap(~ trajectory) +
  scale_colour_viridis_c(option = "inferno") +
  coord_fixed() +
  labs(x = "UMAP1", y = "UMAP2")


temp <- melt(getReducedDim(se, "UMAP-MNN_30"),
             measure.vars = patterns("^slingPseudotime"),
             variable.name = "trajectory", value.name = "pseudotime")
temp$trajectory <- gsub("slingPseudotime_", "Trajectory ", temp$trajectory)

ggplot(temp[order(trajectory, pseudotime, na.last = FALSE)],
       aes(V1, V2)) +
  geom_point(aes(colour = pseudotime), alpha = 0.1) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn, label = cluster_mnn)) +
  facet_wrap(~ trajectory) +
  scale_colour_viridis_c(option = "inferno") +
  coord_fixed() +
  labs(x = "UMAP1", y = "UMAP2")



# Differential genes along PPP trajectory
library(tradeSeq)
icMat <- evaluateK(counts = counts(ppp), sds = ppp, k = 3:10, nGenes = 200, verbose = TRUE)

temp <- fitGAM(counts = counts(ppp),
               pseudotime = slingPseudotime(ppp, na = FALSE),
               cellWeights = slingCurveWeights(ppp))
