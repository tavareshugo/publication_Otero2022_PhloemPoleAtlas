library(data.table)
library(scater)
library(scran)
library(batchelor)
library(ggplot2)
library(patchwork)

# set seed for reproducible results
set.seed(1001)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")


# Read data ---------------------------------------------------------------

# ring data both soft and hard filtered
ring_soft <- readRDS("data/processed/SingleCellExperiment/ring_batches_softfilt.rds")
# ring_hard <- readRDS("data/processed/SingleCellExperiment/ring_batches_hardfilt.rds")


# Find markers ----------------------------------------------

# This runs test for each marker - useful for checking cell types from known genes
markers <- findMarkers(ring_soft,
                       groups = ring_soft$cluster_mnn,
                       block = ring_soft$Sample,
                       test = "wilcox", direction = "up",
                       pval.type = "some", min.prop = 0.9)

# some tidying up
markers <- lapply(markers, function(i){
  as.data.table(i, keep.rownames = "id")
})
markers <- rbindlist(markers, idcol = "cluster", fill = TRUE)
markers[, cluster := factor(as.numeric(cluster))]


# Cell cycle --------------------------------------------------------------

# Fetch cyclin gene IDs
cyclins <- rowData(ring_soft)[grep("CYC[A,B,D]",
                                   rowData(ring_soft)$alternative_name,
                                   ignore.case = TRUE), ]
cyclins <- as.data.table(cyclins) # convert to DT
# clean up the names a bit
cyclins[, alternative_name := gsub(".*/", "", alternative_name)]
cyclins[, cyclin := gsub("^CYC|[0-9].*$", "", alternative_name)]


# Make plot
temp <- getReducedDim(ring_soft, "UMAP-MNN_30", genes = cyclins$ID,
                      melted = TRUE)
temp <- merge(temp, cyclins, by.x = "id", by.y = "ID", all.x = TRUE)

p1 <- ggplot(getReducedDim(ring_soft, "UMAP-MNN_30"), aes(V1, V2, group = cluster_mnn)) +
  geom_point(aes(colour = cluster_mnn), show.legend = FALSE) +
  geom_text(stat = "centroid", aes(label = cluster_mnn), fontface = "bold", size = 4) +
  ggthemes::scale_colour_tableau(palette = "Tableau 20") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  coord_fixed()

p2 <- ggplot(getReducedDim(ring_soft, "UMAP-MNN_30"),
       aes(V1, V2, colour = total < 1000)) +
  geom_point() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  coord_fixed()

p3 <- ggplot(temp, aes(paste(Sample, Barcode), alternative_name)) +
  geom_tile(aes(fill = expr)) +
  facet_grid(cyclin ~ cluster_mnn, scales = "free", space = "free") +
  scale_fill_viridis_c(na.value = "lightgrey") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "Cell", y = "Gene")

(p1 | p2) / p3 + plot_annotation(tag_levels = "A")


# plot test result
ggplot(merge(markers, cyclins, by.x = "id", by.y = "ID"),
             aes(factor(cluster), alternative_name, fill = FDR)) +
  geom_tile() +
  facet_grid(cyclin ~ ., scales = "free", space = "free") +
  scale_fill_gradient(low = "brown", high = "lightgrey", limits = c(0, 1)) +
  labs(x = "Cluster", y = "Gene")


# Cell type annotations ---------------------------------------------------

# read known genes
gene_sets <- readxl::read_xlsx("data/raw/gene_expression_patterns_root_SO.xlsx")
gene_sets <- as.data.table(gene_sets)
gene_sets <- gene_sets[gene_ID %in% rownames(ring_soft), ]

gene_sets <- melt(gene_sets, id.vars = c("gene_name", "gene_name2", "gene_ID", "NOTES"),
                  variable.name = "tissue", value.name = "expressed")
gene_sets <- gene_sets[expressed == "YES"]
gene_sets[, alternative_name := ifelse(is.na(gene_name), gene_ID,
                                       ifelse(is.na(gene_name2), gene_name,
                                              paste(gene_name, gene_name2, sep = "/")))]

# plot test result
ggplot(merge(markers, unique(gene_sets[, .(gene_ID, alternative_name)]),
             by.x = "id", by.y = "gene_ID"),
       aes(cluster, alternative_name, fill = FDR)) +
  geom_tile() +
  # facet_grid(cyclin ~ ., scales = "free", space = "free") +
  scale_fill_gradient(low = "brown", high = "lightgrey", limits = c(0, 1)) +
  labs(x = "Cluster", y = "Gene")

# Make a clustered heatmap with annotations
row_annot <- gene_sets[, .(alternative_name, tissue, expressed)]
row_annot <- dcast(row_annot, alternative_name ~ tissue, value.var = "expressed")
row_annot <- as.matrix(row_annot, rownames = "alternative_name")
row_annot <- row_annot[, c("CC", "PPP", "MSE", "PSE", "PROCAMBIUM", "XPP", "MXY", "PXY", "EPIDERMIS", "CORTEX", "ENDODERMIS", "COLUMELLA", "QC", "LRC")]

temp <- merge(markers, unique(gene_sets[, .(gene_ID, alternative_name)]), by.x = "id", by.y = "gene_ID")
temp <- dcast(temp[, .(alternative_name, cluster, FDR)], alternative_name ~ cluster,
              value.var = "FDR")
temp <- as.matrix(temp, rownames = "alternative_name")
pheatmap::pheatmap(temp,
                   color = colorRampPalette(c("brown", "grey"))(100),
                   annotation_row = as.data.frame(row_annot),
                   annotation_legend = FALSE,
                   treeheight_row = 20, treeheight_col = 20)


# Analysis of subsets -----------------------------------------------------

# Exclude clusters 4 & 8 (cell cycle), 5 (pSE), 7 & 11 (outer layers), and 2 & 6 & 9 (low counts)
# in other words keep only clusters 1, 3, 10
ring_subset <- ring_soft[, which(ring_soft$cluster_mnn %in% c(1, 3, 10))]

# remove old dims
for(i in reducedDimNames(ring_subset)){
  reducedDim(ring_subset, i) <- NULL
  rm(i)
}

# get new set of variable genes
metadata(ring_subset)[["genevar"]] <- modelGeneVar(ring_subset, block = colData(ring_subset)$Sample,
                                                   subset.row = metadata(ring_subset)$nucgenes)
metadata(ring_subset)[["hvgs"]] <- getTopHVGs(metadata(ring_subset)$genevar,
                                              n = 2000)

# create object with MNN correction
ring_subset_mnn <- fastMNN(ring_subset,
                           batch = colData(ring_subset)$Sample,
                           subset.row = metadata(ring_subset)$nucgenes)
reducedDim(ring_subset, "MNN_corrected") <- reducedDim(ring_subset_mnn, "corrected")
rm(ring_subset_mnn); gc()

# Clustering
ring_subset$cluster_mnn_old <- ring_subset$cluster_mnn
graph_subset <-buildSNNGraph(ring_subset, k = 100,
                             use.dimred = "MNN_corrected",
                             type = "jaccard")
ring_subset$cluster_mnn <- factor(igraph::cluster_louvain(graph_subset)$membership)
rm(graph_subset)

# PCA
reducedDim(ring_subset, "PCA") <- calculatePCA(ring_subset,
                                               subset_row = metadata(ring_subset)$hvgs)

# UMAP
reducedDim(ring_subset, "UMAP-MNN_30") <- calculateUMAP(ring_subset,
                                                        dimred = "MNN_corrected",
                                                        n_neighbors = 30)

# visualise
ggplot(getReducedDim(ring_subset, "UMAP-MNN_30"),
       aes(V1, V2)) +
  geom_point(aes(colour = cluster_mnn), show.legend = FALSE) +
  geom_point(stat = "centroid", aes(group = cluster_mnn),
             shape = 21, size = 5, fill = "white", alpha = 0.7) +
  geom_text(stat = "centroid",
            aes(group = cluster_mnn, label = cluster_mnn)) +
  ggthemes::scale_colour_tableau(palette = "Tableau 20") +
  labs(x = "UMAP 1", y = "UMAP 2")

ggplot(colData(ring_subset), aes(cluster_mnn_old, cluster_mnn)) +
  geom_count()


# This runs test for each marker - useful for checking cell types from known genes
markers <- findMarkers(ring_subset, groups = ring_subset$cluster_mnn,
                       block = ring_subset$Sample,
                       test = "wilcox", direction = "up",
                       pval.type = "some", min.prop = 5/7)

# some tidying up
markers <- lapply(markers, function(i){
  as.data.table(i, keep.rownames = "id")
})
markers <- rbindlist(markers, idcol = "cluster", fill = TRUE)
markers[, cluster := factor(as.numeric(cluster))]


# Make a clustered heatmap with annotations
row_annot <- gene_sets[, .(alternative_name, tissue, expressed)]
row_annot <- dcast(row_annot, alternative_name ~ tissue, value.var = "expressed")
row_annot <- as.matrix(row_annot, rownames = "alternative_name")
row_annot <- row_annot[, c("CC", "PPP", "MSE", "PSE", "PROCAMBIUM", "XPP", "MXY", "PXY", "EPIDERMIS", "CORTEX", "ENDODERMIS", "COLUMELLA", "QC", "LRC")]

temp <- merge(markers, unique(gene_sets[, .(gene_ID, alternative_name)]), by.x = "id", by.y = "gene_ID")
temp <- dcast(temp[, .(alternative_name, cluster, FDR)], alternative_name ~ cluster,
              value.var = "FDR")
temp <- as.matrix(temp, rownames = "alternative_name")
pheatmap::pheatmap(temp,
                   color = colorRampPalette(c("brown", "grey"))(100),
                   annotation_row = as.data.frame(row_annot),
                   annotation_legend = FALSE,
                   treeheight_row = 20, treeheight_col = 20)

# even more restricted set and cleaner set
row_annot <- gene_sets[, .(alternative_name, tissue, expressed)]
row_annot <- dcast(row_annot, alternative_name ~ tissue, value.var = "expressed")
row_annot <- as.matrix(row_annot, rownames = "alternative_name")
row_annot <- row_annot[, c("CC", "PPP", "MSE", "PSE")]
# row_annot <- row_annot[, c("CC", "PPP", "MSE", "PSE", "PROCAMBIUM", "XPP", "MXY", "PXY", "EPIDERMIS", "CORTEX", "ENDODERMIS", "COLUMELLA", "QC", "LRC")]

temp <- merge(markers[, .SD[any(FDR < 0.2)], by = id],
              unique(gene_sets[, .(gene_ID, alternative_name)]),
              by.x = "id", by.y = "gene_ID")
temp <- dcast(temp[, .(alternative_name, cluster, FDR)], alternative_name ~ cluster,
              value.var = "FDR")
temp <- as.matrix(temp, rownames = "alternative_name")
pheatmap::pheatmap(temp,
                   color = colorRampPalette(c("brown", "grey"))(100),
                   annotation_row = as.data.frame(row_annot),
                   annotation_legend = FALSE,
                   treeheight_row = 20, treeheight_col = 20)



# Specific informative genes ----------------------------------------------
# Most informative genes (according to microscopy)

# S17 & CALS8 helps identify region of PPP cells
temp <- getReducedDim(ring_soft, "UMAP-MNN_30",
              genes = c("AT2G22850", "AT3G14570"), melted = TRUE)

ggplot(temp[order(expr, na.last = FALSE)], aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  geom_label(stat = "centroid",
            aes(group = cluster_mnn, label = cluster_mnn),
            alpha = 0.5) +
  scale_colour_viridis_c() +
  facet_grid(~ id) +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "Phloem Pole Pericycle markers")

ggplot(temp, aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn, label = cluster_mnn),
             alpha = 0.5) +
  scale_colour_viridis_c() +
  facet_grid(~ id) +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "Phloem Pole Pericycle markers")

ggplot(temp, aes(V1, V2)) +
  geom_point(colour = "lightgrey") +
  ggpointdensity::geom_pointdensity(data = temp[!is.na(expr)]) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn, label = cluster_mnn),
             alpha = 0.5) +
  scale_colour_viridis_c(option = "magma") +
  facet_grid(~ id) +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "Phloem Pole Pericycle markers")

# CC markers
temp <- getReducedDim(ring_soft, "UMAP-MNN_30",
                      genes = c("AT3G12730", "AT5G57350", "AT1G22710"), melted = TRUE)

p1 <- ggplot(temp[order(expr, na.last = FALSE)], aes(V1, V2)) +
  geom_point(aes(colour = expr)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn, label = cluster_mnn),
             alpha = 0.5) +
  scale_colour_viridis_c() +
  facet_grid(~ id) +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "Companion Cell markers")

p2 <- ggplot(temp, aes(V1, V2)) +
  geom_point(colour = "lightgrey") +
  ggpointdensity::geom_pointdensity(data = temp[!is.na(expr)], show.legend = FALSE) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn, label = cluster_mnn),
             alpha = 0.5) +
  scale_colour_viridis_c(option = "magma") +
  facet_grid(~ id) +
  labs(x = "UMAP 1", y = "UMAP 2")

p1 / p2


# Counting approach -------------------------------------------------------

# this is here for completeness, but the approach above seems better
# this is what is used by Denyer et al their Fig 1B

# Get those genes
temp <- getReducedDim(ring_soft, type = "UMAP-MNN_30",
                      genes = unique(gene_sets$gene_ID), melted = TRUE)
temp <- temp[, c("id", "expr", "cluster_mnn")]

# calculate total number of genes in a cluster
temp[, ngenes_cluster := .N, by = .(cluster_mnn)]

# remove cells with no expression for the gene
temp <- temp[expr > 0]

# scale gene expression (z-scores)
temp[, expr_scaled := scale(expr), by = id]

# for each cluster, calculate mean expression and number of expressed genes
temp <- temp[, .(expr_mean = mean(expr),
                 expr_scaled_mean = mean(expr_scaled),
                 nexpr = .N),
             by = .(id, cluster_mnn, ngenes_cluster)]

# merge with annotated table of genes
temp <- merge(temp, unique(gene_sets[, .(alternative_name, gene_ID)]),
              by.x = "id", by.y = "gene_ID", allow.cartesian = TRUE)

ggplot(temp,
       aes(factor(cluster_mnn), alternative_name, size = nexpr/ngenes_cluster*100)) +
  geom_point(aes(colour = expr_scaled_mean)) +
  scale_colour_viridis_c() +
  labs(x = "Cluster", y = "Gene", colour = "Mean\nexpr", size = "Percentage")

