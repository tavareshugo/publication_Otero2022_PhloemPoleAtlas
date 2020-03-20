# Package setup -----------------------------------------------------------

suppressPackageStartupMessages({
  library(DropletUtils)
  library(scater)
  library(scran)
  library(batchelor)
})

set.seed(1001)


# Read data ---------------------------------------------------------------

# read all files in
sce_all <- lapply(
  list.files("data/processed/SingleCellExperiment/",
                          pattern = ".rds", recursive = TRUE, full.names = TRUE),
  function(i){
    message("Reading ", i)
    # read data
    i <- readRDS(i)

    # retain rowData columns of interest
    rowData(i) <- rowData(i)[, 1:10]

    # remove MT and C genes
    i <- i[grep("^AT[C,M]", rowData(i)$ID, invert = TRUE), ]
  })

# combine samples to single object
sce_all <- do.call("cbind", sce_all)


# Filter data -------------------------------------------------------------

# remove genes with no counts at all
sce_all <- sce_all[which(rowSums(counts(sce_all)) > 0),]

# retain cells with <10% MT counts
sce_all <- sce_all[, which(colData(sce_all)$subsets_mitochondria_percent < 10)]

# retain cells with total counts >5% of 99th quantile
# sce_all <- sce_all[, which((sce_all$total) >= (quantile((sce_all$total), 0.99) * 0.05))]

# retain cells with >1k total UMIs (this seems OK based on the distribution)
sce_all <- sce_all[, which(sce_all$total > 1e3)]

# rescale logcounts within each batch (to account for different library sizes)
sce_all <- multiBatchNorm(sce_all, batch = sce_all$Sample)

# make separate object with ring cells only
sce_ring <- sce_all[, !grepl("denyer", sce_all$Sample)]


# MNN correction ----------------------------------------------------------

# create object with MNN correction
sce_all_mnn <- fastMNN(sce_all, batch = colData(sce_all)$Sample)
colnames(colData(sce_all_mnn)) <- "Sample"

sce_ring_mnn <- fastMNN(sce_ring, batch = colData(sce_ring)$Sample)
colnames(colData(sce_ring_mnn)) <- "Sample"


# Dimensionality reduction - all batches ----------------------------------

# get the estimated variance for each gene and top "highly variable genes"
metadata(sce_all)[["genevar"]] <- modelGeneVar(sce_all, block = colData(sce_all)$Sample)
metadata(sce_all)[["hvgs"]] <- getTopHVGs(metadata(sce_all)$genevar, n = 2000)

# PCA
reducedDim(sce_all, "PCA") <- calculatePCA(sce_all, subset_row = metadata(sce_all)$hvgs)
reducedDim(sce_all, "PCAall") <- calculatePCA(sce_all)

# UMAP with different n_neighbors (uwot::umap default is 15)
for(i in c(7, 15, 30, 100)){
  message("UMAP with n_neighbors = ", i)

  reducedDim(sce_all, paste0("UMAP_", i)) <- calculateUMAP(sce_all, dimred = "PCA",
                                                           n_neighbors = i)
  reducedDim(sce_all, paste0("UMAPall_", i)) <- calculateUMAP(sce_all,
                                                              n_neighbors = i)
  reducedDim(sce_all_mnn, paste0("UMAP_", i)) <- calculateUMAP(sce_all_mnn, dimred = "corrected",
                                                           n_neighbors = i)
}

# t-SNE with different perplexity (Rtsne::Rtsne default is 30)
for(i in c(15, 30, 60)){
  message("t-SNE with perplexity = ", i)
  reducedDim(sce_all, paste0("TSNE_", i)) <- calculateTSNE(sce_all, dimred = "PCA",
                                                           perplexity = i)
  reducedDim(sce_all, paste0("TSNEall_", i)) <- calculateTSNE(sce_all ,
                                                              perplexity = i)
  reducedDim(sce_all_mnn, paste0("TSNE_", i)) <- calculateTSNE(sce_all_mnn, dimred = "corrected",
                                                           perplexity = i)
}


# Dimensionality reduction - ring batches ---------------------------------

# get the estimated variance for each gene and top "highly variable genes"
metadata(sce_ring)[["genevar"]] <- modelGeneVar(sce_ring, block = colData(sce_ring)$Sample)
metadata(sce_ring)[["hvgs"]] <- getTopHVGs(metadata(sce_ring)$genevar, n = 2000)

# PCA
reducedDim(sce_ring, "PCA") <- calculatePCA(sce_ring, subset_row = metadata(sce_ring)$hvgs)
reducedDim(sce_ring, "PCAall") <- calculatePCA(sce_ring)

# UMAP with different n_neighbors (uwot::umap default is 15)
# number of neighbours refers to how many neighbours are used for the computation
# higher values for a more global picture and smaller for more local picture
for(i in c(7, 15, 30, 100)){
  message("UMAP with n_neighbors = ", i)

  reducedDim(sce_ring, paste0("UMAP_", i)) <- calculateUMAP(sce_ring, dimred = "PCA",
                                                           n_neighbors = i)
  reducedDim(sce_ring, paste0("UMAPall_", i)) <- calculateUMAP(sce_ring,
                                                              n_neighbors = i)
  reducedDim(sce_ring_mnn, paste0("UMAP_", i)) <- calculateUMAP(sce_ring_mnn, dimred = "corrected",
                                                               n_neighbors = i)
}

# t-SNE with different perplexity (Rtsne::Rtsne default is 30)
# perplexity refers to variance of gaussian distribution that weights distances in high-dim space
for(i in c(15, 30, 60)){
  message("t-SNE with perplexity = ", i)
  reducedDim(sce_ring, paste0("TSNE_", i)) <- calculateTSNE(sce_ring, dimred = "PCA",
                                                           perplexity = i)
  reducedDim(sce_ring, paste0("TSNEall_", i)) <- calculateTSNE(sce_ring ,
                                                              perplexity = i)
  reducedDim(sce_ring_mnn, paste0("TSNE_", i)) <- calculateTSNE(sce_ring_mnn, dimred = "corrected",
                                                               perplexity = i)
}



# Save objects ------------------------------------------------------------

saveRDS(sce_all, "data/processed/SingleCellExperiment/all_batches.rds")
saveRDS(sce_all_mnn, "data/processed/SingleCellExperiment/all_batches_mnn.rds")
saveRDS(sce_ring, "data/processed/SingleCellExperiment/ring_batches.rds")
saveRDS(sce_ring_mnn, "data/processed/SingleCellExperiment/ring_batches_mnn.rds")


#######

# wrap_plots(
#   plotReducedDim(sce_all, "UMAP_7", colour_by = "Sample"),
#   plotReducedDim(sce_all, "UMAP_15", colour_by = "Sample"),
#   plotReducedDim(sce_all, "UMAP_30", colour_by = "Sample"),
#   plotReducedDim(sce_all, "UMAP_100", colour_by = "Sample")
# )
#
# wrap_plots(
#   plotReducedDim(sce_mnn, "UMAP_7", colour_by = "dataset"),
#   plotReducedDim(sce_mnn, "UMAP_15", colour_by = "dataset"),
#   plotReducedDim(sce_mnn, "UMAP_30", colour_by = "dataset"),
#   plotReducedDim(sce_mnn, "UMAP_100", colour_by = "dataset")
# )
#
#
# ggplot(cbind(colData(sce_mnn), reducedDim(sce_mnn, "UMAP_7")),
#        aes(V1, V2)) +
#   geom_point(aes(colour = dataset)) +
#   facet_grid(~ dataset)
#
#
# sce_mnn$dataset <- gsub("denyer.*", "Denyer 2019", sce_mnn$Sample)



# Diffusion Map has sigma of a gaussian and also number of nearest neighbours
# destiny::DiffusionMap default sigma = "local"; also try "global"
# default k is determined with some heuristic; try 2000 and 6000?

