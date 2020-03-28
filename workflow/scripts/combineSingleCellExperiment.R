# Package setup -----------------------------------------------------------

suppressPackageStartupMessages({
  library(DropletUtils)
  library(scater)
  library(scran)
  library(batchelor)
  library(zinbwave)
  # library(glmpca)
})

set.seed(1001)


# Read data ---------------------------------------------------------------

# read all files in
sce_all <- lapply(
  list.files("data/processed/SingleCellExperiment/",
                          pattern = "_sce.rds", recursive = TRUE, full.names = TRUE),
  function(i){
    message("Reading ", i)
    # read data
    i <- readRDS(i)

    # retain rowData columns of interest
    rowData(i) <- rowData(i)[, 1:10]

    # add prefix to cell names (the same barcodes may occur in different batches)
    colnames(i) <- paste0(unique(i$Sample),  "_", colnames(i))

    # remove MT and C genes
    # i <- i[grep("^AT[C,M]", rownames(i), invert = TRUE), ]

    return(i)
    })

# combine samples to single object
sce_all <- do.call("cbind", sce_all)

# remove dimensionality reduction slots
reducedDim(sce_all, "PCA") <- NULL
reducedDim(sce_all, "UMAP") <- NULL
reducedDim(sce_all, "TSNE") <- NULL

# sanity checks
if(!all(gsub(".*_", "", colnames(sce_all)) == colData(sce_all)$Barcode)) stop("colData corruped!")
if(!all(rownames(sce_all) == rowData(sce_all)$ID)) stop("rowData corruped!")


# Filter data -------------------------------------------------------------

#### Cell filtering
# retain cells with <10% MT counts
sce_all <- sce_all[, which(sce_all$subsets_mitochondria_percent < 10)]

# retain cells with total counts >5% of 99th quantile
# sce_all <- sce_all[, which((sce_all$total) >= (quantile((sce_all$total), 0.99) * 0.05))]

# retain cells with >1k total UMIs (this seems OK based on the joint distribution)
sce_all <- sce_all[, which(sce_all$total > 1e3)]

# rescale logcounts within each batch (to account for different library sizes)
sce_all <- multiBatchNorm(sce_all, batch = sce_all$Sample)
metadata(sce_all) <- list() # reset metadata, because multiBatchNorm makes it weird


#### filter genes
# remove genes with no counts at all
sce_all <- sce_all[which(rowSums(counts(sce_all)) > 0),]

# get ID of nuclear genes
metadata(sce_all)$nucgenes <- rownames(sce_all)[grep("AT[1,2,3,4,5]", rownames(sce_all))]


#### ring data
# make separate object with ring cells only
sce_ring <- sce_all[, !grepl("denyer", sce_all$Sample)]

# remove genes with no counts
sce_ring <- sce_ring[which(rowSums(counts(sce_ring)) > 0),]

# get ID of nuclear genes
metadata(sce_ring)$nucgenes <- rownames(sce_ring)[grep("AT[1,2,3,4,5]", rownames(sce_ring))]


# MNN correction ----------------------------------------------------------

# create object with MNN correction
sce_all_mnn <- fastMNN(sce_all, batch = colData(sce_all)$Sample,
                       subset.row = metadata(sce_all)$nucgenes)
# assay(sce_all, "MNN_reconstructed") <- assay(sce_all_mnn, "reconstructed")
reducedDim(sce_all, "MNN_corrected") <- reducedDim(sce_all_mnn, "corrected")

sce_ring_mnn <- fastMNN(sce_ring, batch = colData(sce_ring)$Sample,
                        subset.row = metadata(sce_ring)$nucgenes)
# assay(sce_ring, "MNN_reconstructed") <- assay(sce_ring_mnn, "reconstructed")
reducedDim(sce_ring, "MNN_corrected") <- reducedDim(sce_ring_mnn, "corrected")

rm(sce_all_mnn, sce_ring_mnn)


# Dimensionality reduction - all batches ----------------------------------

# get the estimated variance for each gene and top "highly variable genes"
metadata(sce_all)[["genevar"]] <- modelGeneVar(sce_all, block = colData(sce_all)$Sample,
                                                subset.row = metadata(sce_all)$nucgenes)
metadata(sce_all)[["hvgs"]] <- getTopHVGs(metadata(sce_all)$genevar, n = 2000)

# PCA
reducedDim(sce_all, "PCA") <- calculatePCA(sce_all, subset_row = metadata(sce_all)$hvgs)

# UMAP with different n_neighbors (uwot::umap default is 15)
# number of neighbours refers to how many neighbours are used for the computation
# higher values for a more global picture and smaller for more local picture
for(i in c(7, 15, 30, 100)){
  message("UMAP with n_neighbors = ", i)

  reducedDim(sce_all, paste0("UMAP_", i)) <- calculateUMAP(sce_all, dimred = "PCA",
                                                            n_dimred = 10,
                                                            n_neighbors = i)
  reducedDim(sce_all, paste0("UMAPall_", i)) <- calculateUMAP(sce_all, dimred = "PCA",
                                                               n_neighbors = i)
  reducedDim(sce_all, paste0("UMAP-MNN_", i)) <- calculateUMAP(sce_all, dimred = "MNN_corrected",
                                                                n_neighbors = i)
}

# t-SNE with different perplexity (Rtsne::Rtsne default is 30)
# perplexity refers to variance of gaussian distribution that weights distances in high-dim space
for(i in c(15, 30, 60)){
  message("t-SNE with perplexity = ", i)

  reducedDim(sce_all, paste0("TSNE_", i)) <- calculateTSNE(sce_all, dimred = "PCA",
                                                            n_dimred = 10,
                                                            n_neighbors = i)
  reducedDim(sce_all, paste0("TSNEall_", i)) <- calculateTSNE(sce_all, dimred = "PCA",
                                                               n_neighbors = i)
  reducedDim(sce_all, paste0("TSNE-MNN_", i)) <- calculateTSNE(sce_all, dimred = "MNN_corrected",
                                                                n_neighbors = i)
}

# DiffusionMap
reducedDim(sce_all, "DIFFMAP") <- calculateDiffusionMap(sce_all, dimred = "PCA",
                                                         n_dimred = 10)
reducedDim(sce_all, "DIFFMAPall") <- calculateDiffusionMap(sce_all, dimred = "PCA")
reducedDim(sce_all, "DIFFMAP-MNN") <- calculateDiffusionMap(sce_all,
                                                             dimred = "MNN_corrected")


# Dimensionality reduction - ring batches ---------------------------------

# get the estimated variance for each gene and top "highly variable genes"
metadata(sce_ring)[["genevar"]] <- modelGeneVar(sce_ring, block = colData(sce_ring)$Sample,
                                                subset.row = metadata(sce_ring)$nucgenes)
metadata(sce_ring)[["hvgs"]] <- getTopHVGs(metadata(sce_ring)$genevar, n = 2000)

# PCA
reducedDim(sce_ring, "PCA") <- calculatePCA(sce_ring, subset_row = metadata(sce_ring)$hvgs)

# UMAP with different n_neighbors (uwot::umap default is 15)
# number of neighbours refers to how many neighbours are used for the computation
# higher values for a more global picture and smaller for more local picture
for(i in c(7, 15, 30, 100)){
  message("UMAP with n_neighbors = ", i)

  reducedDim(sce_ring, paste0("UMAP_", i)) <- calculateUMAP(sce_ring, dimred = "PCA",
                                                            n_dimred = 10,
                                                            n_neighbors = i)
  reducedDim(sce_ring, paste0("UMAPall_", i)) <- calculateUMAP(sce_ring, dimred = "PCA",
                                                            n_neighbors = i)
  reducedDim(sce_ring, paste0("UMAP-MNN_", i)) <- calculateUMAP(sce_ring, dimred = "MNN_corrected",
                                                               n_neighbors = i)
}

# t-SNE with different perplexity (Rtsne::Rtsne default is 30)
# perplexity refers to variance of gaussian distribution that weights distances in high-dim space
for(i in c(15, 30, 60)){
  message("t-SNE with perplexity = ", i)

  reducedDim(sce_ring, paste0("TSNE_", i)) <- calculateTSNE(sce_ring, dimred = "PCA",
                                                            n_dimred = 10,
                                                            n_neighbors = i)
  reducedDim(sce_ring, paste0("TSNEall_", i)) <- calculateTSNE(sce_ring, dimred = "PCA",
                                                               n_neighbors = i)
  reducedDim(sce_ring, paste0("TSNE-MNN_", i)) <- calculateTSNE(sce_ring, dimred = "MNN_corrected",
                                                                n_neighbors = i)
}

# DiffusionMap
reducedDim(sce_ring, "DIFFMAP") <- calculateDiffusionMap(sce_ring, dimred = "PCA",
                                                         n_dimred = 10)
reducedDim(sce_ring, "DIFFMAPall") <- calculateDiffusionMap(sce_ring, dimred = "PCA")
reducedDim(sce_ring, "DIFFMAP-MNN") <- calculateDiffusionMap(sce_ring,
                                                             dimred = "MNN_corrected")


# ZINB-WAVE correction ----------------------------------------------------

# We do this on the ring data only (as it's the focus of this work)
sce_wave <- zinbwave(sce_ring,
                     X = ~ Sample,
                     K = 10,
                     zeroinflation = FALSE, verbose = TRUE,
                     which_assay = "counts", which_genes = metadata(sce_ring)$hvgs,
                     normalizedValues = TRUE, residuals = TRUE, observationalWeights = TRUE,
                     BPPARAM = BiocParallel::MulticoreParam(parallel::detectCores()))

# # Try glmpca
# metadata(sce_ring)$glmpca <- glmpca(counts(sce_ring[metadata(sce_ring)$nucgenes, ]),
#                                     L = 10,
#                                     X = cbind(Sample = sce_ring$Sample),
#                                     fam = "nb")


# Save objects ------------------------------------------------------------

saveRDS(sce_wave, "data/processed/SingleCellExperiment/ring_batches_wave.rds")
saveRDS(sce_all, "data/processed/SingleCellExperiment/all_batches.rds")
saveRDS(sce_ring, "data/processed/SingleCellExperiment/ring_batches.rds")


