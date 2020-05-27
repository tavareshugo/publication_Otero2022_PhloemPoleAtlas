# Package setup -----------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(DropletUtils)
  library(scater)
  library(scran)
  library(batchelor)
  # library(zinbwave)
  # library(glmpca)
})

set.seed(1001)


# Capture user input ------------------------------------------------------

option_list = list(
  make_option(c("--min_total_umi_per_cell"),
              action = "store",
              default = 1,
              type = "integer",
              help = "Minimum total UMI to retain a cell."),
  make_option(c("--min_genes_per_cell"),
              action = "store",
              default = 1,
              type = "integer",
              help = "Minimum detected genes to retain a cell."),
  make_option(c("--min_cells_gene_detected_in"),
              action = "store",
              default = 1,
              type = "integer",
              help = "Minimum number of cells that a gene has to be expressed in to be retained."),
  make_option(c("--min_gene_counts"),
              action = "store",
              default = 1,
              type = "integer",
              help = "Minimum number of counts to consider a gene as expressed for the --min_gene_detected filter."),
  make_option(c("--max_mito_pct"),
              action = "store",
              default = 10,
              type = "integer",
              help = "Maximum percentage of counts assigned to mitochondrial genes allowed to retain a cell."),
  make_option(c("--n_hvgs"),
              action = "store",
              default = 1000,
              type = "integer",
              help = "Number of highly variable genes to use for PCA (and stored in metadata)."),
  make_option(c("--out"),
              action = "store",
              default = NA,
              type = 'character',
              help = "Output file name."),
  make_option("--cores",
              action = "store",
              default = NA,
              type = "integer",
              help = "Number of cores to use. Will use all available by default.")
)
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = c(1, Inf))

options <- opt$options
files <- opt$args

# # for testing
# options$min_total_umi_per_cell <- 1000
# options$min_cells_gene_detected_in <- 5
# options$min_gene_counts <- 5
# options$max_mito_pct <- 10
# options$n_hvgs <- 2000
# options$out <- "data/processed/SingleCellExperiment/test.rds"
# options$cores <- 2
# files <- list.files(".", "_sce.rds", recursive = TRUE, full.names = TRUE)
# files <- files[!grepl("denyer", files)]

# check inputs
if(any(!file.exists(files))) stop("Some input files do not exist!")

if(is.na(options$cores) | options$cores > parallel::detectCores()){
  options$cores <- BiocParallel::MulticoreParam(parallel::detectCores())
} else {
  options$cores <- BiocParallel::MulticoreParam(options$cores)
}


# Read data ---------------------------------------------------------------

# read all files in
sce_all <- lapply(
  files,
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

    message("  - With ", nrow(i), " genes and ", ncol(i), " cells.")
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
message("Combined data before filtering: ", nrow(sce_all), " genes and ", ncol(sce_all), " cells.")

#### Cell filtering
# mitochondrial percent filter
sce_all <- sce_all[, which(sce_all$subsets_mitochondria_percent <= options$max_mito_pct)]

# retain cells with total counts >5% of 99th quantile
# sce_all <- sce_all[, which((sce_all$total) >= (quantile((sce_all$total), 0.99) * 0.05))]

# filter on total UMIs
sce_all <- sce_all[, which(sce_all$total >= options$min_total_umi_per_cell)]

# filter on detected genes
sce_all[, which(colSums(counts(sce_all) >= options$min_gene_counts) >= options$min_genes_per_cell)]

# rescale logcounts within each batch (to account for different library sizes)
sce_all <- multiBatchNorm(sce_all, batch = sce_all$Sample)
metadata(sce_all) <- list() # reset metadata, because multiBatchNorm makes it weird


#### filter genes
# remove genes that don't pass minimal thresholds
sce_all <- sce_all[which(rowSums(counts(sce_all) >= options$min_gene_counts) >= options$min_cells_gene_detected_in), ]

# get ID of nuclear genes
metadata(sce_all)$nucgenes <- rownames(sce_all)[grep("AT[1,2,3,4,5]", rownames(sce_all))]

message("After filtering: ", nrow(sce_all), " genes and ", ncol(sce_all), " cells.")


# Highly variable genes ---------------------------------------------------

# get the estimated variance for each gene and top "highly variable genes"
metadata(sce_all)[["genevar"]] <- modelGeneVar(sce_all, block = colData(sce_all)$Sample,
                                               subset.row = metadata(sce_all)$nucgenes)
metadata(sce_all)[["hvgs"]] <- getTopHVGs(metadata(sce_all)$genevar,
                                          n = options$n_hvgs)


# MNN correction ----------------------------------------------------------
message("Applying MNN batch correction...")

# create object with MNN correction
sce_all_mnn <- fastMNN(sce_all, batch = colData(sce_all)$Sample,
                       subset.row = metadata(sce_all)$nucgenes,
                       BPPARAM = options$cores)
reducedDim(sce_all, "MNN_corrected") <- reducedDim(sce_all_mnn, "corrected")

rm(sce_all_mnn)


# # ZINB-WAVE correction ----------------------------------------------------
# message("Applying ZINB-WAVE correction...")
#
# # getting the error reported in this issue: https://github.com/drisso/zinbwave/issues/45
# # getting arount it by enforcing a matrix
# assay(sce_all, "counts_matrix") <- as.matrix(counts(sce_all))
#
# # We do this on the ring data only (as it's the focus of this work)
# sce_all <- zinbwave(sce_all,
#                     X = ~ Sample,
#                     K = 10,
#                     zeroinflation = FALSE, verbose = TRUE,
#                     which_assay = "counts_matrix",
#                     which_genes = metadata(sce_all)$hvgs,
#                     normalizedValues = TRUE,
#                     residuals = TRUE,
#                     observationalWeights = TRUE,
#                     BPPARAM = BiocParallel::MulticoreParam(3))
#
# assay(sce_all, "counts_matrix") <- NULL   # remove it from the object

# # glmpca
# metadata(sce_ring)$glmpca <- glmpca(counts(sce_ring[metadata(sce_ring)$nucgenes, ]),
#                                     L = 10,
#                                     X = cbind(Sample = sce_ring$Sample),
#                                     fam = "nb")

# sctransform normalisation -----------------------------------------------

vst_norm <- vst(umi = counts(sce_all),
                cell_attr = colData(sce_all),
                batch_var = "Sample",
                return_corrected_umi = TRUE)

assay(sce_all, "SCT") <- vst_norm$umi_corrected
assay(sce_all, "SCTlog1p") <- log1p(assay(sce_all, "SCT"))


# Dimensionality reduction ----------------------------------------------

# PCA
reducedDim(sce_all, "PCA") <- calculatePCA(sce_all,
                                           subset_row = metadata(sce_all)$hvgs)

# UMAP with different n_neighbors (uwot::umap default is 15)
# number of neighbours refers to how many neighbours are used for the computation
# higher values for a more global picture and smaller for more local picture
for(i in c(7, 15, 30, 100)){
  message("UMAP with n_neighbors = ", i)

  reducedDim(sce_all, paste0("UMAP_", i)) <- calculateUMAP(sce_all,
                                                           dimred = "PCA",
                                                           n_dimred = 10,
                                                           n_neighbors = i)
  reducedDim(sce_all, paste0("UMAPall_", i)) <- calculateUMAP(sce_all,
                                                              dimred = "PCA",
                                                              n_neighbors = i)
  reducedDim(sce_all, paste0("UMAP-MNN_", i)) <- calculateUMAP(sce_all,
                                                               dimred = "MNN_corrected",
                                                                n_neighbors = i)
  # reducedDim(sce_all, paste0("UMAP-SCT_", i)) <- calculateUMAP(sce_all,
  #                                                              dimred = "SCTlog1p",
  #                                                              n_neighbors = i)
}

# t-SNE with different perplexity (Rtsne::Rtsne default is 30)
# perplexity refers to variance of gaussian distribution that weights distances in high-dim space
for(i in c(15, 30, 60)){
  message("t-SNE with perplexity = ", i)

  reducedDim(sce_all, paste0("TSNE_", i)) <- calculateTSNE(sce_all, dimred = "PCA",
                                                            n_dimred = 10,
                                                            n_neighbors = i)
  reducedDim(sce_all, paste0("TSNEall_", i)) <- calculateTSNE(sce_all,
                                                              dimred = "PCA",
                                                               n_neighbors = i)
  reducedDim(sce_all, paste0("TSNE-MNN_", i)) <- calculateTSNE(sce_all,
                                                               dimred = "MNN_corrected",
                                                                n_neighbors = i)
  # reducedDim(sce_all, paste0("TSNE-ZINB_", i)) <- calculateTSNE(sce_all,
  #                                                              dimred = "zinbwave",
  #                                                              n_neighbors = i)
}

# DiffusionMap
message("Diffusion maps...")
reducedDim(sce_all, "DIFFMAP") <- calculateDiffusionMap(sce_all, dimred = "PCA",
                                                        ncomponents = 10,
                                                        n_dimred = 10)

reducedDim(sce_all, "DIFFMAPall") <- calculateDiffusionMap(sce_all,
                                                           ncomponents = 10,
                                                           dimred = "PCA")

reducedDim(sce_all, "DIFFMAP-MNN") <- calculateDiffusionMap(sce_all,
                                                            ncomponents = 10,
                                                            dimred = "MNN_corrected")


# Clustering --------------------------------------------------------------

# Graph-based clustering (similar to Seurat method)
graph_mnn <- buildSNNGraph(sce_all, k = 100, use.dimred = "MNN_corrected",
                                             type = "jaccard")
sce_all$cluster_mnn <- factor(igraph::cluster_louvain(graph_mnn)$membership)

graph_pca <- buildSNNGraph(sce_all, k = 100, use.dimred = "PCA",
                                             type = "jaccard")
sce_all$cluster_pca <- factor(igraph::cluster_louvain(graph_pca)$membership)


# Save objects ------------------------------------------------------------
message("Saving output...")
saveRDS(sce_all, options$out)
