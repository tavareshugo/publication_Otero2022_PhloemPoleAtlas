# Package setup -----------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(DropletUtils)
  library(scater)
  library(scran)
  library(batchelor)
  library(sctransform)
  # library(zinbwave)
  # library(glmpca)
})

set.seed(1590572757)


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
              default = 100,
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
# opt <- parse_args(OptionParser(option_list=option_list))
# options <- opt
# options$min_total_umi_per_cell <- 1
# options$min_genes_per_cell <- 300
# options$min_cells_gene_detected_in <- 3
# options$min_gene_counts <- 1
# options$max_mito_pct <- 10
# options$n_hvgs <- 1000
# options$out <- "data/processed/SingleCellExperiment/test.rds"
# options$cores <- 2
# files <- list.files("data/processed/SingleCellExperiment/", "_sce.rds", recursive = TRUE, full.names = TRUE)
# files <- files[!grepl("denyer", files)]
# files <- files[!grepl("mutant|wt", files)]

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

# remove non-nuclear genes
metadata(sce_all)$nucgenes <- rownames(sce_all)[grep("AT[1,2,3,4,5]", rownames(sce_all))]
sce_all <- sce_all[metadata(sce_all)$nucgenes, ]

# mitochondrial percent filter
sce_all <- sce_all[, which(sce_all$subsets_mitochondria_percent <= options$max_mito_pct)]

# iterate through these filters until convergence
before <- 1; after <- 0 # to initiate the while loop
while (sum(before - after) > 0) {
  # rows and cols before filtering
  before <- dim(sce_all)

  # filter on total UMIs
  sce_all <- sce_all[, which(colSums(counts(sce_all)) >= options$min_total_umi_per_cell)]

  # filter on detected genes
  sce_all <- sce_all[, which(colSums(counts(sce_all) >= options$min_gene_counts) >= options$min_genes_per_cell)]

  # remove genes that don't pass minimal thresholds
  sce_all <- sce_all[which(rowSums(counts(sce_all) >= options$min_gene_counts) >= options$min_cells_gene_detected_in), ]

  # rows and cols after filtering
  after <- dim(sce_all)
}

# rescale logcounts within each batch (to account for different library sizes)
sce_all <- multiBatchNorm(sce_all, batch = sce_all$Sample)
metadata(sce_all) <- list() # reset metadata, because multiBatchNorm makes it weird

message("After filtering: ", nrow(sce_all), " genes and ", ncol(sce_all), " cells.")


# Highly variable genes ---------------------------------------------------

# get the estimated variance for each gene and top "highly variable genes"
metadata(sce_all)[["genevar"]] <- modelGeneVar(sce_all,
                                               block = colData(sce_all)$Sample,
                                               subset.row = metadata(sce_all)$nucgenes)
metadata(sce_all)[["hvgs"]] <- getTopHVGs(metadata(sce_all)$genevar,
                                          n = options$n_hvgs)


# sctransform normalisation -----------------------------------------------

# variance-stabilising transform
# (setting min_cells = 1 as the filtering is done earlier)
vst_norm <- vst(umi = counts(sce_all),
                cell_attr = colData(sce_all),
                batch_var = "Sample",
                min_cells = 1,
                return_cell_attr = TRUE,
                return_corrected_umi = TRUE)

# add the model matrices to assay slot
assay(sce_all, "vst") <- vst_norm$umi_corrected
assay(sce_all, "logvst") <- log1p(assay(sce_all, "vst"))

# get the most variable (nuclear) genes according to this model
temp <- vst_norm$gene_attr[metadata(sce_all)$nucgenes, ]
metadata(sce_all)[["hvgs_vst"]] <- rownames(temp)[order(-temp$residual_variance)][1:options$n_hvgs]

rm(temp, vst_norm) # free-up memory


# MNN correction ----------------------------------------------------------
message("Applying MNN batch correction...")

# create object with MNN correction - on simple logcounts
sce_all_mnn <- fastMNN(sce_all,
                       batch = colData(sce_all)$Sample,
                       subset.row = metadata(sce_all)$nucgenes,
                       assay.type = "logcounts",
                       BPPARAM = options$cores)
reducedDim(sce_all, "MNN_logcounts") <- reducedDim(sce_all_mnn, "corrected")

# create object with MNN correction - on VST logcounts
sce_all_mnn <- fastMNN(sce_all,
                       batch = colData(sce_all)$Sample,
                       subset.row = metadata(sce_all)$nucgenes,
                       assay.type = "logvst",
                       BPPARAM = options$cores)
reducedDim(sce_all, "MNN_logvst") <- reducedDim(sce_all_mnn, "corrected")

rm(sce_all_mnn) # free-up memory


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


# Dimensionality reduction ----------------------------------------------
message("Running PCA...")

# PCA - on simple logcounts
reducedDim(sce_all, "PCA_logcounts") <- calculatePCA(sce_all,
                                           exprs_values = "logcounts",
                                           subset_row = metadata(sce_all)$hvgs)

# identify elbow and PC which explains at least 70% variance
# pc_variance <- attr(reducedDim(sce_all, "PCA_logcounts"), "percentVar")
# attr(reducedDim(sce_all, "PCA_logcounts"), "elbow") <- PCAtools::findElbowPoint(pc_variance)
# attr(reducedDim(sce_all, "PCA_logcounts"), "pc70") <-  min(which(cumsum(pc_variance) >= 70))

# PCA - on VST logcounts
reducedDim(sce_all, "PCA_logvst") <- calculatePCA(sce_all,
                                                  exprs_values = "logvst",
                                                  subset_row = metadata(sce_all)$hvgs_vst)
# pc_variance <- attr(reducedDim(sce_all, "PCA_logvst"), "percentVar")
# attr(reducedDim(sce_all, "PCA_logvst"), "elbow") <- PCAtools::findElbowPoint(pc_variance)
# attr(reducedDim(sce_all, "PCA_logvst"), "pc70") <-  min(which(cumsum(pc_variance) >= 70))

# UMAP with different n_neighbors (uwot::umap default is 15)
# number of neighbours refers to how many neighbours are used for the computation
# higher values for a more global picture and smaller for more local picture
set.seed(1590572757)
for(i in c(7, 15, 30, 100)){
  message("UMAP with n_neighbors = ", i)

  # using 10 PCs as we do not expect more than 10 cell types
  # in fact probably less than that, so this should include majority of biological signal
  reducedDim(sce_all, paste0("UMAP", i, "_logcounts")) <-
    calculateUMAP(sce_all,
                  dimred = "PCA_logcounts",
                  n_dimred = 10,
                  n_neighbors = i)

  reducedDim(sce_all, paste0("UMAP", i, "_logvst")) <-
    calculateUMAP(sce_all,
                  dimred = "PCA_logvst",
                  n_dimred = 10,
                  n_neighbors = i)

  # MNN-corrected data on logcounts
  reducedDim(sce_all, paste0("UMAP", i, "_MNN_logcounts")) <-
    calculateUMAP(sce_all,
                  dimred = "MNN_logcounts",
                  n_neighbors = i)

  # MNN-corrected data on logvst
  reducedDim(sce_all, paste0("UMAP", i, "_MNN_logvst")) <-
    calculateUMAP(sce_all,
                  dimred = "MNN_logvst",
                  n_neighbors = i)
}

# t-SNE with different perplexity (Rtsne::Rtsne default is 30)
# perplexity refers to variance of gaussian distribution that weights distances in high-dim space
set.seed(1590572757)
for(i in c(15, 30, 60)){
  message("t-SNE with perplexity = ", i)

  # using PCs explaining at least 70% of variance
  reducedDim(sce_all, paste0("TSNE", i, "_logcounts")) <-
    calculateTSNE(sce_all,
                  dimred = "PCA_logcounts",
                  n_dimred = 10,
                  n_neighbors = i)

  reducedDim(sce_all, paste0("TSNE", i, "_logvst")) <-
    calculateTSNE(sce_all,
                  dimred = "PCA_logvst",
                  n_dimred = 10,
                  n_neighbors = i)

  # MNN-corrected data on logcounts
  reducedDim(sce_all, paste0("TSNE", i, "_MNN_logcounts")) <-
    calculateTSNE(sce_all,
                  dimred = "MNN_logcounts",
                  n_neighbors = i)

  # MNN-corrected data on logvst
  reducedDim(sce_all, paste0("TSNE", i, "_MNN_logvst")) <-
    calculateTSNE(sce_all,
                  dimred = "MNN_logvst",
                  n_neighbors = i)
}

# DiffusionMap
set.seed(1590572757)
message("Diffusion maps...")
reducedDim(sce_all, "DIFFMAP_logcounts") <-
  calculateDiffusionMap(sce_all,
                        dimred = "PCA_logcounts",
                        ncomponents = 10,
                        n_dimred = 10)

reducedDim(sce_all, "DIFFMAP_logvst") <-
  calculateDiffusionMap(sce_all,
                        dimred = "PCA_logvst",
                        ncomponents = 10,
                        n_dimred = 10)

reducedDim(sce_all, "DIFFMAP_MNN_logcounts") <-
  calculateDiffusionMap(sce_all,
                        ncomponents = 10,
                        dimred = "MNN_logcounts")

reducedDim(sce_all, "DIFFMAP_MNN_logvst") <-
  calculateDiffusionMap(sce_all,
                        ncomponents = 10,
                        dimred = "MNN_logvst")


# Clustering --------------------------------------------------------------
set.seed(1590572757)

# Graph-based clustering (essentially similar to Seurat's implementation)
for(i in c("PCA", "MNN")){
  for(ii in c("logcounts", "logvst")){
    message("clustering with ", paste(i, ii, sep = "_"))

    # compute graph
    graph <- buildSNNGraph(sce_all, k = 100, type = "jaccard",
                           use.dimred = paste(i, ii, sep = "_"))

    # add cluster membership to colData
    sce_all[[tolower(paste("cluster", i, ii, sep = "_"))]] <-
      factor(igraph::cluster_louvain(graph)$membership)
  }
}


# Save objects ------------------------------------------------------------
message("Saving output...")
saveRDS(sce_all, options$out)
