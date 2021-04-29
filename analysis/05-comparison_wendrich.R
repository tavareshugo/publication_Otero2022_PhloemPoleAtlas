# Setup -------------------------------------------------------------------

library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(scran)
library(batchelor)

# options for filtering
options <- list(min_total_umi_per_cell = 1,
                min_genes_per_cell = 2000,
                min_cells_gene_detected_in = 100,
                min_gene_counts = 1,
                max_mito_pct = 10,
                n_hvgs = 2000,
                cores =  BiocParallel::MulticoreParam(2))


# Read data ---------------------------------------------------------------

sce_files <- list.files("data/external/weidrich_et_al/GSE141730_RAW/", full.names = TRUE)
sce_files <- c(
  sce_files,
  paste0("data/intermediate/counts_cellranger/",
         c("APL", "MAKR5", "MAKR5diff", "PEARdel", "S17", "sAPL"),
         "/outs/filtered_feature_bc_matrix.h5")
)

sce_all <- lapply(sce_files,
       function(i){
         message("reading: ", i)
         sce <- read10xCounts(i)
         sce$Sample <- gsub("_filtered.*", "", basename(i))
         colnames(sce) <- paste0(unique(sce$Sample),  "_cell", 1:ncol(sce))

         # make sure gene names are all uppercase
         rownames(sce) <- toupper(rownames(sce))

         # this is for consistency in rowData across datasets
         rowData(sce) <- DataFrame(ID = rownames(sce), Symbol = rownames(sce))
         sce
       })

# subset to genes common to all datasets
common_genes <- Reduce(intersect, lapply(sce_all, rownames))
sce_all <- lapply(sce_all, function(i) i[common_genes, ])
sce_all <- do.call("cbind", sce_all)


# Quality control metrics -----------------------------------------------
message("Adding QC metrics.")

# quality control metrics for cells
sce_all <- addPerCellQC(sce_all,
                    subsets = list(mitochondria = grep("^ATM", rownames(sce_all)),
                                   chloroplast = grep("^ATC", rownames(sce_all))))

# quality control metrics for genes
sce_all <- addPerFeatureQC(sce_all)


# Filter data -------------------------------------------------------------
message("Combined data before filtering: ", nrow(sce_all), " genes and ", ncol(sce_all), " cells.")

# remove non-nuclear genes
metadata(sce_all)$nucgenes <- rownames(sce_all)[grep("AT[1,2,3,4,5]", rownames(sce_all))]
sce_all <- sce_all[metadata(sce_all)$nucgenes, ]

# mitochondrial percent filter
sce_all <- sce_all[, which(sce_all$subsets_mitochondria_percent <= options$max_mito_pct)]

# iterate through these filters until convergence
before <- 1; after <- 0 # to initiate the while loop
i <- 0
while (sum(before - after) > 0) {
  i <- i + 1; message("iteration: ", i)
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

# update nuclear genes ID
metadata(sce_all)$nucgenes <- rownames(sce_all)[grep("AT[1,2,3,4,5]", rownames(sce_all))]

message("After filtering: ", nrow(sce_all), " genes and ", ncol(sce_all), " cells.")


# Highly variable genes ---------------------------------------------------

# get the estimated variance for each gene and top "highly variable genes"
metadata(sce_all)[["genevar"]] <- modelGeneVar(sce_all,
                                               block = colData(sce_all)$Sample,
                                               subset.row = metadata(sce_all)$nucgenes)
metadata(sce_all)[["hvgs"]] <- getTopHVGs(metadata(sce_all)$genevar,
                                          n = options$n_hvgs)



# MNN correction ----------------------------------------------------------
message("Applying MNN batch correction...")

# create object with MNN correction - on simple logcounts
sce_all_mnn <- fastMNN(sce_all,
                       batch = colData(sce_all)$Sample,
                       subset.row = metadata(sce_all)$nucgenes,
                       assay.type = "logcounts")
reducedDim(sce_all, "MNN_logcounts") <- reducedDim(sce_all_mnn, "corrected")

# create object with MNN correction - on VST logcounts
sce_all_mnn <- fastMNN(sce_all,
                       batch = colData(sce_all)$Sample,
                       subset.row = metadata(sce_all)$nucgenes,
                       assay.type = "logvst",
                       BPPARAM = options$cores)
reducedDim(sce_all, "MNN_logvst") <- reducedDim(sce_all_mnn, "corrected")

rm(sce_all_mnn) # free-up memory
