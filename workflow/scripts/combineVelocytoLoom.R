library(LoomExperiment)

# Read and combine loom files ---------------------------------------------

samples <- c("APL", "MAKR5", "MAKR5diff", "PEARdel", "S17", "sAPL")

# Import LOOM files
looms <- lapply(samples, function(sample){
  le <- import(paste0("data/intermediate/velocyto/", sample, ".loom"))

  # Make sure cell IDs are same format as in sce object
  le$CellID <- gsub(paste0(sample, ":"), paste0(sample, "_"), le$CellID)
  le$CellID <- gsub("x$", "-1", le$CellID)

  le
})


# combine experiments
looms <- Reduce(cbind, looms)

# add rownames and colnames
rownames(looms) <- rowData(looms)$Accession
colnames(looms) <- looms$CellID
# export(looms, "data/intermediate/velocyto/ring_combined_nofilt.loom")


# Filter ------------------------------------------------------------------

# read filtered sce object
# ring_sce <- readRDS("data/processed/SingleCellExperiment/ring_batches_hardfilt.rds")
ring_sce <- readRDS("data/processed/SingleCellExperiment/ring_batches_strictfilt.rds")

# check intersection between these two objects
row_ids <- intersect(rownames(ring_sce), rownames(looms))
col_ids <- intersect(colnames(ring_sce), colnames(looms))

# filter to retain only those cells/genes
looms <- looms[row_ids, col_ids]
ring_sce <- ring_sce[row_ids, col_ids] # ensure genes/cells ordering

# add clustering information
looms$cluster_mnn_logvst <- ring_sce$cluster_mnn_logvst

# save all cells
export(looms, "data/intermediate/velocyto/ring_strictfilt.loom")

# export UMAP (to project velocities on)
x <- reducedDim(ring_sce, "UMAP30_MNN_logvst")
colnames(x) <- c("UMAP1", "UMAP2")
write.csv(x, "data/intermediate/velocyto/umap_for_scvelo.csv")


# # filter to retain only the initially determined 1000 HGVs
# export(looms[which(rowData(looms)$Accession %in% metadata(ring_sce)$hvgs), ],
#        "data/intermediate/velocyto/ring_combined_1000hvgs.loom")


# # Separate cell types -----------------------------------------------------
#
# # based on previous annotations, we separate different cell types to get individual lineages
# # cluster 9 is included in all sets (cycling cells, presumably initials)
#
# # clusters for each set
# annotated_clusters <- list(
#   cc_clusters = c(9, 7, 12),
#   ppp_clusters = c(9, 1, 4, 10),
#   se_clusters = c(9, 6, 5)
# )
#
# # subset and save each
# for(i in names(annotated_clusters)){
#   message("exporting: ", i)
#   export(looms[, which(looms$cluster_mnn_logvst %in% annotated_clusters[[i]])],
#          paste0("data/intermediate/velocyto/ring_", i, ".loom"))
# }
#
# # clusters for each set - without cycling cells
# annotated_clusters <- list(
#   cc_clusters = c(7, 12),
#   ppp_clusters = c(1, 4, 10),
#   se_clusters = c(6, 5)
# )
#
# # subset and save each
# for(i in names(annotated_clusters)){
#   message("exporting: ", i)
#   export(looms[, which(looms$cluster_mnn_logvst %in% annotated_clusters[[i]])],
#          paste0("data/intermediate/velocyto/ring_", i, "_nocycling.loom"))
# }
#
#
# # save set without outer layers
# export(looms[, which(looms$cluster_mnn_logvst != 11)],
#        paste0("data/intermediate/velocyto/ring_combined_noouter.loom"))
