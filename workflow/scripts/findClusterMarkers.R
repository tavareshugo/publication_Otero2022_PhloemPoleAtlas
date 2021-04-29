suppressPackageStartupMessages({
  library(data.table)
  library(scran)
})


# User Arguments ----------------------------------------------------------

# Fetch input and output file names
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  stop("Need to supply input file name")
} else if (!file.exists(args)) {
  stop("Input file does not exist: ", args)
} else {
  infile <- args[1]
}


# Read data ---------------------------------------------------------------
message("Reading ", infile)

sce <- readRDS(infile)

# retain nuclear genes only
sce <- sce[metadata(sce)$nucgenes, ]


# Find markers in 8/15 clusters -----------------------------------------
message("Finding markers...")

# run test for each marker
markers8in15 <- findMarkers(sce,
                            assay.type = "logvst",
                            groups = sce$cluster_mnn_logvst,
                            block = sce$Sample,
                            test = "wilcox", direction = "up",
                            pval.type = "some",
                            # initially ran this with 0.9 but that was a bit stringent
                            min.prop = 8/15,
                            # set fold-change
                            lfc = log(1.5))

# tidy up
markers8in15 <- lapply(markers8in15, function(i){
  as.data.table(i, keep.rownames = "id")
})
markers8in15 <- rbindlist(markers8in15, idcol = "cluster", fill = TRUE)
markers8in15[, cluster := factor(as.numeric(cluster))]


# Find markers in 8/15 clusters -----------------------------------------
message("Finding markers...")

# run test for each marker
markers12in15 <- findMarkers(sce,
                            assay.type = "logvst",
                            groups = sce$cluster_mnn_logvst,
                            block = sce$Sample,
                            test = "wilcox", direction = "up",
                            pval.type = "some",
                            # initially ran this with 0.9 but that was a bit stringent
                            min.prop = 12/15,
                            # set fold-change
                            lfc = log(1.5))

# tidy up
markers12in15 <- lapply(markers12in15, function(i){
  as.data.table(i, keep.rownames = "id")
})
markers12in15 <- rbindlist(markers12in15, idcol = "cluster", fill = TRUE)
markers12in15[, cluster := factor(as.numeric(cluster))]


# Add TF information ------------------------------------------------------

# TFs
tfs <- fread("data/external/transcription_factors/PlantTFDB_tidy.csv")

markers8in15[, is_tf := id %in% tfs$gene]
markers12in15[, is_tf := id %in% tfs$gene]


# Write Result ------------------------------------------------------------
message("Writing result as a CSV file...")

base_filename <- tools::file_path_sans_ext(basename(infile))
fwrite(markers8in15,
       paste0("data/processed/gene_sets/", base_filename, "_cluster_markers_8in15.csv"),
       sep = ",")
fwrite(markers12in15,
       paste0("data/processed/gene_sets/", base_filename, "_cluster_markers_12in15.csv"),
       sep = ",")
