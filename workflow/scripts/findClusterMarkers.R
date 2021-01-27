suppressPackageStartupMessages({
  library(data.table)
  library(scran)
})


# User Arguments ----------------------------------------------------------

# Fetch input and output file names
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Need to supply input and output file names")
} else if (!file.exists(args[1])) {
  stop("Input file does not exist: ", args[1])
} else {
  infile <- args[1]
  outfile <- args[2]
}


# Read data ---------------------------------------------------------------
message("Reading ", infile)

sce <- readRDS(infile)

# retain nuclear genes only
sce <- sce[metadata(sce)$nucgenes, ]


# Find markers ------------------------------------------------------------
message("Finding markers...")

# run test for each marker
markers <- findMarkers(sce,
                       assay.type = "logvst",
                       groups = sce$cluster_mnn_logvst,
                       block = sce$Sample,
                       test = "wilcox", direction = "up",
                       pval.type = "some",
                       # initially ran this with 0.9 but that was a bit stringent
                       min.prop = 0.5,
                       # set fold-change
                       lfc = log(1.5))

# tidy up
markers <- lapply(markers, function(i){
  as.data.table(i, keep.rownames = "id")
})
markers <- rbindlist(markers, idcol = "cluster", fill = TRUE)
markers[, cluster := factor(as.numeric(cluster))]



# Add TF information ------------------------------------------------------

# TFs
tfs <- fread("data/external/transcription_factors/PlantTFDB_tidy.csv")

markers[, is_tf := id %in% tfs$gene]


# Write Result ------------------------------------------------------------
message("Writing result as a CSV file...")

fwrite(markers, outfile, sep = ",")
