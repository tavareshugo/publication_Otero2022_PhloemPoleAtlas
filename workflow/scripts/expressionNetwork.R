# Setup -------------------------------------------------------------------

if(!("bigSCale" %in% rownames(installed.packages()))){
  remotes::install_github("iaconogi/bigSCale2")
}

library(SingleCellExperiment)
library(bigSCale)

# sce object
sce <- readRDS("data/processed/SingleCellExperiment/ring_batches_hardfilt.rds")


# GRN ---------------------------------------------------------------------

grn <- compute.network(assay(sce, "counts"),
                gene.names = rownames(sce),
                clustering = "direct")


# Save output -------------------------------------------------------------

saveRDS(grn, "data/intermediate/expression_network/ring_hardfilt_network.rds")
