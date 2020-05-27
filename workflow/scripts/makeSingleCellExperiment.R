# Setup -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(biomaRt)
  library(data.table)
  library(DropletUtils)
  library(scater)
  library(scran)
  library(ggplot2)
  library(ggpointdensity)
})

# set seed for reproducible results
set.seed(1001)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))

# ggplot S3 method for BioConductor::DataFrame
ggplot.DFrame <- function(data, ...) ggplot(as.data.frame(data), ...)



# Capture user input ------------------------------------------------------

option_list = list(
  make_option(c("--cellranger"),
              action = "store",
              default = NA,
              type = 'character',
              help = "The directory with cellranger output."),
  make_option(c("--sample"),
              action = "store",
              default = NA,
              type = "character",
              help = "The sample name, which will be used to name output files."),
  make_option(c("--outdir"),
                action = "store",
                default = NA,
                type = 'character',
                help = "An output directory."),
  make_option("--fdr",
                action = "store",
                default = 0.001,
                type = 'double',
                help = "The FDR threshold used for identifying empty droplets."),
  make_option("--cores",
              action = "store",
              default = NA,
              type = 'integer',
              help = "Number of cores to use. Will use all available by default.")
  )
opt <- parse_args(OptionParser(option_list=option_list))

# opt$cellranger <- "data/intermediate/counts_cellranger/sAPL/"
# opt$outdir <- "data/processed/SingleCellExperiment/sAPL/"
# opt$fdr <- 0.001
# opt$sample <- basename(dirname(dirname(opt$cellranger)))
# opt$cores <- 3

# user input
if(is.na(opt$cellranger) | !dir.exists(opt$cellranger)){
  stop("Input directory does not exist.")
} else {
  opt$cellranger <- file.path(opt$cellranger, "outs/raw_feature_bc_matrix/")
}

if(is.na(opt$sample)){
  message("Inferring sample name automatically from directory name")
  opt$sample <- basename(dirname(dirname(opt$cellranger)))
}

if(is.na(opt$outdir)){
  stop("Please provide an output directory.")
} else if(!dir.exists(opt$outdir)){
  message("Output directory does not exist. Will be created.")
  dir.create(opt$outdir, recursive = TRUE)
}

if(opt$fdr < 0 | opt$fdr > 1){
  stop("FDR has to be between 0 and 1.")
}

if(is.na(opt$cores) | opt$cores > parallel::detectCores()){
  opt$cores <- BiocParallel::MulticoreParam(parallel::detectCores())
} else {
  opt$cores <- BiocParallel::MulticoreParam(opt$cores)
}

message("Running job with following options:")
print(opt)


# Read droplet data -------------------------------------------------------

# Read 10x data
sce <- read10xCounts(opt$cellranger, sample.names = opt$sample, col.names = TRUE)

# remove droplets with no data at all
sce <- sce[, which(colSums(counts(sce)) > 0)]

# remove genes with no data at all
# not sure I should do this, because later on want to combine datasets with MNN
# sce <- sce[which(rowSums(counts(sce)) > 0),]



# Add gene metadata -------------------------------------------------------

# Specify the plants_mart
m <- useMart("plants_mart", host="plants.ensembl.org")

# Define dataset
m <- useMart("plants_mart", host="plants.ensembl.org", dataset="athaliana_eg_gene")

# get all gene information
gene_meta <- getBM(attributes=c("ensembl_gene_id",
                                "tair_symbol",
                                "external_gene_name",
                                "chromosome_name",
                                "start_position", "end_position"),
                    mart=m)
setDT(gene_meta)

# Collapse alternative names for each ID
gene_meta <- gene_meta[ , .(alternative_name = paste(c(tair_symbol, external_gene_name),
                                                     collapse = "/")) ,
           by = .(ensembl_gene_id, chromosome_name, start_position, end_position)]
gene_meta[, alternative_name := ifelse(alternative_name == "/", NA, alternative_name)]
gene_meta$alternative_name <- gsub("^/", "", gene_meta$alternative_name)
gene_meta$alternative_name <- gsub("/$", "", gene_meta$alternative_name)

setnames(gene_meta,
         old = c("ensembl_gene_id", "chromosome_name"),
         new = c("ID", "chrom"))

# get genes with GO terms of interest
cell_cycle <- getBM(attributes=c("ensembl_gene_id"),
                    mart=m,
                    filters = "go", values = "GO:0007049")$ensembl_gene_id

mitosis <- getBM(attributes=c("ensembl_gene_id"),
                 mart=m,
                 filters = "go", values = "GO:0000278")$ensembl_gene_id

phloem_dev <- getBM(attributes=c("ensembl_gene_id"),
                    mart=m,
                    filters = "go", values = "GO:0010088")$ensembl_gene_id

# add to rowData
rowData(sce)$go_cell_cycle <- rowData(sce)$ID %in% cell_cycle
rowData(sce)$go_mitosis <- rowData(sce)$ID %in% mitosis
rowData(sce)$go_phloem_dev <- rowData(sce)$ID %in% phloem_dev

# add other gene meta information
gene_meta <- merge(gene_meta, rowData(sce), by = "ID", all.y = TRUE)
rownames(gene_meta) <- gene_meta$ID
rowData(sce) <- gene_meta[rownames(sce), ] # ensure correct order

if(!all(rownames(sce) == rowData(sce)$ID)) stop("rowData is corrupted.")

# add chromosome for fluorescent proteins
rowData(sce)$chrom[which(is.na(rowData(sce)$chrom))] <- "FP"


# Identify empty droplets --------------------------------------------------
message("Identifying empty droplets.")

# Run test - including the ambient cells in output to check uniform p-value distribution
empty <- emptyDrops(counts(sce), niters = 1e5, lower = 100, test.ambient = TRUE,
                    BPPARAM = opt$cores)

# store p-value distribution of "null" cells
null_pvals <- empty$PValue[empty$Total <= 100 & empty$Total > 0]

# re-correct FDR by ignoring the "null" cells
empty$PValue[empty$Total <= 100 & empty$Total > 0] <- NA
empty$FDR <- p.adjust(empty$PValue, method = "BH")

# check if any cells are limited by number of iterations and re-run if so
if(sum(empty$FDR > 0.001 & empty$Limited, na.rm = TRUE) > 0){
  nlim <- sum(empty$FDR > 0.001 & empty$Limited, na.rm = TRUE)
  warning("Empty droplet test had ", nlim,
          " p-values limited by the number of iterations of the MC algorithm.\n",
          "Will re-run using 5x more iterations (will take a while...)")
  empty <- emptyDrops(counts(sce), niters = 5e5, lower = 100,
                      BPPARAM = opt$cores)
}
if(sum(empty$FDR > 0.001 & empty$Limited, na.rm = TRUE) > 0){
  stop("Empty droplet test still had limited p-values! Check what's wrong?")
}


# Barcode ranks ----------------------------------------------------------
message("Getting barcode ranks.")

# Get barcode ranks
bcrank <- barcodeRanks(counts(sce))

# add information about being empty cell or not (for plotting later on)
bcrank$is_cell <- empty$FDR <= opt$fdr

# add information to colData
colData(sce)$bc_rank <-  bcrank$rank
colData(sce)$bc_total <-  bcrank$total
colData(sce)$bc_fitted <-  bcrank$fitted


# Quality control metrics -----------------------------------------------
message("Adding QC metrics.")

# remove empty droplets
sce <- sce[, which(empty$FDR <= opt$fdr)]

# quality control metrics for cells
sce <- addPerCellQC(sce,
                    subsets = list(mitochondria = grep("^ATM", rownames(sce)),
                                   chloroplast = grep("^ATC", rownames(sce)),
                                   cellcycle = which(rowData(sce)$go_cell_cycle),
                                   mitosis = which(rowData(sce)$go_mitosis),
                                   phloem = which(rowData(sce)$go_phloem_dev)),
                    BPPARAM = opt$cores)

# quality control metrics for genes
sce <- addPerFeatureQC(sce,
                       BPPARAM = opt$cores)



# Normalise ---------------------------------------------------------------
message("Log-normalising counts.")

# Here we just do a simple log-normalisation
sce <- computeSumFactors(sce,
                         clusters = quickCluster(sce, BPPARAM = opt$cores),
                         min.mean = 0.1)
sce <- logNormCounts(sce)

# Add mean and variance per gene to rowData
rowData(sce)$mean_logcounts <- rowMeans(logcounts(sce))
rowData(sce)$var_logcounts <- rowVars(as.matrix(logcounts(sce)))

# variance-stabilising transform (need to filter on genes with counts for at least one cell)
# m <- as.matrix(counts(sce))
# colnames(m) <- colData(sce)$Barcode
# trans <- sctransform::vst(m, min_cells = 1)



# dimensionality reduction ------------------------------------------------
message("Running PCA, UMAP and TSNE.")

# add dimensionality reduction slots - for exploratory analysis
sce <- runPCA(sce)
sce <- runUMAP(sce)
sce <- runTSNE(sce)
# sce <- runDiffusionMap(sce)



# Save SingleCellExperiment object ----------------------------------------
message("Saving SingleCellExperiment object.")

saveRDS(sce, file.path(opt$outdir, paste0(opt$sample, "_sce.rds")))



# Make QC graphs ----------------------------------------------------------
message("Making output graphs.")

# open PDF graphical device
pdf(file.path(opt$outdir, paste0(opt$sample, "_qc.pdf")))

# Distribution of "null" p-values from emptyDrops()
qplot(null_pvals, geom = "histogram", binwidth = 0.1) +
  xlim(0, 1) +
  labs(x = "P-value",
       subtitle = "Distribution of p-values for 'null' cells used by emptyDrops()")

# barcode rank graph
ggplot(bcrank[!duplicated(bcrank$rank), ],
       aes(x = rank, y = total)) +
  geom_point() +
  geom_hline(aes(yintercept = metadata(bcrank)$inflection, colour ="inflection")) +
  geom_hline(aes(yintercept = metadata(bcrank)$knee, colour = "knee")) +
  geom_hline(aes(yintercept = 100, colour = "limit for emptyDrops test")) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Rank", y = "Total UMI count", colour = "",
       subtitle = paste("Number cells:", ncol(sce)))

# total vs detected
ggplot(colData(sce), aes(total, detected)) +
  geom_point() +
  labs(x = "Total UMIs", y = "Detected genes")


# Distributions
ggplot(colData(sce),
       aes(total)) +
  geom_density(fill = "lightgrey") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  labs(x = "Total UMIs per cell")

ggplot(colData(sce),
       aes(detected)) +
  geom_density(fill = "lightgrey") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  labs(x = "# detected genes per cell")

ggplot(rowData(sce), aes(mean + 0.001)) +
  geom_density(aes(colour = chrom)) +
  scale_x_log10() +
  scale_colour_brewer(palette = "Dark2") +
  annotation_logticks(sides = "b") +
  labs(x = "Mean expression per cell (+0.001)")


# total UMI vs % UMIs in MT genes
ggplot(colData(sce), aes(sum, subsets_mitochondria_percent)) +
  geom_pointdensity(show.legend = FALSE) +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  scale_colour_viridis_c() +
  geom_rug(alpha = 0.1) +
  labs(x = "UMI counts (total)", y = "UMI counts (% from mitochondrial genes)",
       caption = "Each point represents a cell. Colour represents density of points in the graph.")

# total counts vs chloroplast
ggplot(colData(sce), aes(sum, subsets_chloroplast_percent)) +
  geom_pointdensity(show.legend = FALSE) +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  scale_colour_viridis_c() +
  geom_rug(alpha = 0.1) +
  # annotate(geom = "rect", xmin = 1e3, xmax = Inf, ymin = -Inf, ymax = 10,
  #          colour = "red3", fill = NA) +
  labs(x = "UMI counts (total)", y = "UMI counts (% from chloroplast genes)",
       caption = "Each point represents a cell. Colour represents density of points in the graph.")

# chloroplast vs mitocondria pct
ggplot(colData(sce),
       aes(subsets_mitochondria_percent, subsets_chloroplast_percent)) +
  geom_pointdensity(show.legend = FALSE) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  scale_colour_viridis_c() +
  geom_rug(alpha = 0.1) +
  # annotate(geom = "rect", xmin = 1e3, xmax = Inf, ymin = -Inf, ymax = 10,
  #          colour = "red3", fill = NA) +
  labs(x = "UMI counts (% from mitochondrial genes)",
       y = "UMI counts (% from chloroplast genes)",
       caption = "Each point represents a cell. Colour represents density of points in the graph.")


# Make a graph to illustrate log-FC between lost and kept cells
mito_rich <- colData(sce)$subsets_mitochondria_percent > 10
lost <- calculateAverage(counts(sce)[,which(mito_rich), drop = FALSE])
kept <- calculateAverage(counts(sce)[,which(!mito_rich), drop = FALSE])

logged <- edgeR::cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)
qplot(abundance, logFC, geom = "point") +
  geom_hline(yintercept = 0, colour = "dodgerblue") +
  labs(x = "Mean CPM", y = "Log-FC (lost/kept)",
       caption = "Comparing gene expression in group of cells with >10% vs <10% MT UMIs")


# mean-variance plot - log-transformed counts
ggplot(rowData(sce),
       aes(mean_logcounts, var_logcounts)) +
  geom_point() +
  geom_abline() +
  facet_wrap(~ chrom) +
  scale_colour_manual(values = c("black", "dodgerblue")) +
  labs(x = "Mean log-counts", y = "Variance log-counts",
       colour = "Mitochondrial genes",
       caption = "mean-variance relationship of logCounts.")


# Dimensionality reduction graphs
ggplot(cbind(as.data.frame(reducedDim(sce, "PCA")), colData(sce)),
       aes(PC1, PC2)) +
  geom_point(aes(colour = sum, shape = subsets_mitochondria_percent > 10)) +
  scale_colour_viridis_c(trans = "log10") +
  scale_shape_manual(values = c(4, 16)) +
  labs(shape = ">10% MT UMIs",
       title = "PCA", caption = "on logCounts (no cell filtering)")

ggplot(cbind(as.data.frame(reducedDim(sce, "TSNE")), colData(sce)),
       aes(V1, V2)) +
  geom_point(aes(colour = sum, shape = subsets_mitochondria_percent > 10)) +
  scale_colour_viridis_c(trans = "log10") +
  scale_shape_manual(values = c(4, 16)) +
  labs(shape = ">10% MT UMIs",
       title = "TSNE", caption = "on logCounts (no cell filtering)")

ggplot(cbind(as.data.frame(reducedDim(sce, "UMAP")), colData(sce)),
       aes(V1, V2)) +
  geom_point(aes(colour = sum, shape = subsets_mitochondria_percent > 10)) +
  scale_colour_viridis_c(trans = "log10") +
  scale_shape_manual(values = c(4, 16)) +
  labs(shape = ">10% MT UMIs",
       title = "UMAP", caption = "on logCounts (no cell filtering)")


# close graphical device
dev.off()

