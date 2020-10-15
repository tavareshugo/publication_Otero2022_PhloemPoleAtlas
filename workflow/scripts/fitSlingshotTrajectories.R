suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(scran)
  library(scater)
  library(data.table)
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


# Re-cluster --------------------------------------------------------------

# remove outer layer clusters
sce <- sce[, which(sce$cluster_mnn_logvst != 11)]

# remove genes with too few counts
sce <- sce[which(rowSums(counts(sce)) > 5), ] # genes with total 5 counts
sce <- sce[which(rowSums(counts(sce) > 1) > 5), ] # genes in at least 5 cells

# Clustering - finer clusters
graph <- buildSNNGraph(sce, k = 50,
                       use.dimred = "MNN_logvst",
                       type = "jaccard")
sce$clusters_slingshot <- factor(igraph::cluster_louvain(graph)$membership)
rm(graph)

# visualise
sce %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = clusters_slingshot)) +
  geom_label(stat = "centroid",
             aes(group = clusters_slingshot, label = clusters_slingshot),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  theme_void() +
  ggthemes::scale_colour_tableau(palette = "Tableau 20") +
  ggthemes::scale_fill_tableau(palette = "Tableau 20")



# Fit trajectories --------------------------------------------------------
message("Fitting trajectories...")

# slingshot trajectories
sce <- slingshot(sce,
                 clusterLabels = sce$clusters_slingshot,
                 reducedDim = "DIFFMAP_MNN_logvst",
                 start.clus = 11,
                 end.clus = c(15, 6, 7))

# options to tinker with
# reweight = TRUE by default
# reassign = TRUE by default



# Fit GAM model to genes --------------------------------------------------

fit_gam <- function(Y, t){
  ngenes <- nrow(Y)
  gam_pvals <- vector("numeric", length = ngenes)
  gam_rsq <- vector("numeric", length = ngenes)
  gam_devexpl <- vector("numeric", length = ngenes)

  for(i in 1:nrow(Y)){
    fit <- gam(Y[i,] ~ s(t))
    gam_pvals[i] <- summary(fit)$s.table[4]
    gam_rsq[i] <- summary(fit)$r.sq
    gam_devexpl[i] <- summary(fit)$dev.expl
    message("Gene ", i, " of ", ngenes)
  }

  data.frame(gene = rownames(Y),
             pval = gam_pvals,
             rsq = gam_rsq,
             devexpl = gam_devexpl)
}

Y <- as.matrix(assay(sce, "logvst"))

result <- list(
  trajectory1 = fit_gam(Y, sce$slingPseudotime_1),
  trajectory3 = fit_gam(Y, sce$slingPseudotime_3),
  trajectory4 = fit_gam(Y, sce$slingPseudotime_4),
  trajectory5 = fit_gam(Y, sce$slingPseudotime_5)
)

result <- rbindlist(result, idcol = "trajectory")


# Write Result ------------------------------------------------------------
message("Writing result as a CSV file...")

fwrite(result, outfile, sep = ",")
