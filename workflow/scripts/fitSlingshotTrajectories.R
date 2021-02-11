suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(scran)
  library(scater)
  library(data.table)
  library(mgcv)
  library(tradeSeq)
})


# User Arguments ----------------------------------------------------------

# Fetch input and output file names
args <- commandArgs(trailingOnly=TRUE)

if (!file.exists(args[1])) {
  stop("Input file does not exist: ", args[1])
}
infile <- args[1]

# if (length(args) < 2) {
#   stop("Need to supply input and output file names")
# } else if (!file.exists(args[1])) {
#   stop("Input file does not exist: ", args[1])
# } else {
#   infile <- args[1]
# }


# Read data ---------------------------------------------------------------
message("Reading ", infile)

sce <- readRDS(infile)

# remove genes with too few counts
sce <- sce[which(rowSums(counts(sce)) > 5), ] # at least total 5 counts
sce <- sce[which(rowSums(counts(sce) > 1) > 5), ] # > 1 count in at least 5 cells


# Fit trajectories --------------------------------------------------------
message("Fitting trajectories...")

# slingshot trajectories
# start cluster are based on prior annotation
sce <- slingshot(sce,
                 clusterLabels = sce$cluster_mnn_logvst,
                 reducedDim = "DIFFMAP_MNN_logvst",
                 start.clus = c(13))

# options to tinker with
# reweight = TRUE by default
# reassign = TRUE by default

# Fit GAM models with tradeSeq --------------------------------------------
message("Fitting tradeSeq models...")

# get trajectories of interest (biologically meaningful)
cell_weights <- slingCurveWeights(SlingshotDataSet(sce))[, c(1, 2, 3)]
cell_weights <- cell_weights[which(rowSums(cell_weights) > 0), ]
pseudotime <- slingPseudotime(SlingshotDataSet(sce), na = FALSE)[, c(1, 2, 3)]
pseudotime <- pseudotime[rownames(cell_weights), ]

# Fit GAM models
set.seed(7)
tradeseq_sce <- fitGAM(counts = counts(sce[, rownames(pseudotime)]),
                       pseudotime = pseudotime,
                       cellWeights = cell_weights,
                       nknots = 7,
                       verbose = TRUE,
                       parallel = TRUE,
                       BPARAM = BiocParallel::MulticoreParam(parallel::detectCores()))

if(any(!rowData(sce)$tradeSeq$converged)){
  warning(sum(!rowData(sce)$tradeSeq$converged), " models didn't converge.")
}


# # Fit GAM model to genes --------------------------------------------------
#
# # this is essentially identical to what is implemented in tradeSeq
# fit_gam <- function(Y, t, k = -1){
#   ngenes <- nrow(Y)
#   gam_pvals <- vector("numeric", length = ngenes)
#   gam_rsq <- vector("numeric", length = ngenes)
#   gam_devexpl <- vector("numeric", length = ngenes)
#   gam_aic <- vector("numeric", length = ngenes)
#
#   for(i in 1:nrow(Y)){
#     fit <- gam(Y[i,] ~ s(t, bs = "cr", k = k), family = "nb")
#     gam_pvals[i] <- summary(fit)$s.table[4]
#     gam_rsq[i] <- summary(fit)$r.sq
#     gam_devexpl[i] <- summary(fit)$dev.expl
#     gam_aic[i] <- AIC(fit)
#     message("Gene ", i, " of ", ngenes)
#   }
#
#   data.frame(gene = rownames(Y),
#              pval = gam_pvals,
#              rsq = gam_rsq,
#              devexpl = gam_devexpl,
#              aic = gam_aic)
# }
#
# # testing optimal K - based on tradeSeq approach
# Y <- as.matrix(counts(sce)[sample(metadata(sce)$hvgs_vst, 100), ])
# result <- lapply(c(-1, seq(3, 15, by = 2)), function(i){
#   message("k = ", i)
#   out <- list(
#     trajectory1 = fit_gam(Y, sce$slingPseudotime_1, i),
#     trajectory2 = fit_gam(Y, sce$slingPseudotime_2, i),
#     trajectory3 = fit_gam(Y, sce$slingPseudotime_3, i)
#   )
#   out <- rbindlist(out, idcol = "trajectory")
#   out$k <- i
#   out
# })
# result <- rbindlist(result)
#
# result[, pct := aic / aic[k == min(k)], by = .(gene, trajectory)]
# ggplot(result,
#        aes(factor(k), pct)) +
#   geom_line(aes(group = gene)) +
#   geom_point(stat = "summary", fun = "median", colour = "red", size = 5) +
#   facet_wrap(~ trajectory, scales = "free")
#
# # # Fit the model
# # Y <- as.matrix(counts(sce))
# # result <- list(
# #   trajectory1 = fit_gam(Y, sce$slingPseudotime_1),
# #   trajectory2 = fit_gam(Y, sce$slingPseudotime_2),
# #   trajectory3 = fit_gam(Y, sce$slingPseudotime_3)
# # )
# # result <- rbindlist(result, idcol = "trajectory")


# Write Result ------------------------------------------------------------
message("Writing result...")

outfile <- tools::file_path_sans_ext(basename(infile))

saveRDS(
  sce,
  file = paste0("data/processed/trajectories/", outfile, "_slingshot.rds")
)
saveRDS(
  tradeseq_sce,
  file = paste0("data/processed/trajectories/", outfile, "_slingshot_tradeseq.rds")
)

# fwrite(result, "data/processed/trajectories/slingshot_trajectory_gam.csv", sep = ",")
# saveRDS(sce, "data/processed/trajectories/ring_hardfilt_slingshot.rds")
