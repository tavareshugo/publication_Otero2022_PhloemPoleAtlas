if(!("tradeSeq" %in% rownames(installed.packages()))){
  BiocManager::install("tradeSeq")
}

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

# # remove genes with too few counts
# sce <- sce[which(rowSums(counts(sce)) > 5), ] # at least total 5 counts
# sce <- sce[which(rowSums(counts(sce) > 1) > 5), ] # > 1 count in at least 5 cells


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
tradeseq <- fitGAM(counts = counts(sce[, rownames(pseudotime)]),
                   pseudotime = pseudotime,
                   cellWeights = cell_weights,
                   nknots = 7,
                   verbose = TRUE)

# sce <- fitGAM(counts = sce,
#               nknots = 7,
#               verbose = TRUE)


# Trajectory tests --------------------------------------------------------

# general change in time
general_test <- associationTest(tradeseq, global = TRUE, lineages = TRUE)
general_test <- data.table::as.data.table(general_test, keep.rownames = "gene")

# difference between start and end (per lineage)
start_end_lineage <- startVsEndTest(tradeseq, global = TRUE, lineages = TRUE)
start_end_lineage <- data.table::as.data.table(start_end_lineage, keep.rownames = "gene")

# pairwise comparisons for the expression at the end of the trajectory
end_pairwise <- diffEndTest(tradeseq, global = TRUE, pairwise = TRUE)
end_pairwise <- data.table::as.data.table(end_pairwise, keep.rownames = "gene")
#end_res <- diffEndTest(tradeseq)

# test if pattern of expression differs between trajectories
pattern_pairwise <- patternTest(tradeseq, global = TRUE, pairwise = TRUE)
pattern_pairwise <- data.table::as.data.table(pattern_pairwise, keep.rownames = "gene")


# Fit GAM model to genes --------------------------------------------------
# I also fit a model "manually"
message("Fitting GAM models...")

# this is similar to what is implemented in tradeSeq
# but we fit to the normalised logcounts for visualisation purposes
fit_gam <- function(Y, t){
  ngenes <- nrow(Y)
  gam_pvals <- vector("numeric", length = ngenes)
  gam_rsq <- vector("numeric", length = ngenes)
  gam_devexpl <- vector("numeric", length = ngenes)
  gam_aic <- vector("numeric", length = ngenes)
  pred <- vector("list", length = ngenes)

  for(i in 1:nrow(Y)){
    # fit <- gam(Y[i,] ~ s(t), family = "nb") # on raw counts
    fit <- gam(Y[i,] ~ s(t))
    gam_pvals[i] <- summary(fit)$s.table[4]
    gam_rsq[i] <- summary(fit)$r.sq
    gam_devexpl[i] <- summary(fit)$dev.expl
    gam_aic[i] <- AIC(fit)
    newd <- data.frame(t = seq(min(t, na.rm = TRUE), max(t, na.rm = TRUE), length.out = 100))
    pred[[i]] <- data.frame(t = newd$t,
                            .pred = predict(fit, newdata = newd))
    message("Gene ", i, " of ", ngenes)
  }

  data.table(gene = rownames(Y),
             pval = gam_pvals,
             rsq = gam_rsq,
             devexpl = gam_devexpl,
             aic = gam_aic,
             pred = pred)
}

# Fit the model
Y <- as.matrix(logcounts(sce))
result <- list(
  trajectory1 = fit_gam(Y, sce$slingPseudotime_1),
  trajectory2 = fit_gam(Y, sce$slingPseudotime_2),
  trajectory3 = fit_gam(Y, sce$slingPseudotime_3)
)
result <- rbindlist(result, idcol = "trajectory")


# Write Result ------------------------------------------------------------
message("Writing result...")

outfile <- tools::file_path_sans_ext(basename(infile))

saveRDS(
  sce,
  file = paste0("data/processed/trajectories/", outfile, "_slingshot.rds")
)
saveRDS(
  tradeseq,
  file = paste0("data/processed/trajectories/", outfile, "_tradeseq.rds")
)
saveRDS(
  result,
  file = paste0("data/processed/trajectories/", outfile, "_gams.rds")
)
fwrite(general_test, "data/processed/trajectories/tradeseq_associationTest.csv", sep = ",")
fwrite(start_end_lineage, "data/processed/trajectories/tradeseq_startVsEndTest.csv", sep = ",")
fwrite(end_pairwise, "data/processed/trajectories/tradeseq_diffEndTest.csv", sep = ",")
fwrite(pattern_pairwise, "data/processed/trajectories/tradeseq_patternTest.csv", sep = ",")
# fwrite(result, "data/processed/trajectories/slingshot_trajectory_gam.csv", sep = ",")
# saveRDS(sce, "data/processed/trajectories/ring_hardfilt_slingshot.rds")


# # Test GAM models --------------------------------------------------
#
# # this is essentially identical to what is implemented in tradeSeq
# fit_gam <- function(Y, t, k = -1){
#   ngenes <- nrow(Y)
#   gam_pvals <- vector("numeric", length = ngenes)
#   gam_rsq <- vector("numeric", length = ngenes)
#   gam_devexpl <- vector("numeric", length = ngenes)
#   gam_aic <- vector("numeric", length = ngenes)
#   pred <- vector("list", length = ngenes)
#
#   for(i in 1:nrow(Y)){
#     fit <- gam(Y[i,] ~ s(t, bs = "cr", k = k), family = "nb")
#     gam_pvals[i] <- summary(fit)$s.table[4]
#     gam_rsq[i] <- summary(fit)$r.sq
#     gam_devexpl[i] <- summary(fit)$dev.expl
#     gam_aic[i] <- AIC(fit)
#     newd <- data.frame(t = seq(min(t, na.rm = TRUE), max(t, na.rm = TRUE), length.out = 100))
#     pred[[i]] <- data.frame(t = newd$t,
#                             y = predict(fit, newdata = newd, type = "link")/log(10))
#     message("Gene ", i, " of ", ngenes)
#   }
#
#   data.table(gene = rownames(Y),
#              pval = gam_pvals,
#              rsq = gam_rsq,
#              devexpl = gam_devexpl,
#              aic = gam_aic,
#              pred = pred)
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
