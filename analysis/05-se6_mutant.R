library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(tidyverse)
library(readxl)
library(ggpointdensity)
library(patchwork)

# set ggplot2 theme
theme_set(theme_classic() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")

# min-max scaling
min_max <- function(x, na.rm = FALSE){
  (x - min(x, na.rm = na.rm))/(max(x, na.rm = na.rm) - min(x, na.rm = na.rm))
}

# tidy predictions from tradeseq
predict_tradeseq <- function(x, gene, scale = FALSE){
  gene <- gene[gene %in% rownames(x)]
  pred <- predictSmooth(x, gene = gene)
  if(scale) pred <- pred %>% t() %>% scale() %>% t()

  pred %>%
    as_tibble(rownames = "gene") %>%
    pivot_longer(-gene, values_to = ".pred") %>%
    separate(name, c("trajectory", "t"), sep = "_", convert = TRUE)
}


# Read data ---------------------------------------------------------------

# SCE object after slingshot
sling <- readRDS("data/processed/trajectories/ring_batches_strictfilt_slingshot.rds")

# results from tradeSeq
tradeseq <- readRDS("data/processed/trajectories/ring_batches_strictfilt_tradeseq.rds")

# results from tradeSeqtests
association <- read_csv("data/processed/trajectories/tradeseq_associationTest.csv")
startend <- read_csv("data/processed/trajectories/tradeseq_startVsEndTest.csv")
diffend <- read_csv("data/processed/trajectories/tradeseq_diffEndTest.csv")
pattern <- read_csv("data/processed/trajectories/tradeseq_patternTest.csv")

# results from fitting GAM models
gamfits <- readRDS("data/processed/trajectories/ring_batches_strictfilt_gams.rds") %>%
  as_tibble()

# all tradeseq test results together
gene_tests <- bind_rows(
  pattern %>% select(gene, waldStat, pvalue) %>% mutate(test = "pattern"),
  diffend %>% select(gene, waldStat, pvalue) %>% mutate(test = "diffend"),
  startend %>% select(gene, waldStat, pvalue) %>% mutate(test = "startend")
)

# genes with ring expression (from microscopy)
ring_genes <- read_csv("data/raw/ring_genes.csv") %>%
  select(name = gene_name, gene = gene_id, domain) %>%
  mutate(name = ifelse(is.na(name), gene, name))

# sixtuple mutant results from previous work
se6 <- read_excel("data/raw/SE6vsCol_DEG.xlsx",
                  sheet = "SE6vsCol.DEG",
                  range = cell_cols("A:F")) %>%
  rename(gene = Gene_id, se6 = readcount_SE6, wt = readcount_Col)

# fold-change distribution (these were filtered for padj < 0.05 already)
se6 %>%
  ggplot(aes(log2FoldChange)) +
  geom_density(fill = "steelblue") +
  geom_vline(xintercept = c(-1, 1), lty = 2)


# Get fit trajectories ----------------------------------------------------

# score based on association with any trajectory
gois_association <- association %>%
  filter(waldStat >= quantile(waldStat, 0.95, na.rm = TRUE)) %>%
  left_join(gamfits) %>% filter(rsq > 0.6) %>%
  pull(gene) %>% unique()

gois <- se6 %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
  filter(gene %in% rownames(tradeseq)) %>%
  pull(gene)

# no particular skew in test statistic distributions compared to overall
association %>%
  ggplot(aes(waldStat)) +
  geom_density(aes(colour = gene %in% gois)) +
  scale_x_log10()

# intersect these two
gois <- gois[gois %in% gois_association]

set.seed(20210223)
gene_clusters <- predictSmooth(tradeseq, gois) %>%
  t() %>% scale() %>% t() %>%
  kmeans(centers = 10, nstart = 100, iter.max = 100) %>%
  `$`(cluster) %>%
  enframe("gene", "cluster") %>%
  left_join(predict_tradeseq(tradeseq, gois, scale = TRUE)) %>%
  group_by(cluster) %>%
  # mutate(cluster = paste0(LETTERS[cluster], "\n(", n_distinct(gene), ")"))
  mutate(ngenes = n_distinct(gene)) %>%
  ungroup() %>%
  mutate(cluster = LETTERS[cluster])

# pdf("documents/pdf for figures/trajectory_gene_clusters.pdf",
#     width = 3, height = 4)
gene_clusters %>%
  # full_join(cluster_order, by = "cluster") %>%
  # mutate(cluster = paste0(LETTERS[10-order], "\n(", ngenes, ")")) %>%
  mutate(trajectory = case_when(trajectory == "lineage1" ~ "PPP",
                                trajectory == "lineage2" ~ "CC",
                                trajectory == "lineage3" ~ "PSE",
                                TRUE ~ NA_character_)) %>%
  # mutate(cluster = fct_reorder(cluster, order)) %>%
  ggplot(aes(t, .pred)) +
  geom_hline(yintercept = 0, colour = "steelblue") +
  geom_line(aes(group = interaction(gene, trajectory)),
            alpha = 0.2) +
  geom_line(stat = "summary", fun = median, colour = "brown", size = 1) +
  facet_grid(cluster ~ trajectory, scales = "free_y") +
  theme_void(base_size = 12)
# dev.off()



##### deprecated #####

# Get LFC between clusters ----------------------------------------------

markers <- scran::findMarkers(ring,
                              assay.type = "logvst",
                              groups = ring$cluster_mnn_logvst,
                              block = ring$Sample,
                              test = "t", direction = "up",
                              pval.type = "some",
                              min.prop = 0.5,
                              # set fold-change
                              lfc = log(1.5))

# tidy up
markers <- markers %>%
  lapply(function(i) as_tibble(i, rownames = "id")) %>%
  bind_rows(.id = "cluster") %>%
  mutate(cluster = factor(cluster))

comparisons <- combn(1:max(as.numeric(markers$cluster)), 2)
comparisons <- paste(comparisons[1, ], comparisons[2, ], sep = "_")

markers %>%
  select(cluster, id, starts_with("logFC.")) %>%
  pivot_longer(starts_with("logFC.")) %>%
  mutate(name = str_remove(name, "logFC.")) %>%
  mutate(comparison = paste(cluster, name, sep = "_")) %>%
  select(id, comparison, value) %>%
  drop_na(value) %>%
  filter(comparison %in% comparisons) %>%
  ggplot(aes(value)) +
  geom_density(aes(group = comparison)) +
  coord_cartesian(xlim = c(-0.5, 0.5))

# make a matrix with genes having at least one comparison with LFC > 1
lfc <- markers %>%
  # retain genes below FDR threshold
  group_by(id) %>%
  filter(any(FDR < 0.01)) %>%
  ungroup() %>%
  # reshape and tidy
  select(cluster, id, starts_with("logFC.")) %>%
  pivot_longer(starts_with("logFC.")) %>%
  mutate(name = str_remove(name, "logFC.")) %>%
  mutate(comparison = paste(cluster, name, sep = "_")) %>%
  select(id, comparison, value) %>%
  drop_na(value) %>%
  filter(comparison %in% comparisons) %>%
  # retain genes with LFC > 1
  group_by(id) %>%
  filter(any(abs(value) > 1)) %>%
  ungroup() %>%
  # turn into matrix
  pivot_wider(names_from = comparison) %>%
  column_to_rownames("id") %>%
  as.matrix()


# Up-regulated LFC heatmap --------------------------------------------

# top DEGs
up_se6 <- se6 %>%
  filter(padj < 0.01 & log2FoldChange >= 1) %>%
  pull(Gene_id)

# make heatmap matrix
heatmap_mat <- lfc[rownames(lfc) %in% up_se6, ]
gene_clust <- heatmap_mat %>%
  dist() %>%
  hclust() %>%
  cutree(7) %>%
  factor() %>%
  enframe(value = "cluster") %>%
  column_to_rownames("name")
cell_clust <- heatmap_mat %>%
  t() %>%
  dist() %>%
  hclust() %>%
  cutree(7) %>%
  factor() %>%
  enframe(value = "cluster") %>%
  column_to_rownames("name")
pheatmap(heatmap_mat[order(gene_clust$cluster), order(cell_clust$cluster)],
         cutree_rows = 7,
         cutree_cols = 7,
         annotation_col = cell_clust,
         annotation_row = gene_clust,
         main = "LFC > 1 in SE6",
         height = 10, width = 7)

interesting_genes <- rownames(gene_clust)[gene_clust$cluster %in% c(3, 4, 7)]
ring %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = interesting_genes, melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  # geom_label(stat = "centroid",
  #            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
  #            alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression") +
  theme_void() +
  coord_equal()


# Down-regulated LFC heatmap --------------------------------------------

# top DEGs
down_se6 <- se6 %>%
  filter(padj < 0.01 & log2FoldChange <= -1) %>%
  pull(Gene_id)

# make heatmap matrix
heatmap_mat <- lfc[rownames(lfc) %in% down_se6, ]
gene_clust <- heatmap_mat %>%
  dist() %>%
  hclust() %>%
  cutree(7) %>%
  factor() %>%
  enframe(value = "cluster") %>%
  column_to_rownames("name")
cell_clust <- heatmap_mat %>%
  t() %>%
  dist() %>%
  hclust() %>%
  cutree(7) %>%
  factor() %>%
  enframe(value = "cluster") %>%
  column_to_rownames("name")
pheatmap(heatmap_mat[order(gene_clust$cluster), order(cell_clust$cluster)],
         cutree_rows = 7,
         cutree_cols = 7,
         annotation_col = cell_clust,
         annotation_row = gene_clust,
         main = "LFC < -1 in SE6",
         height = 10, width = 7)

interesting_genes <- rownames(gene_clust)[gene_clust$cluster == 5]
ring %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = interesting_genes, melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  # geom_label(stat = "centroid",
  #            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
  #            alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression") +
  theme_void() +
  coord_equal()


# Up-regulated correlation with reference ---------------------------

# curated list of genes with known ring pattern
curated_ring <- c("AT3G16330", "AT3G21770", "AT2G02230",
                  "AT1G52140", "AT1G29160")

# this was produced in 04-cluster_markers
cor_reference <- read_csv("data/processed/gene_sets/correlation_with_reference_genes.csv")

# this is how the curated genes themselves are correlated with each other
cor_reference %>%
  filter(id %in% curated_ring) %>%
  ggplot(aes(id, reference_gene)) +
  geom_tile(aes(fill = correlation))

# get a subset of potentially interesting genes
interesting_genes <- cor_reference %>%
  filter(id != reference_gene & !id %in% curated_ring) %>%
  filter(id %in% up_se6) %>%
  group_by(id) %>%
  summarise(rank_min = min(rank)) %>%
  slice_min(rank_min, n = 15) %>%
  pull(id) %>% unique()

cor_reference %>%
  filter(id != reference_gene & !id %in% curated_ring) %>%
  group_by(id) %>%
  summarise(rank_sum = sum(rank),
            rank_mean = mean(rank),
            cor_mean = mean(correlation),
            cor_prod = prod(correlation)) %>%
  ggplot(aes(cor_prod, cor_mean)) +
  geom_point(aes(colour = id %in% interesting_genes)) +
  facet_wrap(~ id %in% up_se6)

ring %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = c(curated_ring, interesting_genes), melted = TRUE,
                exprs_values = "logcounts") %>%
  arrange(!is.na(expr)) %>%
  ggplot(aes(cluster_mnn_logvst, expr)) +
  geom_violin(scale = "width", fill = "steelblue") +
  facet_wrap(~ id, scales = "free") +
  labs(x = "Cluster", y = "logcounts",
       title = "Log2FC(SE6) > 1") +
  theme_classic()

ring %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = c(curated_ring, interesting_genes), melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  # geom_label(stat = "centroid",
  #            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
  #            alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression",
       title = "Log2FC(SE6) > 1") +
  theme_void() +
  coord_equal()


# Down-regulated correlation with reference ---------------------------

# get interesting genes
interesting_genes <- cor_reference %>%
  filter(id != reference_gene & !id %in% curated_ring) %>%
  filter(id %in% down_se6) %>%
  group_by(id) %>%
  summarise(rank_min = min(rank)) %>%
  slice_min(rank_min, n = 15, with_ties = FALSE) %>%
  pull(id) %>% unique()

cor_reference %>%
  filter(id != reference_gene & !id %in% curated_ring) %>%
  group_by(id) %>%
  summarise(rank_sum = sum(rank),
            rank_mean = mean(rank),
            cor_mean = mean(correlation),
            cor_prod = prod(correlation)) %>%
  ggplot(aes(cor_prod, cor_mean)) +
  geom_point(aes(colour = id %in% interesting_genes)) +
  facet_wrap(~ id %in% down_se6)

ring %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = c(curated_ring, interesting_genes), melted = TRUE,
                exprs_values = "logcounts") %>%
  arrange(!is.na(expr)) %>%
  ggplot(aes(cluster_mnn_logvst, expr)) +
  geom_violin(scale = "width", fill = "steelblue") +
  facet_wrap(~ id, scales = "free") +
  labs(x = "Cluster", y = "logcounts",
       title = "Log2FC(SE6) < -1") +
  theme_classic()

ring %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = c(curated_ring, interesting_genes), melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  # geom_label(stat = "centroid",
  #            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
  #            alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression",
       title = "Log2FC(SE6) < -1") +
  theme_void() +
  coord_equal()


# Up-regulated intersect ---------------------------------------------

# top DEGs
up_se6 <- se6 %>%
  filter(padj < 0.01 & log2FoldChange >= 1) %>%
  pull(Gene_id)

# top DEGs in single-cell data
degs_sc <- markers %>%
  filter(FDR < 0.01) %>%
  pivot_longer(starts_with("logFC.")) %>%
  filter(!cluster %in% c(3, 8, 11) & !name %in% paste0("logFC.", c(3, 8, 11))) %>%
  filter(abs(value) > 1) %>%
  distinct(id) %>%
  pull(id)

# intersect size
list(SE6 = up_se6, SC = degs_sc) %>%
  fromList() %>%
  upset()

up_intersect <- intersect(up_se6, degs_sc)

heatmap_mat <- ring %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = up_intersect, melted = TRUE,
                exprs_values = "logcounts", zeros_to_na = FALSE) %>%
  filter(!cluster_mnn_logvst %in% c(3, 8, 11)) %>%
  group_by(id, cluster_mnn_logvst) %>%
  summarise(expr = mean(expr)) %>%
  pivot_wider(names_from = "cluster_mnn_logvst",
              values_from = expr) %>%
  column_to_rownames("id") %>%
  as.matrix() %>%
  # scale
  t() %>% scale() %>% t()

gene_clust <- heatmap_mat %>%
  dist() %>%
  hclust() %>%
  cutree(7) %>%
  factor() %>%
  enframe(value = "cluster") %>%
  column_to_rownames("name")
pheatmap(heatmap_mat[order(gene_clust$cluster), ],
         cutree_rows = 7,
         annotation_row = gene_clust,
         main = "LFC > 1 in SE6",
         height = 10, width = 7)



# Down-regulated intersect ---------------------------------------------

# top DEGs
down_se6 <- se6 %>%
  filter(padj < 0.01 & log2FoldChange <= -1) %>%
  pull(Gene_id)

# top DEGs in single-cell data
degs_sc <- markers %>%
  filter(FDR < 0.01) %>%
  pivot_longer(starts_with("logFC.")) %>%
  filter(!cluster %in% c(3, 8, 11) & !name %in% paste0("logFC.", c(3, 8, 11))) %>%
  filter(abs(value) > 1) %>%
  distinct(id) %>%
  pull(id)

# intersect size
list(SE6 = down_se6, SC = degs_sc) %>%
  fromList() %>%
  upset()

down_intersect <- intersect(down_se6, degs_sc)

heatmap_mat <- ring %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = down_intersect, melted = TRUE,
                exprs_values = "logcounts", zeros_to_na = FALSE) %>%
  filter(!cluster_mnn_logvst %in% c(3, 8, 11)) %>%
  group_by(id, cluster_mnn_logvst) %>%
  summarise(expr = mean(expr)) %>%
  pivot_wider(names_from = "cluster_mnn_logvst",
              values_from = expr) %>%
  column_to_rownames("id") %>%
  as.matrix() %>%
  # scale
  t() %>% scale() %>% t()

gene_clust <- heatmap_mat %>%
  dist() %>%
  hclust() %>%
  cutree(7) %>%
  factor() %>%
  enframe(value = "cluster") %>%
  column_to_rownames("name")
pheatmap(heatmap_mat[order(gene_clust$cluster), ],
         cutree_rows = 7,
         annotation_row = gene_clust,
         main = "LFC < -1 in SE6",
         show_rownames = FALSE,
         height = 10, width = 7)



# Trajectory patterns -----------------------------------------------------
library(mgcv)

sling <- read_rds("data/processed/trajectories/ring_hardfilt_slingshot.rds")
traj <- read_csv("data/processed/trajectories/slingshot_trajectory_gam.csv") %>%
  filter(trajectory %in% paste0("trajectory", c(1, 3, 5)))

# distribution of r-squared values from GAM model
traj %>%
  ggplot(aes(rsq, trajectory)) +
  ggridges::geom_density_ridges(fill = "lightgrey", alpha = 0.5) +
  ggridges::theme_ridges()

# Intersection for genes with rsq > 0.7
traj %>%
  filter(rsq > 0.7) %>%
  with(split(gene, trajectory)) %>%
  fromList() %>%
  upset(order.by = "freq")

# get that subset of genes
top_genes <- traj %>%
  filter(rsq > 0.7) %>%
  distinct(gene, trajectory)

# getting predicted expression pattern

# prepare data for the model
traj_fit <- lapply(c("1", "3", "5"), function(pseudotime){
  # vector of genes relevant for this trajectory
  relevant_genes <- top_genes %>%
    filter(trajectory == paste0("trajectory", pseudotime)) %>%
    pull(gene) %>% unique()

  # subset the single cell object
  sce <- sling[relevant_genes, ]

  # get the pseudo-time
  t <- sce[[paste0("slingPseudotime_",pseudotime)]]
  names(t) <- colnames(sce)
  t <- t[which(!is.na(t))]
  t <- (t - min(t))/(max(t) - min(t)) # scale between 0-1

  # get the expression data
  Y <- as.matrix(assay(sce[, names(t)], "logvst"))
  Y <- Y %>% t() %>% scale() %>% t()

  # equally spaced timepoints for prediction
  newdata <- tibble(t = seq(0, 1, 0.02))

  # produce matrix of predictions
  out <- matrix(nrow = length(newdata$t), ncol = nrow(Y))
  colnames(out) <- rownames(Y)
  for(i in rownames(Y)){
    message(round(match(i, rownames(Y))/nrow(Y)*100, 2), "%")
    fit <- gam(Y[i, names(t)] ~ s(t))
    out[, i] <- predict(fit, newdata = newdata)
  }

  # return the matrix of predictions
  out
})

names(traj_fit) <- paste0("trajectory", c("1", "3", "5"))

traj_clusters <- traj_fit %>%
  map_dfr(function(x){
    clust <- x %>% t() %>% dist() %>% hclust() %>%
      cutree(k = 6) %>%
      enframe(name = "gene", value = "gene_clust")

    # clust <- 1 - cor(x)
    # clust <- clust %>% as.dist() %>% hclust() %>%
    #   cutree(k = 6) %>%
    #   enframe(name = "gene", value = "gene_clust")

    # clust <- x %>% t() %>% cluster::pam(k = 6)
    # clust <- clust$clustering %>% enframe(name = "gene", value = "gene_clust")

    cbind(t = seq(0, 1, 0.02), x) %>%
      as_tibble() %>%
      pivot_longer(-t, names_to = "gene", values_to = "pred") %>%
      full_join(clust, by = "gene")
  }, .id = "trajectory")

traj_clusters %>%
  ggplot(aes(t, pred)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "brown") +
  geom_line(aes(group = gene), alpha = 0.5) +
  facet_grid(trajectory ~ gene_clust, scales = "free") +
  theme_classic()


# Trajectory 1 ------------------------------------------------------------

p1 <- sling %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(pseudotime = (slingPseudotime_1 - min(slingPseudotime_1, na.rm = TRUE))/(max(slingPseudotime_1, na.rm = TRUE) - min(slingPseudotime_1, na.rm = TRUE))) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = pseudotime), alpha = 0.1) +
  geom_point(stat = "centroid", aes(group = cluster_mnn_logvst), shape = 21, size = 5,
             fill = "white", alpha = 0.7) +
  geom_text(stat = "centroid",
            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst)) +
  scale_colour_viridis_c(option = "inferno", na.value = "lightgrey") +
  coord_fixed() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_void() +
  theme(legend.position = "none")

p1.2 <- sling %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(pseudotime = (slingPseudotime_1 - min(slingPseudotime_1, na.rm = TRUE))/(max(slingPseudotime_1, na.rm = TRUE) - min(slingPseudotime_1, na.rm = TRUE))) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = cluster_mnn_logvst), alpha = 0.1) +
  geom_point(stat = "centroid", aes(group = cluster_mnn_logvst), shape = 21, size = 5,
             fill = "white", alpha = 0.7) +
  geom_text(stat = "centroid",
            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst)) +
  ggthemes::scale_colour_tableau("Tableau 20") +
  coord_fixed() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_void() +
  theme(legend.position = "none")


p2 <- traj_clusters %>%
  filter(trajectory == "trajectory1") %>%
  mutate(gene_clust = LETTERS[gene_clust]) %>%
  mutate(gene_up = ifelse(gene %in% up_se6, pred, NA),
         gene_down = ifelse(gene %in% down_se6, pred, NA)) %>%
  group_by(gene_clust, t) %>%
  mutate(n_up = sum(!is.na(gene_up)),
         n_down = sum(!is.na(gene_down)),
         n = n()) %>%
  ungroup() %>%
  mutate(clust_label = paste0(gene_clust, ") ", n_up, " Up; ", n_down, " Down")) %>%
  ggplot(aes(t, pred, group = gene)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(alpha = 0.5) +
  geom_line(aes(y = gene_down, colour = "Down")) +
  geom_line(aes(y = gene_up, colour = "Up")) +
  facet_wrap( ~ clust_label) +
  theme_classic() +
  scale_colour_manual(values = c("steelblue", "firebrick")) +
  labs(x = "Pseudotime", y = "Average Expression", colour = "SE6\nmutant")

temp <- sling %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(pseudotime = (slingPseudotime_1 - min(slingPseudotime_1, na.rm = TRUE))/(max(slingPseudotime_1, na.rm = TRUE) - min(slingPseudotime_1, na.rm = TRUE)))

p2 <- traj_clusters %>%
  filter(trajectory == "trajectory1") %>%
  mutate(gene_clust = LETTERS[gene_clust]) %>%
  mutate(gene_up = ifelse(gene %in% up_se6, pred, NA),
         gene_down = ifelse(gene %in% down_se6, pred, NA)) %>%
  group_by(gene_clust, t) %>%
  mutate(n_up = sum(!is.na(gene_up)),
         n_down = sum(!is.na(gene_down)),
         n = n()) %>%
  ungroup() %>%
  mutate(clust_label = paste0(gene_clust, ") ", n_up, " Up; ", n_down, " Down")) %>%
  ggplot(aes(t, pred, group = gene)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(alpha = 0.5) +
  geom_line(aes(y = gene_down), colour = "steelblue") +
  geom_line(aes(y = gene_up), colour = "firebrick") +
  geom_point(data = temp, aes(x = pseudotime, y = -2, colour = cluster_mnn_logvst), shape = "|", inherit.aes = FALSE, size = 2) +
  facet_wrap( ~ clust_label) +
  ggthemes::scale_colour_tableau("Tableau 20") +
  theme_classic() + theme(legend.position = "none") +
  labs(x = "Pseudotime", y = "Average Expression", colour = "SE6\nmutant")

((p1 / p1.2) | p2) +
  plot_annotation(title = "PPP trajectory") +
  plot_layout(widths = c(1, 2))
