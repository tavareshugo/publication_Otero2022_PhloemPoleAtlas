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



# Visualise ---------------------------------------------------------------

# pseudotime projected on UMAP
sling %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  pivot_longer(matches("^slingPseudotime"),
               values_to = "pseudotime", names_to = "trajectory") %>%
  mutate(trajectory = gsub("slingPseudotime_", "Trajectory ", trajectory)) %>%
  arrange(trajectory, !is.na(pseudotime), pseudotime) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = pseudotime), alpha = 0.1) +
  geom_point(stat = "centroid", aes(group = cluster_mnn_logvst), shape = 21, size = 5,
             fill = "white", alpha = 0.7) +
  geom_text(stat = "centroid",
            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst)) +
  facet_grid(~ trajectory) +
  scale_colour_viridis_c(option = "inferno", na.value = "lightgrey") +
  coord_fixed() +
  theme_void()

# number of cells in each trajectory/cluster
cell_weights <- slingCurveWeights(SlingshotDataSet(sling))
cell_trajs <- tibble(
  cell = colnames(sling),
  trajectory = colnames(cell_weights)[max.col(cell_weights)]
)

sling %>%
  colData() %>%
  as_tibble(rownames = "cell") %>%
  full_join(cell_trajs, by = "cell") %>%
  select(cluster_mnn_logvst, trajectory,
         matches("slingPseudotime")) %>%
  pivot_longer(matches("sling")) %>%
  mutate(name = str_replace(name, "slingPseudotime_", "curve")) %>%
  filter(trajectory == name) %>%
  count(cluster_mnn_logvst, name) %>%
  group_by(cluster_mnn_logvst) %>%
  mutate(pct = n/sum(n)*100) %>%
  ungroup() %>%
  ggplot(aes(cluster_mnn_logvst, pct)) +
  geom_col(aes(fill = name %>% str_remove("curve"))) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Cluster", y = "%", fill = "Trajectory")


# Probability assignment --------------------------------------------------

cell_weights <- slingCurveWeights(SlingshotDataSet(sling), as.probs = TRUE)
cell_trajs <- tibble(
  cell = colnames(sling),
  trajectory = colnames(cell_weights)[max.col(cell_weights)]
)

cell_weights %>%
  as_tibble(rownames = "cell") %>%
  full_join(sling %>% colData() %>% as_tibble(rownames = "cell"),
            by = "cell") %>%
  select(matches("curve"), cluster_mnn_logvst) %>%
  pivot_longer(matches("curve")) %>%
  mutate(name = str_replace(name, "curve", "Trajectory ")) %>%
  ggplot(aes(value)) +
  geom_density(aes(colour = name)) +
  # facet_wrap(~ cluster_mnn_logvst, scales = "free") +
  facet_grid(cluster_mnn_logvst ~ name, scales = "free") +
  scale_colour_brewer(palette = "Dark2")

cell_weights %>%
  as_tibble(rownames = "cell") %>%
  full_join(sling %>% colData() %>% as_tibble(rownames = "cell"),
            by = "cell") %>%
  select(matches("curve"), cluster_mnn_logvst) %>%
  pivot_longer(matches("curve")) %>%
  mutate(name = str_replace(name, "curve", "Trajectory ")) %>%
  filter(cluster_mnn_logvst %in% c(3, 4, 5, 8)) %>%
  ggplot(aes(value)) +
  geom_density(aes(colour = name)) +
  facet_wrap(~ cluster_mnn_logvst, scales = "free_y", nrow = 1) +
  # facet_grid(cluster_mnn_logvst ~ name, scales = "free") +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Probability of Cell Assigned to Trajectory", y = "") +
  theme(axis.text.y = element_blank())


# Genes with differential patterns ----------------------------------------

# score based on association with any trajectory
gois <- association %>%
  # arrange(desc(waldStat)) %>%
  # slice(1:500) %>%
  filter(waldStat >= quantile(waldStat, 0.95, na.rm = TRUE)) %>%
  left_join(gamfits) %>% filter(rsq > 0.6) %>%
  pull(gene) %>% unique()

set.seed(20210223)
gene_clusters <- predictSmooth(tradeseq, gois) %>%
  t() %>% scale() %>% t() %>%
  kmeans(centers = 9, nstart = 100, iter.max = 100) %>%
  `$`(cluster) %>%
  enframe("gene", "cluster") %>%
  left_join(predict_tradeseq(tradeseq, gois, scale = TRUE)) %>%
  group_by(cluster) %>%
  # mutate(cluster = paste0(LETTERS[cluster], "\n(", n_distinct(gene), ")"))
  mutate(ngenes = n_distinct(gene)) %>%
  ungroup() %>%
  mutate(cluster = LETTERS[cluster])

# cluster order for pretty plot
cluster_order <- tribble(
  ~cluster, ~order,
  "A", 1,
  "G", 2,
  "H", 3,
  "C", 4,
  "I", 5,
  "B", 6,
  "D", 7,
  "E", 8,
  "F", 9
)

# re-name the clusters with new order
gene_clusters <- gene_clusters %>%
  full_join(cluster_order, by = "cluster") %>%
  mutate(cluster = paste0(LETTERS[10-order], "\n(", ngenes, ")"))

pdf("documents/pdf for figures/trajectory_gene_clusters.pdf",
    width = 3, height = 4)
gene_clusters %>%
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
dev.off()

# Ring genes --------------------------------------------------------------

# how many ring genes are in the list used for clustering?
gene_clusters %>%
  inner_join(ring_genes, by = "gene") %>%
  distinct(gene, name, cluster, domain)

# check their patterns
predict_tradeseq(tradeseq, ring_genes$gene, scale = TRUE) %>%
  inner_join(ring_genes, by = "gene") %>%
  ggplot(aes(t, .pred)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(group = interaction(gene, trajectory), colour = trajectory)) +
  facet_grid( ~ name, scales = "free_y") +
  scale_colour_brewer(palette = "Dark2") +
  theme_void(base_size = 16)

# check the rank of ring genes
association %>%
  select(gene, waldStat) %>%
  mutate(rank = rank(waldStat),
         quantile = ecdf(waldStat)(waldStat)) %>%
  filter(gene %in% ring_genes$gene)

# check their expression
sling %>%
  getReducedDim("UMAP30_MNN_logvst", genes = ring_genes$gene,
                zeros_to_na = FALSE) %>%
  filter(id %in% gene_clusters$gene) %>%
  select(id, matches("slingP"), expr) %>%
  pivot_longer(matches("slingP"), names_to = "trajectory", values_to = "t") %>%
  drop_na() %>%
  mutate(trajectory = str_replace(trajectory, "slingPseudotime_", "Traj")) %>%
  filter(trajectory %in% paste0("Traj", 1:3)) %>%
  left_join(ring_genes, by = c("id" = "gene")) %>%
  ggplot(aes(t, expr)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(se = FALSE) +
  facet_grid(trajectory ~ name) +
  theme_minimal() +
  labs(x = "pseudotime", y = "logcounts")

# figure for paper
p1 <- sling %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = ring_genes$gene) %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  select(cell, cluster_mnn_logvst, matches("slingPseudotime"), id, expr) %>%
  pivot_longer(cols = matches("slingPseudotime"),
               names_to = "traj", values_to = "pseudotime") %>%
  drop_na() %>%
  mutate(traj = str_replace(traj, "slingPseudotime_", "Traj")) %>%
  filter(traj %in% paste0("Traj", 1:3)) %>%
  left_join(ring_genes, by = c("id" = "gene")) %>%
  ggplot(aes(pseudotime, expr)) +
  geom_point(aes(colour = cluster_mnn_logvst), size = 0.5) +
  geom_smooth(se = FALSE, colour = "black",
              method = mgcv::gam, formula = y ~ s(x, bs = "cr", k = 7)) +
  facet_grid(traj ~ name) +
  ggthemes::scale_colour_tableau("Tableau 20") +
  guides(colour = "none") +
  labs(colour = "Cluster", y = "logcounts") +
  theme_minimal()

p2 <- sling %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = ring_genes$gene) %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  select(cell, cluster_mnn_logvst, matches("slingPseudotime"), id, expr) %>%
  pivot_longer(cols = matches("slingPseudotime"),
               names_to = "traj", values_to = "pseudotime") %>%
  drop_na() %>%
  mutate(traj = str_replace(traj, "slingPseudotime_", "Traj")) %>%
  filter(traj %in% paste0("Traj", 1:3)) %>%
  left_join(ring_genes, by = c("id" = "gene")) %>%
  ggplot(aes(pseudotime, expr)) +
  geom_smooth(aes(colour = traj),
              se = FALSE,
              method = mgcv::gam, formula = y ~ s(x, bs = "cr", k = 7)) +
  facet_wrap( ~ name) +
  scale_colour_brewer(palette = "Dark2") +
  labs(colour = "Trajectory", y = "logcounts") +
  theme_minimal()

pdf("documents/pdf for figures/trajectories_ring_genes.pdf",
    width = 8.5, height = 3)
p1 + theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(0, 0.2, by = 0.1))
dev.off()

pdf("documents/pdf for figures/trajectories_ring_genes_compact.pdf",
    width = 4.68, height = 4.68)
p2 + theme(panel.grid.minor = element_blank(),
           legend.position = c(.95, 0),
           legend.justification = c(1, 0)) +
  scale_x_continuous(breaks = seq(0, 0.2, by = 0.1))
dev.off()


# SE6 mutant --------------------------------------------------------------

# sixtuple mutant results from previous work
se6 <- read_excel("data/raw/SE6vsCol_DEG.xlsx",
                  sheet = "SE6vsCol.DEG",
                  range = cell_cols("A:F")) %>%
  rename(gene = Gene_id, se6 = readcount_SE6, wt = readcount_Col) %>%
  mutate(direction = ifelse(log2FoldChange < 0, "Down", "Up"))

se6_gois <- se6 %>%
  # filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
  filter(gene %in% rownames(tradeseq)) %>%
  pull(gene) %>% unique()

# how many genes in each cluster
se6 %>%
  # filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
  filter(gene %in% rownames(tradeseq)) %>%
  right_join(gene_clusters) %>%
  distinct(gene, cluster, direction) %>%
  count(cluster, direction) %>%
  mutate(direction = ifelse(is.na(direction), "Not DE", direction)) %>%
  mutate(direction = fct_relevel(direction, "Not DE")) %>%
  group_by(cluster) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(pct, cluster)) +
  geom_col(aes(fill = direction)) +
  scale_fill_manual(values = c("grey", "steelblue", "firebrick")) +
  labs(x = "% genes in cluster", y = "Cluster", fill = "",
       subtitle = "DEGs in PEAR sixtuple mutant")

# estimate expected proportion of genes in each group based on a permutation test
# first cross-reference the se6 list with the gene clusters
se6_in_clusters <- se6 %>%
  select(gene, direction) %>%
  right_join(gene_clusters, by = "gene") %>%
  distinct(gene, cluster, direction) %>%
  mutate(direction = ifelse(is.na(direction), "Not DE", direction)) %>%
  mutate(direction = fct_relevel(direction, "Not DE"),
         cluster = str_remove(cluster, "\n.*"))

# do the permutation
set.seed(20210302)
se6_samples <- replicate(
  n = 1000,
  se6_in_clusters %>%
    # shuffle clusters
    mutate(cluster = sample(cluster)) %>%
    count(cluster, direction) %>%
    group_by(cluster) %>%
    mutate(pct = n/sum(n)*100) %>%
    ungroup(),
  simplify = FALSE
) %>%
  bind_rows(.id = "iter")

# summarise permutation results
se6_samples <- se6_samples %>%
  filter(direction != "Not DE") %>%
  mutate(pct = ifelse(direction == "Down", pct*-1, pct)) %>%
  group_by(cluster, direction) %>%
  summarise(lo = quantile(pct, 0.025),
            mid = median(pct),
            hi = quantile(pct, 0.975))

# visualise
se6_in_clusters %>%
  count(cluster, direction) %>%
  group_by(cluster) %>%
  mutate(pct = n/sum(n)*100) %>%
  ungroup() %>%
  filter(direction != "Not DE") %>%
  mutate(pct = ifelse(direction == "Down", pct*-1, pct),
         cluster = str_remove(cluster, "\n.*")) %>%
  ggplot(aes(y = cluster)) +
  geom_segment(aes(x = 0, xend = pct, y = cluster, yend = cluster),
               colour = "darkgrey") +
  geom_vline(xintercept = 0, size = 1) +
  # geom_point(aes(x = pct), size = 3) +
  geom_label(aes(x = pct, label = n, fill = direction)) +
  geom_linerange(data = se6_samples,
                 aes(xmin = lo, xmax = hi, colour = direction),
                 size = 10, alpha = 0.2) +
  scale_colour_manual(values = c("Down" = "steelblue", "Up" = "firebrick")) +
  scale_fill_manual(values = c("Down" = "steelblue", "Up" = "firebrick")) +
  labs(x = "% DEGs in cluster", y = "Cluster", fill = "",
       subtitle = paste0(sum(se6$gene %in% gene_clusters$gene),
                         " (out of ", length(unique(se6$gene)),
                         ") DEGs from PEAR sixtuple mutant")) +
  scale_x_continuous(breaks = seq(-100, 100, by = 20),
                     labels = c(rev(seq(20, 100, by = 20)), seq(0, 100, by = 20))) +
  guides(fill = "none") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())




# CCC mutant --------------------------------------------------------------

# triple mutant in genes that are presumed targets of PEARs
ccc <- read_xlsx("data/raw/ccc_vs_wt_upregulated.xlsx") %>%
  select(gene = Gene, ccc = ccc.value, wt = Wt.value, log2FoldChange = log2foldchange,
         pvalue, padj = padjust) %>%
  mutate(direction = ifelse(log2FoldChange < 0, "Down", "Up"))

ccc_gois <- ccc %>%
  # filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
  filter(gene %in% rownames(tradeseq)) %>%
  pull(gene) %>% unique()

# how many genes in each cluster
ccc %>%
  # filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
  filter(gene %in% rownames(tradeseq)) %>%
  right_join(gene_clusters) %>%
  distinct(gene, cluster, direction) %>%
  count(cluster, direction) %>%
  mutate(direction = ifelse(is.na(direction), "Not DE", direction)) %>%
  mutate(direction = fct_relevel(direction, "Not DE")) %>%
  group_by(cluster) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(pct, cluster)) +
  geom_col(aes(fill = direction)) +
  scale_fill_manual(values = c("grey", "steelblue", "firebrick")) +
  labs(x = "% genes in cluster", y = "Cluster", fill = "",
       subtitle = "DEGs in PEAR sixtuple mutant")

# estimate expected proportion of genes in each group based on a permutation test
# first cross-reference the ccc list with the gene clusters
ccc_in_clusters <- ccc %>%
  select(gene, direction) %>%
  right_join(gene_clusters, by = "gene") %>%
  distinct(gene, cluster, direction) %>%
  mutate(direction = ifelse(is.na(direction), "Not DE", direction)) %>%
  mutate(direction = fct_relevel(direction, "Not DE"),
         cluster = str_remove(cluster, "\n.*"))

# do the permutation
set.seed(20210302)
ccc_samples <- replicate(
  n = 1000,
  ccc_in_clusters %>%
    # shuffle clusters
    mutate(cluster = sample(cluster)) %>%
    count(cluster, direction) %>%
    group_by(cluster) %>%
    mutate(pct = n/sum(n)*100) %>%
    ungroup(),
  simplify = FALSE
) %>%
  bind_rows(.id = "iter")

# summarise permutation results
ccc_samples <- ccc_samples %>%
  filter(direction != "Not DE") %>%
  mutate(pct = ifelse(direction == "Down", pct*-1, pct)) %>%
  group_by(cluster, direction) %>%
  summarise(lo = quantile(pct, 0.025),
            mid = median(pct),
            hi = quantile(pct, 0.975))

# visualise
ccc_in_clusters %>%
  count(cluster, direction) %>%
  group_by(cluster) %>%
  mutate(pct = n/sum(n)*100) %>%
  ungroup() %>%
  filter(direction != "Not DE") %>%
  mutate(pct = ifelse(direction == "Down", pct*-1, pct),
         cluster = str_remove(cluster, "\n.*")) %>%
  ggplot(aes(y = cluster)) +
  geom_segment(aes(x = 0, xend = pct, y = cluster, yend = cluster),
               colour = "darkgrey") +
  geom_vline(xintercept = 0, size = 1) +
  # geom_point(aes(x = pct), size = 3) +
  geom_label(aes(x = pct, label = n, fill = direction)) +
  geom_linerange(data = ccc_samples,
                 aes(xmin = lo, xmax = hi, colour = direction),
                 size = 10, alpha = 0.2) +
  scale_colour_manual(values = c("Down" = "steelblue", "Up" = "firebrick")) +
  scale_fill_manual(values = c("Down" = "steelblue", "Up" = "firebrick")) +
  labs(x = "% DEGs in cluster", y = "Cluster", fill = "",
       subtitle = paste0(sum(ccc$gene %in% gene_clusters$gene),
                         " (out of ", length(unique(ccc$gene)),
                         ") DEGs from PAPL triple mutant")) +
  scale_x_continuous(breaks = seq(-100, 100, by = 20),
                     labels = c(rev(seq(20, 100, by = 20)), seq(0, 100, by = 20))) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


# check for overall enrichment in these lists
length(intersect(se6_gois, unique(gene_clusters$gene)))
se6_expect <- replicate(1000,
          length(
            intersect(
              sample(rownames(tradeseq), length(se6_gois)),
              sample(rownames(tradeseq), length(unique(gene_clusters$gene)))
            )
          ))
quantile(se6_expect, c(0.025, 0.5, 0.975))

ccc_expect <- replicate(1000,
                        length(
                          intersect(
                            sample(rownames(tradeseq), length(ccc_gois)),
                            sample(rownames(tradeseq), length(unique(gene_clusters$gene)))
                          )
                        ))
quantile(ccc_expect, c(0.025, 0.5, 0.975))



# Network modules ---------------------------------------------------------

modules <- read_csv("~/Downloads/Module1_sub136.csv")

# predict from tradeseq fits
predict_tradeseq(tradeseq, unique(modules$gene), scale = TRUE) %>%
  inner_join(modules, by = "gene") %>%
  ggplot(aes(t, .pred)) +
  geom_line(aes(group = interaction(gene, trajectory)), alpha = 0.3) +
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 1) +
  geom_hline(yintercept = 0, colour = "steelblue") +
  facet_grid(trajectory ~ sub_cluster) +
  scale_colour_brewer(palette = "Dark2") +
  theme_void(base_size = 16)

# my own GAM output (filter for R2)
gamfits %>%
  inner_join(modules, by = "gene") %>%
  group_by(gene) %>%
  filter(any(rsq > 0.5)) %>%
  ungroup() %>%
  select(gene, trajectory, sub_cluster, pred) %>%
  unnest(pred) %>%
  group_by(gene) %>%
  mutate(.pred = (.pred - mean(.pred))/sd(.pred)) %>%
  ungroup() %>%
  ggplot(aes(t, .pred)) +
  geom_line(aes(group = interaction(gene, trajectory)), alpha = 0.1) +
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 1) +
  geom_hline(yintercept = 0, colour = "steelblue") +
  facet_grid(trajectory ~ sub_cluster) +
  scale_colour_brewer(palette = "Dark2") +
  theme_void(base_size = 16)

# colour by Wald statistic
gamfits %>%
  inner_join(modules, by = "gene") %>%
  group_by(gene) %>%
  # filter(any(rsq > 0.5)) %>%
  ungroup() %>%
  # select(gene, trajectory, sub_cluster, pred) %>%
  unnest(pred) %>%
  group_by(gene) %>%
  mutate(.pred = (.pred - mean(.pred))/sd(.pred)) %>%
  ungroup() %>%
  inner_join(association %>% select(gene, waldStat), by = "gene") %>%
  drop_na(waldStat) %>%
  ggplot(aes(t, .pred)) +
  geom_line(aes(group = fct_reorder(interaction(gene, trajectory, sub_cluster),
                                    rsq),
                colour = rsq),
            alpha = 0.8) +
  geom_line(stat = "summary", fun = "median", size = 1) +
  geom_hline(yintercept = 0, colour = "steelblue") +
  facet_grid(trajectory ~ sub_cluster) +
  scale_colour_viridis_c() +
  theme_void(base_size = 16)

# highlight ring genes
gamfits %>%
  inner_join(modules, by = "gene") %>%
  # group_by(gene) %>%
  # filter(any(rsq > 0.5)) %>%
  # ungroup() %>%
  unnest(pred) %>%
  group_by(gene) %>%
  mutate(.pred = (.pred - mean(.pred))/sd(.pred)) %>%
  ungroup() %>%
  left_join(ring_genes) %>%
  mutate(highlight = ifelse(!is.na(name), .pred, NA),
         sub_cluster = paste0("sub-module", sub_cluster)) %>%
  ggplot(aes(t, .pred)) +
  geom_line(aes(group = interaction(gene, trajectory, sub_cluster)),
            alpha = 0.2, colour = "grey48") +
  geom_line(aes(group = interaction(gene, trajectory, sub_cluster),
                y = highlight, colour = name)) +
  geom_line(stat = "summary", fun = "median", size = 1) +
  geom_hline(yintercept = 0, colour = "steelblue") +
  facet_grid(trajectory ~ sub_cluster) +
  theme_void(base_size = 16)


# Optimal number of clusters -----------------------------------------------

# find optimal number of clusters
clust_matrix <- predictSmooth(tradeseq, gois) %>%
  t() %>% scale() %>% t()
kmeans_test <- tibble(k = seq(1, 30, by = 1)) %>%
  mutate(kmeans = map(k, ~ kmeans(clust_matrix, centers = .x,
                                  iter.max = 100, nstart = 10))) %>%
  mutate(kmeans = map(kmeans, broom::glance)) %>%
  unnest(kmeans)

kmeans_test %>%
  arrange(k) %>%
  mutate(tot.withinss = tot.withinss) %>%
  ggplot(aes(k, tot.withinss)) +
  geom_line() +
  geom_label(aes(label = k))




########## DEPRECATED ##########

# Clustering HGVs ---------------------------------------------------------

# clustering all 2000 HVGs
predictSmooth(tradeseq, metadata(sling)$hvgs) %>%
  t() %>% scale() %>% t() %>%
  kmeans(centers = 9, nstart = 10, iter.max = 100) %>%
  `$`(cluster) %>%
  enframe("gene", "cluster") %>%
  left_join(predict_tradeseq(tradeseq, gois)) %>%
  group_by(gene) %>%
  mutate(.pred = (.pred - mean(.pred))/sd(.pred)) %>%
  ungroup() %>%
  ggplot(aes(t, .pred)) +
  geom_line(aes(group = interaction(gene, trajectory)),
            alpha = 0.2) +
  facet_grid(trajectory ~ cluster)


# Genes along trajectory --------------------------------------------------

gamfits %>%
  group_by(trajectory) %>%
  arrange(desc(rsq)) %>%
  slice(1:100) %>%
  mutate(trajectory = str_remove(trajectory, "trajectory")) %>%
  group_by(gene) %>%
  summarise(trajectory = paste(sort(trajectory), collapse = ";")) %>%
  write_csv("data/processed/gene_sets/top100_along_trajectories.csv")

top_genes <- gamfits %>%
  group_by(trajectory) %>%
  arrange(desc(rsq)) %>%
  slice(1:100) %>%
  ungroup()

top_genes <- split(top_genes$gene, top_genes$trajectory)

for(i in c(1:3)){
  heatdata <- assay(sling, "logvst")[top_genes[[paste0("trajectory", i)]],
                                         order(sling[[paste0("slingPseudotime_", i)]],
                                               na.last = NA)]
  heatdata <- apply(as.matrix(heatdata), 1,
                    function(x) (x - min(x))/(max(x) - min(x)))
  heatdata <- t(heatdata)
  clust_annot <- data.frame(cluster = sling$cluster_mnn_logvst)
  rownames(clust_annot) <- colnames(sling)
  pheatmap::pheatmap(as.matrix(heatdata),
                     cluster_cols = FALSE,
                     show_colnames = FALSE,
                     show_rownames = FALSE,
                     annotation_col = clust_annot,
                     main = paste0("Trajectory ", i))

  # clean environment
  rm(heatdata, clust_annot, i)
}

# intersection
UpSetR::upset(UpSetR::fromList(top_genes), order.by = "freq")

# line plots
for(i in c(1:3)){
  heatdata <- assay(sling, "logvst")[top_genes[[paste0("trajectory", i)]],
                                         order(sling[[paste0("slingPseudotime_", i)]],
                                               na.last = NA)]

  # partition genes into clusters
  kmeans_out <- cluster::pam(as.matrix(heatdata), k = 5)
  kmeans_out <- kmeans_out$clustering %>%
    enframe(name = "gene", value = "k_cluster")

  cell_annot <- data.frame(cell = colnames(sling),
                           cell_cluster = sling$cluster_mnn_logvst,
                           trajectory1 = sling[[paste0("slingPseudotime_", i)]])

  p <- heatdata %>%
    as.matrix() %>%
    as_tibble(rownames = "gene") %>%
    pivot_longer(cols = -gene, values_to = "expr", names_to = "cell") %>%
    full_join(kmeans_out, by = "gene") %>%
    left_join(cell_annot, by = "cell") %>%
    ggplot(aes(trajectory1, expr)) +
    geom_smooth(aes(group = gene, colour = cell_cluster),
                se = FALSE,
                method = mgcv::gam,
                formula = y ~ s(x, bs = "cs"),
                colour = "black", alpha = 0.6) +
    facet_grid(k_cluster ~ .) +
    labs(title = paste0("Trajectory ", i))

  print(p)

  rm(heatdata, cell_annot, kmeans_out, i, p)
}


top_pred <- gamfits %>%
  group_by(gene) %>%
  filter(any(rsq > 0.6)) %>%
  ungroup() %>%
  select(trajectory, gene, pred) %>%
  unnest(pred)

top_pred <- top_pred %>%
  pivot_wider(names_from = c("trajectory", "t"), values_from = "y") %>%
  column_to_rownames("gene") %>%
  dist() %>%
  hclust(method = "ward.D2") %>%
  cutree(k = 8) %>%
  enframe("gene", "cluster") %>%
  full_join(top_pred, by = "gene")

top_pred %>%
  # group_by(trajectory) %>%
  # mutate(t = min_max(t, na.rm = TRUE)) %>%
  # ungroup() %>%
  ggplot(aes(t, y)) +
  geom_line(aes(group = interaction(gene, trajectory),
                colour = trajectory)) +
  facet_wrap(~ cluster)


# Cyclins -----------------------------------------------------------------

for(i in c(1, 3:5)){
  heatdata <- assay(sling, "logvst")[grep("CYC[A,B,D]",
                                              rowData(sling)$alternative_name,
                                              ignore.case = TRUE),
                                         order(sling[[paste0("slingPseudotime_", i)]],
                                               na.last = NA)]
  rownames(heatdata) <- rowData(sling)$alternative_name[grep("CYC[A,B,D]", rowData(sling)$alternative_name)]
  heatdata <- apply(heatdata, 1,
                    function(x) (x - min(x))/(max(x) - min(x)))
  heatdata <- t(heatdata)

  clust_annot <- data.frame(cluster = sling$cluster_mnn_logvst)
  rownames(clust_annot) <- colnames(sling)

  pheatmap::pheatmap(as.matrix(heatdata),
                     cluster_cols = FALSE,
                     show_colnames = FALSE,
                     show_rownames = TRUE,
                     annotation_col = clust_annot,
                     main = paste0("Cyclins; Trajectory ", i))

  # clean environment
  rm(heatdata, clust_annot, i)
}


# Cross-reference cluster markers -----------------------------------------

top_genes %>%
  map_dfr(enframe, name = NULL, value = "id", .id = "trajectory") %>%
  full_join(cluster_test, by = "id") %>%
  ggplot(aes(trajectory, AUC.1)) +
  geom_violin(scale = "width") +
  facet_wrap(~ cluster)

# This is similar to what we did in 04-cluster_markers
# for example picking diagnostic markers for PPP clusters
interesting_genes <- cluster_test %>%
  # FDR threshold
  filter(FDR < 0.05) %>%
  # for each gene
  group_by(id) %>%
  # retain genes that are significant in only certain clusters
  filter(n_distinct(cluster) == 2 &
           all(cluster %in% c(1, 4))) %>%
  ungroup() %>%
  # filter AUC threshold for better signal
  filter(summary.AUC > 0.8) %>%
  # arrange in descending order of AUC and get the gene IDs
  arrange(desc(summary.AUC)) %>%
  pull(id) %>% unique()

res %>%
  filter(gene %in% interesting_genes) %>%
  ggplot(aes(trajectory, rsq)) +
  geom_jitter(height = 0, width = 0.2)



# APL ---------------------------------------------------------------------

genes_of_interest <- tribble(
  ~name, ~gene,
  "APL", "AT1G79430",
  "NAC86", "AT5G17260",
  "NEN4", "AT4G39810")

genes_of_interest <- genes_of_interest %>%
  filter(gene %in% rownames(sling))

# plot genes along trajectories, cells coloured by shared trajectories
cbind(colData(sling),
      as.matrix(t(assay(sling, "logcounts")[genes_of_interest$gene, ]))) %>%
  as_tibble() %>%
  mutate(cell_id = paste(Sample, Barcode, sep = "_")) %>%
  select(cell_id,
         matches("slingPseudotime"), matches("AT.G"),
         cluster_mnn_logvst) %>%
  pivot_longer(cols = matches("AT.G"),
               names_to = "gene",
               values_to = "expr") %>%
  full_join(genes_of_interest, by = c("gene")) %>%
  pivot_longer(cols = matches("slingPseudo"),
               names_to = "trajectory", values_to = "pseudotime") %>%
  drop_na(pseudotime) %>%
  mutate(trajectory = str_remove(trajectory, "slingPseudotime_")) %>%
  filter(trajectory %in% c(1, 2, 3)) %>%
  mutate(trajectory = case_when(trajectory == 1 ~ "PPP",
                                trajectory == 2 ~ "CC",
                                trajectory == 3 ~ "SE",
                                TRUE ~ NA_character_)) %>%
  group_by(trajectory) %>%
  mutate(pseudotime = (pseudotime - min(pseudotime))/(max(pseudotime) - min(pseudotime))) %>%
  ungroup() %>%
  mutate(gene = ifelse(is.na(name), gene, name)) %>%
  # find common trajectories across cells
  group_by(cell_id) %>%
  mutate(cell_trajectories = paste(unique(trajectory), collapse = "; ")) %>%
  ungroup() %>%
  ggplot(aes(pseudotime, expr)) +
  geom_point(alpha = 0.5, aes(colour = cell_trajectories)) +
  geom_smooth(se = FALSE, colour = "black") +
  facet_grid(trajectory ~ gene) +
  theme(axis.text.x = element_blank()) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Pseudotime (min-max scaled)", y = "Expression",
       colour = "Cell\ntrajectories")

# plot genes along trajectories, cells coloured by cluster
cbind(colData(sling),
      as.matrix(t(assay(sling, "logcounts")[genes_of_interest$gene, ]))) %>%
  as_tibble() %>%
  mutate(cell_id = paste(Sample, Barcode, sep = "_")) %>%
  select(cell_id,
         matches("slingPseudotime"), matches("AT.G"),
         cluster_mnn_logvst) %>%
  pivot_longer(cols = matches("AT.G"),
               names_to = "gene",
               values_to = "expr") %>%
  full_join(genes_of_interest, by = c("gene")) %>%
  pivot_longer(cols = matches("slingPseudo"),
               names_to = "trajectory", values_to = "pseudotime") %>%
  drop_na(pseudotime) %>%
  mutate(trajectory = str_remove(trajectory, "slingPseudotime_")) %>%
  filter(trajectory %in% c(1, 2, 3)) %>%
  mutate(trajectory = case_when(trajectory == 1 ~ "PPP",
                                trajectory == 2 ~ "CC",
                                trajectory == 3 ~ "SE",
                                TRUE ~ NA_character_)) %>%
  # group_by(trajectory) %>%
  # mutate(pseudotime = (pseudotime - min(pseudotime))/(max(pseudotime) - min(pseudotime))) %>%
  # ungroup() %>%
  mutate(gene = ifelse(is.na(name), gene, name)) %>%
  # find common trajectories across cells
  group_by(cell_id) %>%
  mutate(cell_trajectories = paste(unique(trajectory), collapse = "; ")) %>%
  ungroup() %>%
  ggplot(aes(pseudotime, expr)) +
  geom_point(alpha = 0.5, aes(colour = cluster_mnn_logvst)) +
  geom_smooth(se = FALSE, colour = "black") +
  facet_grid(trajectory ~ gene) +
  theme(axis.text.x = element_blank()) +
  ggthemes::scale_colour_tableau(palette = "Tableau 20") +
  labs(x = "Pseudotime (min-max scaled)", y = "Expression",
       colour = "Cluster")


cell_weights <- slingCurveWeights(SlingshotDataSet(sling))[, c(1, 2, 3)]
cell_weights <- cell_weights[which(rowSums(cell_weights) > 0), ]
pseudotime <- slingPseudotime(SlingshotDataSet(sling), na = FALSE)[, c(1, 2, 3)]
pseudotime <- pseudotime[rownames(cell_weights), ]

cell_trajs <- tibble(cell = rownames(cell_weights),
                     traj = paste0("Traj", max.col(cell_weights)))

sling %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = genes_of_interest$gene) %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  select(cell, cluster_mnn_logvst, matches("slingPseudotime"), id, expr) %>%
  pivot_longer(cols = matches("slingPseudotime"), names_to = "traj", values_to = "pseudotime") %>%
  drop_na() %>%
  mutate(traj = str_replace(traj, "slingPseudotime_", "Traj")) %>%
  inner_join(cell_trajs, by = c("cell", "traj")) %>%
  ggplot(aes(pseudotime, expr)) +
  geom_point(aes(colour = cluster_mnn_logvst)) +
  # geom_smooth(se = FALSE, colour = "black") +
  facet_grid(traj ~ id, scales = "free", space = "free") +
  ggthemes::scale_colour_tableau("Tableau 20")

sling %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = genes_of_interest$gene) %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  select(cell, cluster_mnn_logvst, matches("slingPseudotime"), id, expr) %>%
  pivot_longer(cols = matches("slingPseudotime"), names_to = "traj", values_to = "pseudotime") %>%
  drop_na() %>%
  mutate(traj = str_replace(traj, "slingPseudotime_", "Traj")) %>%
  # group_by(traj) %>%
  # mutate(pseudotime = min_max(pseudotime, TRUE)) %>%
  # ungroup() %>%
  ggplot(aes(pseudotime, expr)) +
  geom_point(aes(colour = cluster_mnn_logvst)) +
  geom_smooth(se = FALSE, colour = "black") +
  facet_grid(traj ~ id, scales = "free", space = "free") +
  ggthemes::scale_colour_tableau("Tableau 20") +
  labs(colour = "Cluster", y = "logcounts")


# GAMs for comparing genes ----------------------------------------------

Y <- as.matrix(assay(sling, "logvst"))
t1 <- sling$slingPseudotime_1
t5 <- sling$slingPseudotime_5

library(mgcv)

target_gene <- "AT2G40890"
fit <- gam(Y[target_gene, ] ~ te(t1) + te(t5) + ti(t1, t5), method = "REML")
fit <- gam(Y[target_gene, ] ~ te(t1, t5), method = "REML")
fit <- gam(Y[target_gene, ] ~ s(t1, t5), method = "REML")
fit <- gam(Y[target_gene, ] ~ s(t5, k = 50), method = "REML")

# I think this could be a way to find genes with different profiles along these trajectories
d <- tibble(
  y = c(Y[target_gene, ], Y[target_gene, ]),
  time = c(t1, t5),
  traj = factor(rep(c("t1", "t5"), each = length(t1)))
) %>%
  drop_na() %>%
  # min-max normalisation
  group_by(traj) %>%
  mutate(time = (time - min(time))/(max(time) - min(time))) %>%
  ungroup()

# full model
fit <- gam(y ~ traj + s(time, by = traj, k = 50),
           data = d, method = "REML")

# model with no interaction
fit0 <- gam(y ~ traj + s(time, k = 50),
           data = d, method = "REML")
AIC(fit) - AIC(fit0)

summary(fit)
gam.check(fit)
plot(fit, pages = 1)
plot(t1, Y[target_gene, ])
plot(t5, Y[target_gene, ])
coef(fit)
vis.gam(fit, plot.type = "contour")
points(jitter(as.numeric(d$traj)), d$time, pch = 19)

res %>%
  filter(trajectory %in% c("trajectory1", "trajectory5")) %>%
  select(trajectory, gene, rsq) %>%
  pivot_wider(values_from = rsq, names_from = trajectory) %>%
  filter(trajectory1 > 0.7 | trajectory5 > 0.7) %>%
  arrange(-abs(trajectory1 - trajectory5))


# Testing GAMs ------------------------------------------------------------

# this is similar to what is implemented in tradeSeq
fit_gam <- function(Y, t, offsets = NULL){
  ngenes <- nrow(Y)
  gam_pvals <- vector("numeric", length = ngenes)
  gam_rsq <- vector("numeric", length = ngenes)
  gam_devexpl <- vector("numeric", length = ngenes)
  gam_aic <- vector("numeric", length = ngenes)
  pred <- vector("list", length = ngenes)

  for(i in 1:nrow(Y)){
    if(is.null(offsets)){
      fit <- mgcv::gam(Y[i,] ~ s(t, bs = "cr", k = 7), family = "nb")
    } else {
      fit <- mgcv::gam(Y[i,] ~ s(t, bs = "cr", k = 7) + offset(offsets), family = "nb")
    }

    gam_pvals[i] <- summary(fit)$s.table[4]
    gam_rsq[i] <- summary(fit)$r.sq
    gam_devexpl[i] <- summary(fit)$dev.expl
    gam_aic[i] <- AIC(fit)
    newd <- data.frame(t = seq(min(t, na.rm = TRUE), max(t, na.rm = TRUE), length.out = 100))
    if(!is.null(offsets)){
      # smooth offsets with time
      newd$offsets <- as.vector(predict(mgcv::gam(offsets ~ s(t)), newdata = newd))
    }
    pred[[i]] <- data.frame(t = newd$t,
                            .pred = predict(fit, newdata = newd, type = "response"))
    message("Gene ", i, " of ", ngenes)
  }

  tibble(gene = rownames(Y),
             pval = gam_pvals,
             rsq = gam_rsq,
             devexpl = gam_devexpl,
             aic = gam_aic,
             pred = pred)
}

target_gene <- "AT2G40890"
target_gene <- metadata(sling)$hvgs[1:10]
Y <- as.matrix(assay(sling[target_gene, ], "counts"))
t1 <- sling$slingPseudotime_1
log_umis <- log(colSums(counts(sling))*sizeFactors(sling))

res <- list(
  no = fit_gam(Y, t1),
  yes = fit_gam(Y, t1, log_umis)
)
res <- bind_rows(res, .id = "offset")

res %>%
  select(gene, pred, offset) %>%
  unnest(pred) %>%
  pivot_wider(names_from = "offset", values_from = ".pred") %>%
  ggplot(aes(no, yes)) +
  geom_path(aes(group = gene, colour = t)) +
  geom_abline(col = "red")

res %>%
  filter(gene %in% sample(target_gene, 10)) %>%
  select(gene, pred, offset) %>%
  unnest(pred) %>%
  ggplot(aes(t, .pred)) +
  geom_path(aes(colour = offset)) +
  facet_wrap( ~ gene)

qplot(t1, log_umis) + geom_smooth(method = "gam", formula = y ~ s(x)) +
  geom_hline(yintercept = mean(log_umis), colour = "brown", size = 1)

# comparing tradeSeq and my manual approach
ysmooth <- predictSmooth(models = tradeseq, gene = target_gene, nPoints = 100)
ysmooth <- ysmooth %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, values_to = ".pred") %>%
  separate(name, c("trajectory", "t"), convert = TRUE) %>%
  filter(trajectory == "lineage1") %>%
  select(-trajectory) %>%
  mutate(offset = "tradeSeq", method = "tradeSeq")

res %>%
  filter(gene %in% target_gene) %>%
  select(gene, pred, offset) %>%
  unnest(pred) %>%
  mutate(method = "manual") %>%
  # mutate(.pred = exp(.pred)) %>%
  bind_rows(ysmooth %>% mutate(.pred = log(.pred))) %>%
  group_by(gene, offset) %>%
  mutate(t = rank(t)) %>%
  ungroup() %>%
  ggplot(aes(t, .pred)) +
  geom_path(aes(colour = interaction(method, offset))) +
  facet_wrap( ~ gene, scales = "free")


expr <- sling %>%
  getReducedDim("UMAP30_MNN_logvst", genes = target_gene, exprs_values = "counts") %>%
  select(t = slingPseudotime_1, gene = id, expr) %>%
  drop_na()

res %>%
  select(gene, pred, offset) %>%
  unnest(pred) %>%
  ggplot(aes(t, log((.pred) + 1))) +
  geom_point(data = expr, aes(y = log(expr + 1)), alpha = 0.1) +
  geom_path(aes(colour = offset)) +
  facet_wrap( ~ gene)


expr <- sling %>%
  getReducedDim("UMAP30_MNN_logvst", genes = target_gene, exprs_values = "logcounts") %>%
  select(t = slingPseudotime_1, gene = id, expr) %>%
  drop_na()

result %>%
  filter(trajectory == "trajectory1") %>%
  select(gene, pred) %>%
  unnest(pred) %>%
  ggplot(aes(t, .pred)) +
  geom_point(data = expr, aes(y = expr), alpha = 0.1) +
  geom_line(colour = "brown") +
  facet_wrap( ~ gene)
