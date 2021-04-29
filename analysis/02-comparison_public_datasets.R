library(SingleCellExperiment)
library(tidyverse)
library(ggpointdensity)
library(patchwork)

# set ggplot2 theme
theme_set(theme_classic() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")

# helper plotting functions
umap_highlight_cluster <- function(sce, cluster){
  sce %>%
    getReducedDim("UMAP30_MNN_logvst") %>%
    mutate(Sample = str_remove(Sample, "MAKR5_")) %>%
    filter(cluster_mnn_logvst == cluster) %>%
    ggplot(aes(V1, V2)) +
    geom_point(data = getReducedDim(sce, "UMAP30_MNN_logvst") %>% select(-source),
               colour = "lightgrey", size = 0.3) +
    ggpointdensity::geom_pointdensity(size = 0.3) +
    scale_colour_viridis_c(option = "magma") +
    coord_equal() + theme_void() +
    guides(colour = "none")
}

umap_highlight_gene <- function(sce, gene){
  sce %>%
    getReducedDim("UMAP30_MNN_logvst",
                  genes = gene,
                  melted = TRUE,
                  exprs_values = "logcounts", zeros_to_na = TRUE) %>%
    arrange(!is.na(expr)) %>%
    group_by(cluster_mnn_logvst, id) %>%
    mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
    ggplot(aes(V1, V2)) +
    geom_point(aes(colour = expr_weighted), size = 0.3) +
    facet_wrap(~ id) +
    scale_colour_viridis_c(na.value = "lightgrey") +
    guides(colour = "none") +
    labs(x = "UMAP1", y = "UMAP2",
         colour = "Cluster-weighted\nNormalised\nExpression") +
    theme_void()
}


# Read data ---------------------------------------------------------------

# very large object - make sure to have enough RAM available
integrated <- readRDS("data/processed/SingleCellExperiment/integrated_batches_strictfilt_leaner.rds")

integrated$source <- ifelse(grepl("wendrich_", integrated$Sample),
                     "Wendrich et al.",
                     ifelse(grepl("shahan_", integrated$Sample), "Shahan et al.",
                            ifelse(grepl("denyer_", integrated$Sample), "Denyer et al.",
                                   "Ours")))

# list of curated genes
curated_genes <- read_csv("data/raw/signature_genes_circle_plot.csv") %>%
  select(tissue = Expression, id = gene_ID, name = gene_name) %>%
  mutate(tissue = str_replace(tissue, " \\+ ", "\n&\n")) %>%
  mutate(tissue = str_replace_all(tissue, " ", "\n")) %>%
  mutate(id = toupper(str_trim(id)))



# EDA ---------------------------------------------------------------------

scater::plotReducedDim(integrated, "UMAP30_MNN_logcounts",
                       colour_by = "cluster_mnn_logcounts",
                       text_by = "cluster_mnn_logcounts")
scater::plotReducedDim(integrated, "UMAP30_MNN_logvst",
                       colour_by = "cluster_mnn_logvst",
                       text_by = "cluster_mnn_logvst")

integrated %>%
  colData() %>%
  as_tibble() %>%
  ggplot(aes(cluster_mnn_logvst, cluster_mnn_logcounts)) +
  geom_count()

# UMAP with clusters
png("documents/pdf for figures/integrated_datasets/umap_clusters.png",
    width = 3, height = 3, units = "in", res = 300)
integrated %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(Sample = str_remove(Sample, "MAKR5_")) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = cluster_mnn_logvst), size = 0.3) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, size = 1, label.padding = unit(0.1, "lines")) +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "cluster") +
  scale_colour_viridis_d() +
  coord_equal() + theme_void() +
  guides(colour = "none")
dev.off()

# UMAP with samples split by origin
png("documents/pdf for figures/integrated_datasets/umap_by_source.png", width = 7.5, height = 2, units = "in", res = 300)
integrated %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  group_by(source) %>%
  mutate(n_source = n()) %>%
  ungroup() %>%
  mutate(pct_source = n_source/n()*100) %>%
  mutate(source = paste0(source, "\n", n_source, " cells (", round(pct_source), "%)")) %>%
  mutate(source = fct_reorder(source, pct_source)) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(integrated, "UMAP30_MNN_logvst") %>% select(-source),
             colour = "lightgrey", size = 0.3) +
  ggpointdensity::geom_pointdensity(size = 0.3) +
  # geom_label(data = getReducedDim(integrated, "UMAP30_MNN_logvst") %>% select(-source),
  #            stat = "centroid",
  #            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
  #            alpha = 0.8, size = 0.5, label.padding = unit(0.1, "lines")) +
  scale_colour_viridis_c(option = "magma") +
  coord_equal() + theme_void() +
  guides(colour = "none") +
  facet_grid(~ source)
dev.off()

# cell type annotation plot
pdf("documents/pdf for figures/integrated_datasets/cell_type_annotation_integrated.pdf",
    width = 7, height = 5)
integrated %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = curated_genes$id, melted = TRUE,
                exprs_values = "logcounts", zeros_to_na = FALSE) %>%
  left_join(curated_genes, by = "id") %>%
  # scale expression for each gene
  group_by(id) %>%
  mutate(expr_center = (expr - mean(expr))/sd(expr)) %>%
  # calculate percentage of cells where the gene is detected per cluster
  group_by(cluster_mnn_logvst, id, name, tissue) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  # tidy up the tissue column for nicer plotting
  filter(tissue != "early") %>%
  mutate(tissue = factor(tissue,
                         levels = c("CELL\nCYCLE\nGENES",
                                    "CC", "CC\n&\nMSE",
                                    "PPP\n&\nCC", "PPP",
                                    "EARLY\nPSE",
                                    "LATE\nPSE", "OUTER\nLAYERS"))) %>%
  # clip values for clearer plotting
  mutate(mean_expr = ifelse(abs(mean_expr) > 4, sign(mean_expr)*4, mean_expr)) %>%
  # plot
  ggplot(aes(name, factor(cluster_mnn_logvst))) +
  geom_point(aes(colour = mean_expr, size = pct_expressing)) +
  facet_grid(~ tissue, scales = "free", space = "free") +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "gray95",
                         limits = c(-4, 4), breaks = seq(-4, 4, by = 2),
                         labels = c("<= -4", -2, 0, 2, ">= 4")) +
  scale_size_continuous(limits = c(0, 100)) +
  guides(shape = "none") +
  labs(y = "Cluster", x = "", size = "% Detected",
       colour = "Z-score\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# individual plots for each gene
for(i in curated_genes$id){
  png(paste0("documents/pdf for figures/integrated_datasets/umap_", i, ".png"), width = 3, height = 2, units = "in", res = 300)
  print(umap_highlight_gene(integrated, i))
  dev.off()
}

# barplot with pct of cells where each gene is detected in
pdf("documents/pdf for figures/integrated_datasets/barplot_pct_detected_per_dataset_vertical.pdf",
    width = 4, height = 6)
integrated %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = curated_genes$id,
                melted = TRUE,
                exprs_values = "logcounts",
                zeros_to_na = TRUE) %>%
  group_by(source, id) %>%
  summarise(detected = sum(!is.na(expr))) %>%
  group_by(id) %>%
  mutate(pct = detected/sum(detected)*100) %>%
  ungroup() %>%
  inner_join(curated_genes, by = "id") %>%
  mutate(source = fct_relevel(source, "Ours")) %>%
  ggplot(aes(pct, name)) +
  geom_col(aes(fill = source)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  # facet_wrap(tissue ~ ., ncol = 1, scales = "free_y",
  #            strip.position = "top") +
  facet_grid(tissue~ ., scales = "free_y", space = "free_y") +
  labs(y = "Gene", x = "% cells detected", fill = "Dataset") +
  theme(strip.text = element_text(size = 8),
        strip.background = element_rect(colour = NA, fill = NA))
dev.off()

# barplot with pct of cells where each gene is detected in
pdf("documents/pdf for figures/integrated_datasets/barplot_pct_detected_per_dataset_horizontal.pdf",
    width = 7, height = 3)
integrated %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = curated_genes$id,
                melted = TRUE,
                exprs_values = "logcounts",
                zeros_to_na = TRUE) %>%
  group_by(source, id) %>%
  summarise(detected = sum(!is.na(expr))) %>%
  group_by(id) %>%
  mutate(pct = detected/sum(detected)*100) %>%
  ungroup() %>%
  inner_join(curated_genes, by = "id") %>%
  mutate(source = fct_relevel(source, "Ours")) %>%
  ggplot(aes(name, pct)) +
  geom_col(aes(fill = source)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  # facet_wrap(tissue ~ ., ncol = 1, scales = "free_y",
  #            strip.position = "top") +
  facet_grid(~ tissue, scales = "free_x", space = "free_x") +
  labs(x = "Gene", y = "% cells detected", fill = "Dataset") +
  theme(strip.text = element_text(size = 8),
        strip.background = element_rect(colour = NA, fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# Highlight clusters ---------------------------------------------------


# highlight different clusters
p1 <- integrated %>%
  umap_highlight_cluster(3) +
  labs(title = "Cluster 3 (dividing cells)")

p2 <- integrated %>%
  umap_highlight_cluster(23) +
  labs(title = "Cluster 23 (CC cells - late?)")

p3 <- integrated %>%
  umap_highlight_cluster(26) +
  labs(title = "Cluster 26 (PPP cells)")

p4 <- integrated %>%
  umap_highlight_cluster(28) +
  labs(title = "Cluster 28 (PSE cells)")

p5 <- integrated %>%
  umap_highlight_cluster(2) +
  labs(title = "Cluster 2 (Epidermal cells)")

p6 <- integrated %>%
  umap_highlight_cluster(18) +
  labs(title = "Cluster 18 (Root cap cells)")

# png("documents/")
wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3)


# what are these clusters?
p1 <- integrated %>%
  umap_highlight_cluster(25) +
  labs(title = "Cluster 25")

p2 <- integrated %>%
  umap_highlight_cluster(6) +
  labs(title = "Cluster 6")

p3 <- integrated %>%
  umap_highlight_cluster(19) +
  labs(title = "Cluster 19")

wrap_plots(p1, p2, p3, ncol = 3)

# proportion of cells in each cluster
integrated %>%
  colData() %>%
  as_tibble() %>%
  count(cluster = cluster_mnn_logvst, source) %>%
  group_by(cluster) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup() %>%
  mutate(source = fct_relevel(source, "Ours")) %>%
  ggplot(aes(pct, cluster)) +
  geom_col(aes(fill = source)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "% cells by dataset", fill = "Dataset")



# Mixing of cells in clusters ---------------------------------------------
library(igraph)

# build shared neighbours graph
# I pre-calculated this due to high memory requirements
# code commented-out is what was run on the compute server
graph <- readRDS("data/processed/SingleCellExperiment/integrated_batches_strictfilt_graph.rds")
V(graph)$source <- gsub(" .*", "", V(graph)$source)

# this is the same graph used to generate the clusters used above
# graph <- scran::buildSNNGraph(integrated, k = 50, type = "jaccard",
#                               use.dimred = "MNN_logvst")

# # get those clusters (and confirm they are indeed the same)
# clusters <- factor(igraph::cluster_louvain(graph)$membership)
# all.equal(clusters, integrated$cluster_mnn_logvst)
#
# # add metadata to graph
# V(graph)$cell <- colnames(integrated)[as.vector(V(graph))]
# V(graph)$cluster <- clusters
# V(graph)$source <- ifelse(grepl("wendrich_", integrated$Sample),
#                           "Wendrich et al.",
#                           ifelse(grepl("shahan_", integrated$Sample), "Shahan et al.",
#                                  ifelse(grepl("denyer_", integrated$Sample), "Denyer et al.",
#                                         "Ours")))

# subset the graph with cluster 3 vertices
subgraph <- induced_subgraph(graph, V(graph)[V(graph)$cluster == 3])
V(subgraph)$cluster <- factor(igraph::cluster_louvain(subgraph)$membership)

tibble(cluster = V(subgraph)$cluster,
       source = V(subgraph)$source) %>%
  count(source, cluster) %>%
  group_by(cluster) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(pct, factor(cluster), fill = source)) +
  geom_col() +
  scale_fill_brewer(palette = "Dark2")

# # visualise (not very informative)
# temp <- delete_edges(subgraph, E(subgraph)[E(subgraph)$weight < 0.3])
# plot(temp,
#      edge.color = "grey",
#      vertex.color = c("seagreen3", "steelblue")[(V(temp)$source == "Ours") + 1],
#      vertex.label = NA,
#      vertex.size = 5)

# check how many edges occur across data sources
graph_ends <- ends(subgraph, E(subgraph))
graph_ends <- tibble(
  end1 = V(subgraph)$source[graph_ends[, 1]],
  end2 = V(subgraph)$source[graph_ends[, 2]],
  weight = E(subgraph)$weight
) %>%
  mutate(combo = paste(end1, end2, sep = "--"))

# edge weight is slightly higher for edges between same-data than across data
graph_ends %>%
  mutate(same = end1 == end2) %>%
  ggplot(aes(weight)) +
  geom_density(aes(colour = same, group = combo),
               fill = "lightgrey", alpha = 0.2, size = 1) +
  scale_x_log10() +
  scale_colour_brewer(palette = "Dark2")

# how many edges with each combo
graph_ends %>%
  mutate(same = end1 == end2) %>%
  count(same, combo) %>%
  mutate(pct = n/sum(n)*100)

# randomize the labels to get a null expectation
get_combos <- function(x){
  graph_ends <- ends(x, E(x))
  graph_ends <- tibble(
    end1 = V(x)$source[graph_ends[, 1]],
    end2 = V(x)$source[graph_ends[, 2]],
  )
  # sum(graph_ends$end1 != graph_ends$end2)/nrow(graph_ends)
  combo <- paste(pmin(graph_ends$end1, graph_ends$end2),
                 pmax(graph_ends$end1, graph_ends$end2),
                 sep = "--")
  table(combo)/nrow(graph_ends)*100
}

iter <- vector("list", 100)
temp <- subgraph
for(i in 1:length(iter)){
  V(temp)$source <- sample(V(temp)$source)
  iter[[i]] <- get_combos(temp)
}
shuffled_combos <- map_dfr(iter, enframe, .id = "iter")
shuffled_combos %>%
  group_by(name) %>%
  summarise(lo = quantile(value, 0.025),
            mid = median(value),
            hi = quantile(value, 0.975))

full_join(
  shuffled_combos %>%
    group_by(name) %>%
    summarise(lo = quantile(value, 0.025),
              mid = median(value),
              hi = quantile(value, 0.975)),
  subgraph %>%
    get_combos() %>%
    enframe(value = "obs")
) %>%
  ggplot(aes(obs, name)) +
  geom_segment(aes(x = 0, xend = obs, yend = name),
               colour = "darkgrey") +
  geom_point(aes(fill = "Observed")) +
  geom_linerange(aes(xmin = lo, xmax = hi, linetype = "95% CI"),
                 size = 5, alpha = 0.2) +
  labs(x = "% edge connectons", y = "",
       subtitle = "Cell-cell edges between datasets in cluster 3") +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.title = element_blank())

subgraph %>%
  get_combos() %>%
  enframe(value = "obs") %>%
  full_join(shuffled_combos) %>%
  replace_na(list(obs = 0)) %>%
  mutate(dif = obs - value) %>%
  ggplot(aes(dif, name)) +
  ggridges::geom_density_ridges(aes(height = ..density..), stat = "density") +
  geom_vline(xintercept = 0)


# another way to do it would have been to generate a network with the observed degree
# distribution: degree.sequence.game()
# but it's slower than shuffling labels


# we could randomize edges with probability based on the edge weights
# but this doesn't preserve the degree, so I'm not sure that would be a good approach
# in any case it's not hugely different from the observed proportions
plot(density(degree(subgraph)), ylim = c(0, 0.02), col = "orange")
purrr::rerun(10,
             subgraph %>%
               rewire(each_edge(1 - E(subgraph)$weight)) %>%
               degree()) %>%
  map(density) %>%
  map(lines)

# degree-preserving randomization - not sure this is correct
iter <- purrr::rerun(100,
                     subgraph %>%
                       rewire(keeping_degseq(loops = FALSE, niter = vcount(subgraph)*20)) %>%
                       get_combos()) %>%
  bind_rows()

iter %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(lo = quantile(value, 0.025)*100, hi = quantile(value, 0.975)*100)


# check for the whole graph
graph_ends <- ends(graph, E(graph))
graph_ends <- tibble(
  end1 = V(graph)$source[graph_ends[, 1]],
  end2 = V(graph)$source[graph_ends[, 2]],
  weight = E(graph)$weight
) %>%
  mutate(combo = paste(end1, end2, sep = "--"))

# edge weight is slightly higher for edges between same-data than across data
graph_ends %>%
  ggplot(aes(weight)) +
  geom_density(aes(colour = combo), fill = "lightgrey", alpha = 0.2, size = 1) +
  scale_x_log10() +
  scale_colour_brewer(palette = "Dark2")



# Early cells -------------------------------------------------------------

# subset the graph with cluster 3 vertices
subgraph <- induced_subgraph(graph, V(graph)[V(graph)$cluster == 3])

# re-cluster to get finer patterns
V(subgraph)$cluster <- factor(igraph::cluster_walktrap(subgraph)$membership)

# subset single cell object with those cells only
clust3 <- integrated[, V(subgraph)$cell]
clust3$cluster <- factor(V(subgraph)$cluster)
clust3 <- clust3[which(rowSums(logcounts(clust3)) > 0), ]

# DOF1.5
wrap_plots(
  umap_highlight_cluster(integrated, 3) + ggtitle("Cluster 3"),
  umap_highlight_gene(clust3, "AT1G29160"), # DOF1.5
  umap_highlight_gene(clust3, "AT5G52870") # MAKR5
)


clust3 %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = cluster), size = 0.5) +
  ggthemes::scale_colour_tableau("Tableau 20")

wrap_plots(
  scater::plotReducedDim(clust3, "UMAP30_MNN_logvst",
                         colour_by = "cluster", text_by = "cluster"),
  scater::plotExpression(clust3, "AT1G29160", x = "cluster"),
  scater::plotExpression(clust3, "AT5G52870", x = "cluster")
)

# get LFC across these clusters
clust3_markers <- scran::findMarkers(clust3,
                                     groups = clust3$cluster,
                                     pval.type = "some",
                                     min.prop = 4/8)
clust3_markers <- clust3_markers %>%
  as.list() %>%
  map_dfr(as_tibble, rownames = "id", .id = "cluster")

# create matrix of pairwise LFC
clust3_markers_mat <- clust3_markers %>%
  # filter(FDR < 0.01) %>%
  select(cluster, id, starts_with("logFC")) %>%
  pivot_longer(starts_with("logFC")) %>%
  mutate(comparison = paste(name, cluster, sep = "_")) %>%
  select(id, comparison, value) %>%
  drop_na() %>%
  pivot_wider(names_from = "comparison", values_from = "value") %>%
  column_to_rownames("id") %>%
  as.matrix()

# LFC correlations
clust3_cor <- cor(t(clust3_markers_mat),
    t(clust3_markers_mat[c("AT1G29160", "AT5G52870"), ]))

# based on ranks for both genes
ranks <- rank(-clust3_cor[, "AT1G29160"])^2 + rank(-clust3_cor[, "AT5G52870"])^2
gois <- names(sort(ranks)[1:10])

png("~/temp/early_dof1.5_makr5.png", width = 1200, height = 1200, res = 300)
wrap_plots(
  map(c("AT1G29160", "AT5G52870", gois),
      ~ umap_highlight_gene(integrated, .x))
)
dev.off()


# based on DOF1.5 only
ranks <- rank(-clust3_cor[, "AT1G29160"])
gois <- names(sort(ranks)[1:10])

png("~/temp/early_dof1.5.png", height = 1200, width = 1200, res = 300)
wrap_plots(
  map(c("AT1G29160", "AT5G52870", gois[1:10]),
      ~ umap_highlight_gene(integrated, .x))
)
dev.off()


# based on MARK5 only
ranks <- rank(-clust3_cor[, "AT5G52870"])
gois <- names(sort(ranks)[1:10])

png("~/temp/early_makr5.png", height = 1200, width = 1200, res = 300)
wrap_plots(
  map(c("AT1G29160", "AT5G52870", gois[1:10]),
      ~ umap_highlight_gene(integrated, .x))
)
dev.off()

# write the correlations to a file
clust3_cor %>%
  as_tibble(rownames = "gene") %>%
  filter(AT1G29160 > 0.5 | AT5G52870 > 0.5) %>%
  write_csv("~/temp/clust3_cor.csv")
