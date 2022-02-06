library(SingleCellExperiment)
library(tidyverse)
library(ggpointdensity)
library(patchwork)

# set ggplot2 theme
theme_set(theme_classic() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")


# Read data ---------------------------------------------------------------

leaf <- readRDS("data/processed/SingleCellExperiment/leaf_batches_strictfilt.rds")

leaf$source <- ifelse(grepl("kim_", leaf$Sample),
                      "Leaf\n(Kim et al.)", "Phloem Pole Atlas")

# list of curated genes
curated_genes_leaf <- read_csv("data/raw/signature_genes_leaf.csv") %>%
  select(tissue = expression_domain, id = gene_id, name = gene_name, set) %>%
  mutate(across(everything(), str_trim)) %>%
  mutate(name = ifelse(is.na(name), id, name)) %>%
  mutate(tissue = str_replace_all(tissue, " ", "\n"),
         id = toupper(id))

# usual list of curated genes
curated_genes <- read_csv("data/raw/signature_genes_circle_plot.csv") %>%
  select(tissue = Expression, id = gene_ID, name = gene_name) %>%
  mutate(tissue = str_replace(tissue, " \\+ ", "\n&\n")) %>%
  mutate(tissue = str_replace_all(tissue, " ", "\n")) %>%
  mutate(id = toupper(str_trim(id)))

# load our phloem ring data
ring <- readRDS("data/processed/SingleCellExperiment/ring_batches_strictfilt.rds")

# add ring clusters to leaf integration
leaf$cluster_ring <- NA
leaf$cluster_ring[match(colnames(ring), colnames(leaf))] <- ring$cluster_mnn_logvst

# quickly check the matching was done correctly
check1 <- leaf$cluster_ring; names(check1) <- colnames(leaf)
check2 <- ring$cluster_mnn_logvst; names(check2) <- colnames(ring)
check1 <- check1[names(check2)]
all(check1 == check2)

# EDA ---------------------------------------------------------------------

scater::plotReducedDim(leaf, "UMAP30_MNN_logcounts",
                       colour_by = "cluster_mnn_logcounts",
                       text_by = "cluster_mnn_logcounts")
scater::plotReducedDim(leaf, "UMAP30_MNN_logvst",
                       colour_by = "cluster_mnn_logvst",
                       text_by = "cluster_mnn_logvst")
leaf %>%
  colData() %>%
  as_tibble() %>%
  ggplot(aes(cluster_mnn_logvst, cluster_mnn_logcounts)) +
  geom_count()

# UMAP split by cluster
leaf %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(Sample = str_remove(Sample, "MAKR5_")) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(leaf, "UMAP30_MNN_logvst") %>%
               select(-cluster_mnn_logvst),
             colour = "lightgrey") +
  geom_point(aes(colour = source), alpha = 0.1) +
  facet_wrap(~ cluster_mnn_logvst)


# UMAP with samples split by origin
leaf %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(Sample = str_remove(Sample, "MAKR5_")) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(leaf, "UMAP30_MNN_logvst") %>% select(-source),
             colour = "lightgrey") +
  geom_point(aes(colour = source), alpha = 0.1) +
  geom_label(data = getReducedDim(leaf, "UMAP30_MNN_logvst") %>% select(-source),
             stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  scale_colour_manual(values = c("seagreen3", "black")) +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "source") +
  coord_equal() + theme_void() +
  theme(panel.spacing = unit(1, "in")) +
  guides(colour = "none") +
  facet_grid(~ source)

# umap with both of them
p1 <- leaf %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  arrange(desc(source)) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = source), alpha = 0.3, size = 0.5) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines"), size = 3) +
  scale_colour_manual(values = c("seagreen3", "black")) +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "source") +
  coord_equal() + theme_void() +
  guides(colour = "none")

pdf("documents/pdf for figures/integrated_leaf/umap_leaf_root.pdf",
    width = 3, height = 3)
p1
dev.off()


# cell proportion in each cluster
p2 <- leaf %>%
  colData() %>%
  as_tibble() %>%
  count(cluster_mnn_logvst, source) %>%
  group_by(cluster_mnn_logvst) %>%
  mutate(pct = n/sum(n)*100,
         ncluster = sum(n)) %>%
  ungroup() %>%
  mutate(cluster = paste0(cluster_mnn_logvst, " (", ncluster, ")")) %>%
  mutate(cluster = fct_reorder(cluster, as.numeric(cluster_mnn_logvst))) %>%
  ggplot(aes(cluster, pct)) +
  geom_col(aes(fill = source)) +
  scale_fill_manual(values = c("seagreen3", "black")) +
  labs(x = "Cluster", y = "% cells")

pdf("documents/pdf for figures/integrated_leaf/pct_cells_per_cluster.pdf",
    width = 6, height = 6)
p2 + coord_flip() + theme_minimal() +
  theme(panel.grid = element_blank())
dev.off()


p3 <- leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = curated_genes_leaf$id, melted = TRUE,
                exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(curated_genes_leaf, by = "id") %>%
  # scale expression for each gene
  group_by(id) %>%
  mutate(expr_center = (expr - mean(expr))/sd(expr)) %>%
  # calculate percentage of cells where the gene is detected per cluster
  group_by(cluster_mnn_logvst, id, name, tissue) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  # plot
  ggplot(aes(name, factor(cluster_mnn_logvst))) +
  geom_point(aes(colour = mean_expr, size = pct_expressing)) +
  facet_grid(~ tissue, scales = "free", space = "free") +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "#f7f7f7") +
  scale_size_continuous(limits = c(0, 100)) +
  guides(shape = "none") +
  labs(y = "Cluster", x = "", size = "% Detected",
       colour = "Z-score\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p3


# Curated genes -----------------------------------------------------------

# set "a" - we actually will not show this
leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = curated_genes_leaf$id, melted = TRUE,
                exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(curated_genes_leaf, by = "id") %>%
  filter(set == "a") %>%
  # scale expression for each gene
  group_by(id) %>%
  mutate(expr_center = (expr - mean(expr))/sd(expr)) %>%
  # calculate percentage of cells where the gene is detected per cluster
  group_by(cluster_mnn_logvst, id, name, tissue) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  # clip the scaled expression
  mutate(mean_expr = ifelse(abs(mean_expr) > 4, sign(mean_expr)*4, mean_expr)) %>%
  # plot
  ggplot(aes(name, factor(cluster_mnn_logvst))) +
  geom_point(aes(colour = mean_expr, size = pct_expressing)) +
  facet_grid(~ tissue, scales = "free", space = "free") +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "gray94",
                         limits = c(-4, 4), breaks = seq(-4, 4, by = 1),
                         labels = c("<= -4", -3, -2, -1, 0, 1, 2, 3, ">= 4")) +
  scale_size_continuous(limits = c(0, 100)) +
  guides(shape = "none") +
  labs(y = "Cluster", x = "", size = "% Detected",
       colour = "Z-score\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# usual genes
p1 <- leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = curated_genes$id, melted = TRUE,
                exprs_values = "logvst", zeros_to_na = FALSE) %>%
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
  # clip the scaled expression
  mutate(mean_expr = ifelse(abs(mean_expr) > 4, sign(mean_expr)*4, mean_expr)) %>%
  # plot
  ggplot(aes(name, factor(cluster_mnn_logvst))) +
  geom_point(aes(colour = mean_expr, size = pct_expressing)) +
  facet_grid(~ tissue, scales = "free", space = "free") +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "gray94",
                         limits = c(-4, 4), breaks = seq(-4, 4, by = 1),
                         labels = c("<= -4", -3, -2, -1, 0, 1, 2, 3, ">= 4")) +
  scale_size_continuous(limits = c(0, 100)) +
  guides(shape = "none") +
  labs(y = "Cluster", x = "", size = "% Detected",
       colour = "Z-score\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# set "b"
p2 <- leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = curated_genes_leaf$id, melted = TRUE,
                exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(curated_genes_leaf, by = "id") %>%
  filter(set == "b") %>%
  # scale expression for each gene
  group_by(id) %>%
  mutate(expr_center = (expr - mean(expr))/sd(expr)) %>%
  # calculate percentage of cells where the gene is detected per cluster
  group_by(cluster_mnn_logvst, id, name, tissue) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  # clip the scaled expression
  mutate(mean_expr = ifelse(abs(mean_expr) > 4, sign(mean_expr)*4, mean_expr)) %>%
  # plot
  ggplot(aes(name, factor(cluster_mnn_logvst))) +
  geom_point(aes(colour = mean_expr, size = pct_expressing)) +
  facet_grid(~ tissue, scales = "free", space = "free") +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "gray94",
                         limits = c(-4, 4), breaks = seq(-4, 4, by = 1),
                         labels = c("<= -4", -3, -2, -1, 0, 1, 2, 3, ">= 4")) +
  scale_size_continuous(limits = c(0, 100)) +
  guides(shape = "none") +
  labs(y = "Cluster", x = "", size = "% Detected",
       colour = "Z-score\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf("documents/pdf for figures/integrated_leaf/cell_type_annotation_leaf.pdf",
    width = 5.7, height = 4.25*2)
(p1 / p2) + plot_layout(guides = "collect")
dev.off()


# plot them all combined
curated_genes_combined <- bind_rows(
  curated_genes_leaf %>%
    filter(set == "b") %>%
    select(-set),
  curated_genes %>%
    filter(name %in% c("S17", "CALS8", "SAPL"))
)

pdf("documents/pdf for figures/integrated_leaf/cell_type_annotation_leaf_combined.pdf",
    width = 5.7, height = 4.3)
leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = curated_genes_combined$id, melted = TRUE,
                exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(curated_genes_combined, by = "id") %>%
  mutate(tissue = factor(tissue,
                         levels = c("BUNDLE\nSHEATH",
                                    "EPIDERMIS",
                                    "MESOPHYLL",
                                    "PP",
                                    "PPP",
                                    "CC"))) %>%
  # scale expression for each gene
  group_by(id) %>%
  mutate(expr_center = (expr - mean(expr))/sd(expr)) %>%
  # calculate percentage of cells where the gene is detected per cluster
  group_by(cluster_mnn_logvst, id, name, tissue) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  # clip the scaled expression
  mutate(mean_expr = ifelse(abs(mean_expr) > 4, sign(mean_expr)*4, mean_expr)) %>%
  # plot
  ggplot(aes(name, factor(cluster_mnn_logvst))) +
  geom_point(aes(colour = mean_expr, size = pct_expressing)) +
  facet_grid(~ tissue, scales = "free", space = "free") +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "gray94",
                         limits = c(-4, 4), breaks = seq(-4, 4, by = 1),
                         labels = c("<= -4", -3, -2, -1, 0, 1, 2, 3, ">= 4")) +
  scale_size_continuous(limits = c(0, 100)) +
  guides(shape = "none") +
  labs(y = "Cluster", x = "", size = "% Detected",
       colour = "Z-score\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8))
dev.off()



# Violin plots of expression ----------------------------------------------

# genes to highlight distribution
genes_to_plot <- curated_genes_combined %>%
  filter(tissue %in% c("PP", "PPP")) %>%
  filter(id %in% rownames(leaf)) %>%
  pull(id)

# get counts for annotating the plot
gene_cts <- leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = genes_to_plot, melted = TRUE,
                exprs_values = "logcounts", zeros_to_na = TRUE) %>%
  left_join(curated_genes_combined, by = "id") %>%
  filter(cluster_mnn_logvst == 6) %>%
  group_by(id, name, source) %>%
  summarise(n = n(), nexpr = sum(!is.na(expr)), y = max(expr, na.rm = TRUE) + 1) %>%
  ungroup() %>%
  mutate(label = paste0(nexpr, "\n(",
                        ifelse(round(nexpr/n*100) == 0,
                               round(nexpr/n*100, 1),
                               round(nexpr/n*100)), "%)"))

pdf("documents/pdf for figures/integrated_leaf/expression_violins.pdf",
    width = 5.7, height = 4.3)
leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = genes_to_plot, melted = TRUE,
                exprs_values = "logcounts", zeros_to_na = TRUE) %>%
  drop_na(expr) %>%
  filter(cluster_mnn_logvst == 6) %>%
  left_join(curated_genes_combined, by = "id") %>%
  ggplot(aes(source, expr)) +
  geom_violin(scale = "width") +
  geom_text(data = gene_cts, aes(x = source, y = -3, label = label),
            vjust = 0) +
  ggbeeswarm::geom_quasirandom(aes(colour = source), size = 0.8) +
  facet_wrap(~ name, nrow = 2, scales = "free_x") +
  scale_fill_manual(values = c("seagreen3", "black")) +
  scale_colour_manual(values = c("seagreen3", "black")) +
  scale_y_continuous(breaks = seq(0, 9, by = 3)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank()) +
  labs(x = "", y = "logcounts")
dev.off()


# XPP ---------------------------------------------------------------------
source("analysis/functions/plots.R")

# a couple of genes that mark the XPP
genes_to_plot <- curated_genes_leaf %>%
  filter(tissue %in% c("XPP")) %>%
  filter(id %in% rownames(leaf)) %>%
  pull(id)

gene_cts <- leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = genes_to_plot, melted = TRUE,
                exprs_values = "logcounts", zeros_to_na = TRUE) %>%
  left_join(curated_genes_leaf, by = "id") %>%
  filter(cluster_mnn_logvst == 6) %>%
  group_by(id, name, source) %>%
  summarise(n = n(), nexpr = sum(!is.na(expr)), y = max(expr, na.rm = TRUE) + 1) %>%
  ungroup() %>%
  mutate(label = paste0(nexpr, "\n(",
                        ifelse(nexpr/n*100 < 1,
                               round(nexpr/n*100, 1),
                               round(nexpr/n*100)), "%)"))

p1 <- leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = genes_to_plot, melted = TRUE,
                exprs_values = "logcounts", zeros_to_na = TRUE) %>%
  drop_na(expr) %>%
  filter(cluster_mnn_logvst == 6) %>%
  left_join(curated_genes_leaf, by = "id") %>%
  ggplot(aes(source, expr)) +
  geom_violin(scale = "width") +
  geom_text(data = gene_cts, aes(x = source, y = -3, label = label),
            vjust = 0) +
  ggbeeswarm::geom_quasirandom(aes(colour = source), size = 0.8) +
  facet_grid(~ name) +
  scale_fill_manual(values = c("seagreen3", "black")) +
  scale_colour_manual(values = c("seagreen3", "black")) +
  scale_y_continuous(breaks = seq(0, 9, by = 3)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank()) +
  labs(x = "", y = "logcounts",
       subtitle = "Expression in Cluster 6")

p2 <- umap_highlight_gene(leaf, "AT1G02460")
p3 <- umap_highlight_gene(leaf, "AT4G30450")
p4 <- umap_highlight_gene(leaf, "AT2G36120")

wrap_plots(p1, p2, p3, p4)

pdf("documents/pdf for figures/integrated_leaf/expression_violins_xpp.pdf",
    width = 5.7, height = 2.5)
p1
dev.off()


# Percent cells expressing these genes in each dataset and each cluster
leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = genes_to_plot, melted = TRUE,
                exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(curated_genes_leaf, by = "id") %>%
  # calculate percentage of cells where the gene is detected per cluster
  group_by(source, cluster_mnn_logvst, id, name) %>%
  summarise(n_expressing = sum(expr > 0), n_cells = n()) %>%
  group_by(cluster_mnn_logvst, id, name) %>%
  mutate(pct_expressing = n_expressing/sum(n_cells)*100) %>%
  ungroup() %>%
  # plot
  ggplot(aes(pct_expressing, factor(cluster_mnn_logvst))) +
  geom_col(aes(fill = source)) +
  facet_grid(~ name) +
  scale_fill_manual(values = c("seagreen3", "black")) +
  theme_minimal() +
  labs(x = "% Cells with detected expression (per cluster)",
       y = "Cluster (LEAF + RING)")


# cluster equivalence
leaf %>%
  colData() %>%
  as_tibble() %>%
  mutate(highlight = cluster_ring == 11) %>%
  ggplot(aes(factor(cluster_mnn_logvst), factor(cluster_ring))) +
  geom_count() +
  theme_minimal() +
  labs(x = "Cluster Leaf + Ring integrated", y = "Cluster ring only data")



# SWEET genes -------------------------------------------------------------

sweets <- c(SWEET2 = "AT3G14770",
            SWEET3 = "AT5G53190",
            SWEET4 = "AT3G28007",
            SWEET5 = "AT5G62850",
            SWEET6 = "AT1G66770",
            SWEET7 = "AT4G10850",
            SWEET8 = "AT5G40260",
            SWEET9 = "AT2G39060",
            SWEET10 = "AT5G50790",
            SWEET11 = "AT3G48740",
            SWEET12 = "AT5G23660",
            SWEET13 = "AT5G50800",
            SWEET14 = "AT4G25010",
            SWEET15 = "AT5G13170",
            SWEET16 = "AT3G16690",
            SWEET17 = "AT4G15920") %>%
  enframe("name", "id")

# annotation based on sweets
leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = sweets$id, melted = TRUE,
                exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(sweets, by = "id") %>%
  # scale expression for each gene
  group_by(id) %>%
  mutate(expr_center = (expr - mean(expr))/sd(expr)) %>%
  # calculate percentage of cells where the gene is detected per cluster
  group_by(cluster_mnn_logvst, id, name) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  # plot
  ggplot(aes(name, factor(cluster_mnn_logvst))) +
  geom_point(aes(colour = mean_expr, size = pct_expressing)) +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "#f7f7f7") +
  scale_size_continuous(limits = c(0, 100)) +
  guides(shape = "none") +
  labs(y = "Cluster", x = "", size = "% Detected",
       colour = "Z-score\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# annotate separately for each source
# annotation based on sweets
# this gives a distorted perspective because of different numbers of cells in each dataset
leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = sweets$id, melted = TRUE,
                exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(sweets, by = "id") %>%
  # scale expression for each gene
  group_by(id, source) %>%
  mutate(expr_center = (expr - mean(expr))/sd(expr)) %>%
  # calculate percentage of cells where the gene is detected per cluster
  group_by(cluster_mnn_logvst, id, name, source) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  drop_na() %>%
  # plot
  ggplot(aes(source, cluster_mnn_logvst)) +
  geom_point(aes(colour = mean_expr, size = pct_expressing, shape = source)) +
  facet_wrap(~ name, nrow = 1, strip.position = "bottom") +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "grey") +
  scale_size_continuous(limits = c(0, 100)) +
  labs(y = "Cluster", x = "", size = "% Detected",
       colour = "Z-score\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = sweets$id, melted = TRUE,
                exprs_values = "logvst") %>%
  left_join(sweets, by = "id") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  # geom_label(stat = "centroid",
  #            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
  #            alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ name) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression") +
  theme_void()



# Galactinol --------------------------------------------------------------

galactinol <- c(
  GOLS1 = "AT2G47180",
  GOLS2 = "AT1G56600",
  GOLS3 = "AT1G09350",
  GOLS4 = "AT1G60470",
  GOLS5 = "AT5G23790",
  GOLS6 = "AT4G26250",
  GOLS7 = "AT1G60450",
  GOLS8 = "AT3G28340",
  GOLS9 = "AT3G06260",
  GOLS10 = "AT5G30500"
) %>%
  enframe("name", "id")

# annotation based on galactinol
leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = galactinol$id, melted = TRUE,
                exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(galactinol, by = "id") %>%
  # scale expression for each gene
  group_by(id) %>%
  mutate(expr_center = (expr - mean(expr))/sd(expr)) %>%
  # calculate percentage of cells where the gene is detected per cluster
  group_by(cluster_mnn_logvst, id, name) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  # plot
  ggplot(aes(name, factor(cluster_mnn_logvst))) +
  geom_point(aes(colour = mean_expr, size = pct_expressing)) +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "#f7f7f7") +
  scale_size_continuous(limits = c(0, 100)) +
  guides(shape = "none") +
  labs(y = "Cluster", x = "", size = "% Detected",
       colour = "Z-score\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# painted on the UMAP
leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = galactinol$id, melted = TRUE,
                exprs_values = "logcounts", zeros_to_na = TRUE) %>%
  left_join(galactinol, by = "id") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted), size = 0.3) +
  facet_wrap(~ name) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression") +
  theme_void() + coord_equal()

# annotate separately for each source
# annotation based on galactinol
leaf %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = galactinol$id, melted = TRUE,
                exprs_values = "logvst", zeros_to_na = FALSE) %>%
  left_join(galactinol, by = "id") %>%
  # scale expression for each gene
  group_by(id, source) %>%
  mutate(expr_center = (expr - mean(expr))/sd(expr)) %>%
  # calculate percentage of cells where the gene is detected per cluster
  group_by(cluster_mnn_logvst, id, name, source) %>%
  mutate(pct_expressing = (sum(expr > 0)/n())*100) %>%
  summarise(mean_expr = mean(expr_center),
            pct_expressing = unique(pct_expressing)) %>%
  ungroup() %>%
  drop_na() %>%
  # plot
  ggplot(aes(source, cluster_mnn_logvst)) +
  geom_point(aes(colour = mean_expr, size = pct_expressing, shape = source)) +
  facet_wrap(~ name, nrow = 1, strip.position = "bottom") +
  scale_colour_gradient2(low = "#313695", high = "#a50026", mid = "grey") +
  scale_size_continuous(limits = c(0, 100)) +
  labs(y = "Cluster", x = "", size = "% Detected",
       colour = "Z-score\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_blank())


# Check cell mixing -------------------------------------------------------
library(igraph)

# build shared neighbours graph
# this is the same graph used to generate the clusters used above
graph <- scran::buildSNNGraph(leaf, k = 50, type = "jaccard",
                              use.dimred = "MNN_logvst")

# get those clusters (and confirm they are indeed the same)
clusters <- factor(igraph::cluster_louvain(graph)$membership)
all.equal(clusters, leaf$cluster_mnn_logvst)

# add metadata to graph
V(graph)$cell <- colnames(leaf)[as.vector(V(graph))]
V(graph)$cluster <- clusters
V(graph)$source <- ifelse(grepl("Root", leaf$source), "Root", "Leaf")

# subset the graph with cluster 6 vertices
clust6_graph <- induced_subgraph(graph, V(graph)[clusters == 6])

# visualise (not very informative)
temp <- delete_edges(clust6_graph, E(clust6_graph)[E(clust6_graph)$weight < 0.3])
plot(temp,
     edge.color = "grey",
     vertex.color = c("seagreen3", "black")[(V(temp)$source == "Leaf") + 1],
     vertex.label = NA,
     vertex.size = 5)

# check how many edges occur across data sources
graph_ends <- ends(clust6_graph, E(clust6_graph))
graph_ends <- tibble(
  end1 = V(clust6_graph)$source[graph_ends[, 1]],
  end2 = V(clust6_graph)$source[graph_ends[, 2]],
  weight = E(clust6_graph)$weight
) %>%
  mutate(combo = paste(end1, end2, sep = "--"))

# edge weight is slightly higher for edges between same-data than across data
graph_ends %>%
  ggplot(aes(weight)) +
  geom_density(aes(colour = combo), fill = "lightgrey", alpha = 0.2, size = 1) +
  scale_x_log10() +
  scale_colour_brewer(palette = "Dark2")

# how many edges with each combo
graph_ends %>%
  count(combo) %>%
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

iter <- vector("numeric", 1000)
temp <- clust6_graph
for(i in 1:length(iter)){
  V(temp)$source <- sample(V(temp)$source)
  iter[i] <- get_combos(temp)
}
quantile(iter, c(0.025, 0.5, 0.975))


### new approach
iter <- vector("list", 100)
temp <- clust6_graph
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
  clust6_graph %>%
    get_combos() %>%
    enframe(value = "obs")
) %>%
  mutate(obs = ifelse(is.na(obs), 0, obs)) %>%
  ggplot(aes(obs, name)) +
  geom_segment(aes(x = 0, xend = obs, yend = name), colour = "darkgrey") +
  geom_point() +
  geom_linerange(aes(xmin = lo, xmax = hi), size = 5, alpha = 0.2)

clust6_graph %>%
  get_combos() %>%
  enframe(value = "obs") %>%
  full_join(shuffled_combos) %>%
  replace_na(list(obs = 0)) %>%
  mutate(dif = (obs - value)) %>%
  ggplot(aes(name, dif)) +
  geom_violin(scale = "width") +
  geom_vline(xintercept = 0)


# another way to do it would have been to generate a network with the observed degree
# distribution: degree.sequence.game()
# but it's slower than shuffling labels


# we could randomize edges with probability based on the edge weights
# but this doesn't preserve the degree, so I'm not sure that would be a good approach
# in any case it's not hugely different from the observed proportions
plot(density(degree(clust6_graph)), ylim = c(0, 0.02), col = "orange")
purrr::rerun(10,
             clust6_graph %>%
               rewire(each_edge(1 - E(clust6_graph)$weight)) %>%
               degree()) %>%
  map(density) %>%
  map(lines)

# degree-preserving randomization - not sure this is correct
iter <- purrr::rerun(100,
                     clust6_graph %>%
                       rewire(keeping_degseq(loops = FALSE, niter = vcount(clust6_graph)*20)) %>%
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
