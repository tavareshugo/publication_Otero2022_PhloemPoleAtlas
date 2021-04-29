library(SingleCellExperiment)
library(tidyverse)
library(ggpointdensity)
library(patchwork)

# set ggplot2 theme
theme_set(theme_classic() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")


# Read data ---------------------------------------------------------------

wendrich <- readRDS("data/processed/SingleCellExperiment/wendrich_batches_strictfilt.rds")
all <- readRDS("data/processed/SingleCellExperiment/integrated_batches_strictfilt_leaner.rds")

wendrich$source <- ifelse(grepl("wendrich_", wendrich$Sample),
                          "Wendrich et al.", "Ours")

# new list of curated genes (provided by Sofia)
curated_genes <- read_csv("data/raw/signature_genes_circle_plot.csv") %>%
  select(tissue = Expression, id = gene_ID, name = gene_name) %>%
  mutate(tissue = str_replace(tissue, " \\+ ", "\n&\n")) %>%
  mutate(tissue = str_replace_all(tissue, " ", "\n")) %>%
  mutate(id = toupper(str_trim(id)))


# EDA ---------------------------------------------------------------------

scater::plotReducedDim(wendrich, "UMAP30_MNN_logcounts",
                       colour_by = "cluster_mnn_logcounts",
                       text_by = "cluster_mnn_logcounts")
scater::plotReducedDim(wendrich, "UMAP30_MNN_logvst",
                       colour_by = "cluster_mnn_logvst",
                       text_by = "cluster_mnn_logvst")

wendrich %>%
  colData() %>%
  as_tibble() %>%
  ggplot(aes(cluster_mnn_logvst, cluster_mnn_logcounts)) +
  geom_count()

# UMAP with samples split by origin
wendrich %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(Sample = str_remove(Sample, "MAKR5_")) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(wendrich, "UMAP30_MNN_logvst") %>% select(-source),
             colour = "lightgrey") +
  geom_point(aes(colour = source), alpha = 0.1) +
  geom_label(data = getReducedDim(wendrich, "UMAP30_MNN_logvst") %>% select(-source),
             stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  scale_colour_manual(values = c("seagreen3", "steelblue")) +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "source") +
  coord_equal() + theme_void() +
  theme(panel.spacing = unit(1, "in")) +
  guides(colour = "none") +
  facet_grid(~ source)


wendrich %>%
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
  # tidy up the tissue column for nicer plotting
  filter(tissue != "early") %>%
  mutate(tissue = factor(tissue,
                         levels = c("CELL\nCYCLE\nGENES",
                                    "CC", "CC\n&\nMSE",
                                    "PPP\n&\nCC", "PPP",
                                    "EARLY\nPSE",
                                    "LATE\nPSE", "OUTER\nLAYERS"))) %>%
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

# checking dividing cells
scater::plotReducedDim(wendrich, "UMAP30_MNN_logvst",
                       colour_by = "AT4G32830",
                       text_by = "cluster_mnn_logvst")
