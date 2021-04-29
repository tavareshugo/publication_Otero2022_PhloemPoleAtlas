library(SingleCellExperiment)
library(tidyverse)
library(ggpointdensity)
library(patchwork)

# set ggplot2 theme
theme_set(theme_classic() + theme(text = element_text(size = 14)))

# source util functions
source("analysis/functions/utils.R")



# Read Data ---------------------------------------------------------------

# hard filtering
ring <- readRDS("data/processed/SingleCellExperiment/ring_batches_strictfilt.rds")

# order clusters more logically
ring$cluster_mnn_logvst <- factor(ring$cluster_mnn_logvst,
                              levels = c("8", "10", "9",
                                         "5", "3", "1",
                                         "7", "4", "14", "11",
                                         "13",
                                         "12", "2", "6",
                                         "15"))

# cluster-by-cluster test
cluster_test <- read.csv("data/processed/gene_sets/ring_strictfilt_cluster_markers.csv",
                         stringsAsFactors = FALSE) %>%
  mutate(cluster = factor(cluster))

# curated list of genes
curated_genes <- list(`Outer Layers` = tribble(~id, ~name,
                                               "AT1G79580", "SMB",
                                               "AT1G79840", "GL2",
                                               "AT5G14750", "WER"),
                      `SE (early)` = tribble(~id, ~name,
                                             "AT1G05470", "CVP2",
                                             "AT1G54330", "NAC020",
                                             "AT2G37590", "PEAR1",
                                             "AT5G02460", "PEAR2"),
                      `SE (late)` = tribble(~id, ~name,
                                            "AT1G06490", "CALS7",
                                            "AT5G17260", "NAC086"),
                      `PPP` = tribble(~id, ~name,
                                      "AT2G22850", "S17",
                                      "AT3G14570", "CALS8"),
                      `CC` = tribble(~id, ~name,
                                     "AT2G38640", "AT2G38640",
                                     "AT3G12730", "SAPL",
                                     "AT1G22710", "SUC2",
                                     "AT5G02600", "NAKR1",
                                     "AT5G57350", "AHA3"),
                      `early` = tribble(~id, ~name,
                                   "AT1G29160", "DOF1.5",
                                   "AT2G36400", "GRF3",
                                   "AT5G52870", "MAKR5")) %>%
  bind_rows(.id = "tissue")

# cyclins
cyclins <- rowData(ring)[grep("CYC[A,B,D]",
                              rowData(ring)$alternative_name,
                              ignore.case = TRUE), ]
cyclins <- as_tibble(cyclins) # convert to tbl

# clean up the names a bit
cyclins <- cyclins %>%
  mutate(alternative_name = gsub(".*/", "", alternative_name)) %>%
  mutate(cyclin = gsub("^CYC|[0-9].*$", "", alternative_name))

# add only a few as example to curated genes
curated_genes <- cyclins %>%
  filter(alternative_name %in% c("CYCB1-1", "CYCB1-2", "CYCB2-1", "CYCB2-2")) %>%
  mutate(tissue = "Cyclins") %>%
  select(tissue, id = ID, name = alternative_name) %>%
  bind_rows(curated_genes)

# list of curated genes
curated_genes <- read_csv("data/raw/signature_genes_circle_plot.csv") %>%
  select(tissue = Expression, id = gene_ID, name = gene_name) %>%
  mutate(tissue = str_replace(tissue, " \\+ ", "\n&\n")) %>%
  mutate(tissue = str_replace_all(tissue, " ", "\n")) %>%
  mutate(id = toupper(str_trim(id)))


# Cell type --------------------------------------------------------------

# UMAP with clusters
ggplot(getReducedDim(ring, "UMAP30_MNN_logvst"),
       aes(V1, V2)) +
  geom_point(aes(colour = cluster_mnn_logvst),
             show.legend = FALSE) +
  geom_point(stat = "centroid", aes(group = cluster_mnn_logvst),
             shape = 21, size = 5, fill = "white", alpha = 0.7) +
  geom_text(stat = "centroid",
            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst)) +
  ggthemes::scale_colour_tableau(palette = "Tableau 20") +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Cluster") +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  coord_equal() + theme_void()

# plot for paper
p1 <- ring %>%
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("documents/pdf for figures/cell_type_annotation.pdf",
    width = 4.68*1.5, height = 3*1.5)
p1 + theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


