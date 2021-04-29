# Setup -------------------------------------------------------------------

library(scater)
library(tidyverse)
library(readxl)
library(UpSetR)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))

# source custom functions
source("analysis/functions/utils.R")


# Read data ---------------------------------------------------------------

# get wilcoxon test results for each cluster
cluster_test <- read_csv("data/processed/gene_sets/ring_hardfilt_cluster_markers.csv",
                         guess_max = Inf) %>%
  mutate(cluster = factor(cluster))

# make transcription factor list for convenience
tfs <- cluster_test %>%
  filter(is_tf) %>%
  pull(id) %>% unique()

# sce
all_batches <- readRDS("data/processed/SingleCellExperiment/all_batches_hardfilt.rds")


# Read PEAR targets -------------------------------------------------------

# Read gene lists from excel file
pear_targets <- list(
  `PEAR2 1h` = read_xlsx("data/raw/Upregulated PEAR targets.xlsx",
                       sheet = "PEAR1-GR and PEAR2-GR",
                       range = "A3:A14"),
  `PEAR1 2h` = read_xlsx("data/raw/Upregulated PEAR targets.xlsx",
                       sheet = "PEAR1-GR and PEAR2-GR",
                       range = "A20:A232"),
  `PEAR2 2h` = read_xlsx("data/raw/Upregulated PEAR targets.xlsx",
                       sheet = "PEAR1-GR and PEAR2-GR",
                       range = "A236:A671"),
  `PEAR2 direct` = read_xlsx("data/raw/Upregulated PEAR targets.xlsx",
                         sheet = "putative direct targets",
                         range = "A3:A87"),
  `PEAR1 direct` = read_xlsx("data/raw/Upregulated PEAR targets.xlsx",
                              sheet = "putative direct targets",
                              range = "A92:A112"),
  `PEAR1-TMO5_double_ox 1h` = read_xlsx("data/raw/Upregulated PEAR targets.xlsx",
                              sheet = "PEAR1-TMO5 double ox",
                              range = "A3:A21"),
  `PEAR1-TMO5_double_ox 2h` = read_xlsx("data/raw/Upregulated PEAR targets.xlsx",
                                         sheet = "PEAR1-TMO5 double ox",
                                         range = "A27:A159"),
  `PEAR1-TMO5_double_ox 3h` = read_xlsx("data/raw/Upregulated PEAR targets.xlsx",
                                        sheet = "PEAR1-TMO5 double ox",
                                        range = "A164:A474"),
  `PEAR1-TMO5_double_ox 4h` = read_xlsx("data/raw/Upregulated PEAR targets.xlsx",
                                        sheet = "PEAR1-TMO5 double ox",
                                        range = "A479:A1006")
) %>%
  bind_rows(.id = "set") %>%
  mutate(AGI = toupper(AGI))

# PAPL targets
papl_targets <- read_xlsx("data/raw/ccc_vs_wt_upregulated.xlsx") %>%
  mutate(Gene = toupper(Gene))

# all together
pear_papl <- papl_targets %>%
  mutate(set = "ccc vs WT") %>%
  # only genes with increased expression in mutant
  filter(log2foldchange > log2(1.5)) %>%
  select(set, AGI = Gene) %>%
  bind_rows(pear_targets)


# Explore intersections ---------------------------------------------------

# intersections of all datasets
pear_papl %>%
  with(split(AGI, set)) %>%
  fromList() %>%
  upset(nsets = 10, order.by = "freq", text.scale = 2)

# focus on PAPL and PEAR 2h only
pear_papl %>%
  filter(set %in% c("ccc vs WT", "PEAR2 2h")) %>%
  with(split(AGI, set)) %>%
  fromList() %>%
  upset(nsets = 10, order.by = "freq", text.scale = 2)

targets <- pear_papl %>%
  filter(set %in% c("ccc vs WT", "PEAR2 2h")) %>%
  group_by(AGI) %>%
  filter(n() == 2) %>%
  distinct(AGI) %>%
  pull(AGI)

# Add PAPLs (AT1G29160, AT2G34140), PEAR1 (AT2G37590) and PEAR2 (AT5G02460)
targets <- unique(c("AT1G29160", "AT2G34140", "AT2G37590", "AT5G02460", targets))

# plot all genes for browsing (there's a lot of them, so exporting to a png)
png("temp.png", width = 480*2, height = 2*200*length(targets)/4, res = 200)
all_batches %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = targets,
                melted = TRUE) %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  # weight the expression by number of cells expressing it
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  group_by(id) %>%
  # min-max normalisation for scaling purposes
  mutate(expr_weighted = (expr_weighted - min(expr_weighted, na.rm = TRUE))/max(expr_weighted, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  theme_void() +
  labs(colour = "scaled\nexpression") +
  facet_wrap(~ id, ncol = 4)
dev.off()




# # same for expression violins
# png("temp.png", width = 480*2, height = 2*200*length(targets)/2, res = 200)
# plotExpression(object = all_batches,
#                features = targets,
#                x = "cluster_mnn_logvst", "logcounts",
#                ncol = 2)
# dev.off()


# Compare with single-cell ------------------------------------------------

# cross-reference
cluster_test %>%
  filter(FDR < 0.05 & summary.AUC > 0.7) %>%
  full_join(pear_targets, by = c("id" = "AGI")) %>%
  select(cluster, id, set) %>%
  mutate(set = ifelse(is.na(set), "None", set)) %>%
  with(split(id, set)) %>%
  fromList() %>%
  upset(nsets = 90, order.by = "freq")

cluster_test %>%
  filter(FDR < 0.05 & summary.AUC > 0.7) %>%
  full_join(pear_targets, by = c("id" = "AGI")) %>%
  select(cluster, id, set) %>%
  filter(set == "PEAR2 2h") %>%
  mutate(cluster = ifelse(is.na(cluster), "None", cluster)) %>%
  with(split(id, cluster)) %>%
  fromList() %>%
  upset(nsets = 90, order.by = "freq")


# what proportion of these is in our data at all?
# i.e. they are expressed in at least a few cells
pear_targets %>%
  group_by(set) %>%
  summarise(nexpr = sum(AGI %in% cluster_test$id),
            n = n()) %>%
  mutate(prop = nexpr/n)

# need to check if ~90% is more than expected by chance
# there's a total of 32548 nuclear genes
nrow(ring)/32548

# a weird way to do this
expressed_ids <- unique(cluster_test$id)
nonexpressed_ids <- rownames(unfilt)[(!(rownames(unfilt) %in% cluster_test$id))]
pear2targets <- pear_targets %>% filter(set == "PEAR2 2h") %>%
  pull(AGI) %>% unique()

matrix(c(
  sum(pear2targets %in% expressed_ids),
  sum(pear2targets %in% nonexpressed_ids),
  sum(!(expressed_ids %in% pear2targets)),
  sum(!(nonexpressed_ids %in% pear2targets))
), byrow = TRUE, nrow = 2)


