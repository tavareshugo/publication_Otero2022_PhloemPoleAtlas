# Setup -------------------------------------------------------------------

library(scater)
library(tidyverse)
library(readxl)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))

source("analysis/functions/utils.R")


# Download data ---------------------------------------------------------------

if(!file.exists("data/external/brady_et_al/ST12.xls")){
  # download
  download.file("https://science.sciencemag.org/highwire/filestream/588629/field_highwire_adjunct_files/1/1146265BradySTables.zip",
                destfile = "data/external/brady_et_al/sup_tables.zip")

  # unzip
  unzip("data/external/brady_et_al/sup_tables.zip",
        exdir = "data/external/brady_et_al/")

  # retain only files of interest
  file.rename("data/external/brady_et_al/Supplementary Tables_FINAL/ST12.xls",
              "data/external/brady_et_al/ST12.xls")

  # clear files
  unlink("data/external/brady_et_al/Supplementary Tables_FINAL", recursive = TRUE)
  unlink("data/external/brady_et_al/sup_tables.zip")
}


# Tidy Brady data ------------------------------------------------------

# read both sheets
brady <- list(
  longitudinal = read_xls("data/external/brady_et_al/ST12.xls",
                          sheet = "LONGITUDINAL", skip = 1) %>%
    pivot_longer(c(-Probe, -Gene),
                 names_to = "slice", values_to = "expr"),
  radial = read_xls("data/external/brady_et_al/ST12.xls",
                    sheet = "RADIAL", skip = 1) %>%
    pivot_longer(c(-Probe, -Gene),
                 names_to = "slice", values_to = "expr")
) %>%
  bind_rows(.id = "set") %>%
  rename_with(tolower) %>%
  # some probes don't have gene
  drop_na(gene)

# annotations for radial sections
radial_annot <- tribble(
  ~slice, ~tissue,
  "AGL42_MEAN", "QC",
  "APL_MEAN", "Phloem + CC",
  "COBL9_MEAN", "Hair cell",
  "CORTEX_MEAN", "Cortex",
  "gl2_MEAN", "Non-hair cell",
  "JO121_MEAN", "Xylem pole pericyle",
  "J0571_MEAN", "Ground tissue",
  "J2661_MEAN", "Mature pericycle",
  "LRC_MEAN", "LRC",
  "pet111_MEAN", "Columella",
  "rm1000_MEAN", "Lateral root",
  "S17_MEAN", "Phloem Pole Pericycle",
  "S18_MEAN", "Maturing xylem",
  "S4_MEAN", "Developing xylem",
  "S32_MEAN", "All protophloem",
  "scr5_MEAN", "Endodermis",
  "SUC2_MEAN", "Phloem CC",
  "wol_MEAN", "Stele",
  "xylem_2501_MEAN", "Stele"
)

# Read SCE data -----------------------------------------------------------

# sce object
ring <- readRDS("data/processed/SingleCellExperiment/ring_batches_strictfilt.rds")


# Longitudinal sections ----------------------------------------------------
# assign each cell to a section in Brady using correlation

# matrix of expression
lon_matrix <- brady %>%
  filter(set == "longitudinal") %>%
  pivot_wider(names_from = slice, values_from = expr) %>%
  # retain only genes with a single probe
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  # retain genes that also occur in our dataset
  filter(gene %in% rownames(ring)) %>%
  # convert to matrix
  select(-set, -probe) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# correlation matrix
lon_cor <- cor(as.matrix(assay(ring, "logvst")[rownames(lon_matrix), ]),
               lon_matrix,
               method = "spearman")

# get the Brady section with maximum correlation for each cell
lon_classes <- tibble(cell = rownames(lon_cor),
                      section = colnames(lon_cor)[max.col(lon_cor)],
                      cor = rowMax(lon_cor))

# Plot
ring %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  full_join(lon_classes, by = "cell") %>%
  mutate(section = str_remove(section, "SB_MEAN")) %>%
  mutate(section = as.numeric(str_remove(section, "L"))) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(ring, "UMAP30_MNN_logvst"), colour = "lightgrey") +
  geom_point(aes(colour = factor(section))) +
  scale_colour_viridis_d() +
  coord_equal() + theme_void() +
  labs(colour = "Brady et. al\nSection")

ring %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  full_join(lon_classes, by = "cell") %>%
  mutate(section = str_remove(section, "SB_MEAN")) %>%
  mutate(section = as.numeric(str_remove(section, "L"))) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(ring, "UMAP30_MNN_logvst"), colour = "lightgrey") +
  ggpointdensity::geom_pointdensity() +
  facet_grid(~ section) +
  scale_colour_viridis_c(option = "B") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")




# Radial section ----------------------------------------------------------

# matrix of expression
rad_matrix <- brady %>%
  filter(set == "radial") %>%
  pivot_wider(names_from = slice, values_from = expr) %>%
  # retain only genes with a single probe
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  # retain genes that also occur in our dataset
  filter(gene %in% rownames(ring)) %>%
  # convert to matrix
  select(-set, -probe) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# correlation matrix
rad_cor <- cor(as.matrix(assay(ring, "logvst")[rownames(rad_matrix), ]),
               rad_matrix,
               method = "spearman")

# get the Brady section with maximum correlation for each cell
rad_classes <- tibble(cell = rownames(rad_cor),
                      section = colnames(rad_cor)[max.col(rad_cor)],
                      cor = rowMax(rad_cor))

# Plot
getReducedDim(ring, "UMAP30_MNN_logvst") %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  full_join(rad_classes, by = "cell") %>%
  # mutate(section = ifelse(cor < 0.6, "unassigned", section)) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(ring, "UMAP30_MNN_logvst"), colour = "lightgrey") +
  geom_point(aes(colour = cor)) +
  scale_colour_viridis_c() +
  theme_void() +
  labs(colour = "Spearman rho") +
  facet_wrap(~ section)

# point density might be a better representation in this case
getReducedDim(ring, "UMAP30_MNN_logvst") %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  full_join(rad_classes, by = "cell") %>%
  full_join(radial_annot, by = c("section" = "slice")) %>%
  mutate(section = paste0(section, "\n", tissue)) %>%
  # mutate(section = ifelse(cor < 0.6, "unassigned", section)) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(ring, "UMAP30_MNN_logvst"), colour = "lightgrey") +
  ggpointdensity::geom_pointdensity() +
  scale_colour_viridis_c(option = "inferno") +
  theme_void() +
  labs(colour = "Point Density") +
  facet_wrap(~ section)


# Exclude "stele" sections
rad_classes <- tibble(cell = rownames(rad_cor),
                      section = colnames(rad_cor[, 1:17])[max.col(rad_cor[, 1:17])],
                      cor = rowMax(rad_cor[, 1:17]))

# point density might be a better representation in this case
getReducedDim(ring, "UMAP30_MNN_logvst") %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  full_join(rad_classes, by = "cell") %>%
  full_join(radial_annot, by = c("section" = "slice")) %>%
  mutate(section = paste0(section, "\n", tissue)) %>%
  # mutate(section = ifelse(cor < 0.6, "unassigned", section)) %>%
  ggplot(aes(V1, V2)) +
  geom_point(data = getReducedDim(ring, "UMAP30_MNN_logvst"), colour = "lightgrey") +
  ggpointdensity::geom_pointdensity() +
  scale_colour_viridis_c(option = "inferno") +
  theme_void() +
  labs(colour = "Point Density") +
  facet_wrap(~ section)

# Where is the xylem? Genes taken from Denyer Sup Table 1
getReducedDim(ring, "UMAP30_MNN_logvst",
              genes = c("AT1G20850", "AT4G23690", "AT4G35350", "AT5G03170"),
              melted = TRUE) %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  group_by(id) %>%
  mutate(expr_weighted = expr_weighted/max(expr_weighted, na.rm = TRUE)) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  # geom_label(stat = "centroid",
  #            aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
  #            alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression") +
  theme_void()



# Top 100 [deprecated] ----------------------------------------------------

# get the top 100 genes with highest variance
# maybe I should do by CV?
# should probably log-transform these data?
top100 <- brady %>%
  # remove probes that hit more than one gene
  filter(!str_detect(gene, ";")) %>%
  # calculate variance per gene
  group_by(set, gene) %>%
  summarise(var = var(expr)) %>%
  # get the gene ranks
  group_by(set) %>%
  mutate(rank = rank(-var)) %>%
  # get top 100
  filter(rank <= 100) %>%
  ungroup()

brady %>%
  inner_join(top100) %>%
  mutate(slice = fct_inorder(slice)) %>%
  filter(set == "longitudinal") %>%
  group_by(set, gene, slice) %>%
  filter(n() == 1) %>% ungroup() %>%
  ggplot(aes(slice, expr)) +
  geom_line(aes(group = gene)) +
  facet_wrap(~ set, scales = "free")

# not sure how to normalise these data...
brady %>%
  mutate(expr = log10(expr)) %>%
  group_by(set, gene) %>%
  summarise(expr_mean = mean(expr), expr_var = var(expr)) %>%
  ungroup() %>%
  ggplot(aes(expr_mean, expr_var)) +
  geom_point() +
  geom_abline() +
  facet_wrap(~ set)

# cluster by correlation
m <- brady %>%
  filter(set == "longitudinal") %>%
  inner_join(top100) %>%
  group_by(set, gene, slice) %>%
  filter(n() == 1) %>% ungroup() %>%
  select(gene, slice, expr) %>%
  pivot_wider(names_from = slice, values_from = expr) %>%
  column_to_rownames("gene") %>%
  as.matrix()

plot(hclust(as.dist(1 - cor(t(m), method = "spearman"))))

# visualise expression of a few genes
ring %>%
  getReducedDim("UMAP30_MNN_logvst",
                genes = top100$gene[1:9], melted = TRUE,
                exprs_values = "logvst") %>%
  arrange(!is.na(expr)) %>%
  group_by(cluster_mnn_logvst, id) %>%
  mutate(expr_weighted = expr * (sum(!is.na(expr))/n())) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = expr_weighted)) +
  geom_label(stat = "centroid",
             aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst),
             alpha = 0.8, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ id) +
  scale_colour_viridis_c(na.value = "lightgrey") +
  labs(x = "UMAP1", y = "UMAP2",
       colour = "Cluster-weighted\nNormalised\nExpression") +
  theme_void()

brady %>%
  inner_join(top100) %>%
  mutate(slice = fct_inorder(slice)) %>%
  filter(set == "longitudinal") %>%
  group_by(set, gene, slice) %>%
  filter(n() == 1) %>% ungroup() %>%
  filter(gene %in% top100$gene[1:9]) %>%
  ggplot(aes(slice, expr)) +
  geom_line(aes(group = gene)) +
  geom_point() +
  facet_wrap(gene ~ set, scales = "free")

brady %>%
  inner_join(top100) %>%
  mutate(slice = fct_inorder(slice)) %>%
  filter(set == "radial") %>%
  group_by(set, gene, slice) %>%
  filter(n() == 1) %>% ungroup() %>%
  ggplot(aes(slice, expr)) +
  geom_line(aes(group = gene)) +
  facet_wrap( ~ set, scales = "free") +
  coord_flip()
