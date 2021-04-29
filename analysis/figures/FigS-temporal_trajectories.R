# Setup -------------------------------------------------------------------

library(SingleCellExperiment)
library(tidyverse)
library(readxl)
library(patchwork)

# set ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 14)))

source("analysis/functions/utils.R")

# min-max scaling
min_max <- function(x, na.rm = FALSE){
  (x - min(x, na.rm = na.rm))/(max(x, na.rm = na.rm) - min(x, na.rm = na.rm))
}

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

# slingshot
sling <- readRDS("data/processed/trajectories/ring_batches_strictfilt_slingshot.rds")


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
rm(lon_cor, lon_matrix)


# Figure ------------------------------------------------------------------

# Plot
p1 <- ring %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  full_join(lon_classes, by = "cell") %>%
  mutate(section = str_remove(section, "SB_MEAN")) %>%
  mutate(section = as.numeric(str_remove(section, "L"))) %>%
  arrange(section) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = factor(section)), size = 0.5, alpha = 0.3) +
  guides(colour = guide_legend(override.aes = list(alpha=1, size = 2))) +
  scale_colour_viridis_d(option = "inferno", na.value = "lightgrey") +
  coord_equal() + theme_void() +
  # theme(legend.position = c(0, 0), legend.justification = c(0, -.2)) +
  labs(colour = "Brady et al.\nSection")

# empty plot just for filling space
p0 <- ring %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  mutate(cell = paste(Sample, Barcode, sep = "_")) %>%
  full_join(lon_classes, by = "cell") %>%
  mutate(section = str_remove(section, "SB_MEAN")) %>%
  mutate(section = as.numeric(str_remove(section, "L"))) %>%
  arrange(section) %>%
  ggplot(aes(V1, V2)) +
  coord_equal() + theme_void()

p2 <- sling %>%
  getReducedDim("UMAP30_MNN_logvst") %>%
  pivot_longer(matches("^slingPseudotime"),
               values_to = "pseudotime", names_to = "trajectory") %>%
  mutate(trajectory = gsub("slingPseudotime_", "Trajectory ", trajectory)) %>%
  arrange(trajectory, !is.na(pseudotime), pseudotime) %>%
  group_by(trajectory) %>%
  # mutate(pseudotime = min_max(pseudotime, na.rm = TRUE)) %>%
  ggplot(aes(V1, V2)) +
  geom_point(aes(colour = pseudotime), alpha = 0.1, size = 0.5) +
  # geom_point(stat = "centroid", aes(group = cluster_mnn_logvst), shape = 21, size = 5,
  #            fill = "white", alpha = 0.7) +
  # geom_text(stat = "centroid",
  #           aes(group = cluster_mnn_logvst, label = cluster_mnn_logvst)) +
  facet_wrap(~ trajectory, nrow = 2) +
  scale_colour_viridis_c(option = "inferno", na.value = "lightgrey") +
  coord_fixed() +
  labs(x = "UMAP1", y = "UMAP2", colour = "Pseudotime") +
  theme_void() +
  theme(legend.position = c(0.9, 0.15), legend.justification = c(1, 0))

((p1 | p0) / p2) + plot_layout(heights = c(1, 2))


# brady
pdf("documents/pdf for figures/temporal_annotation_brady.pdf",
    width = 3, height = 2.5)
p1
dev.off()

png("documents/pdf for figures/temporal_annotation_brady.png",
    width = 3, height = 2.5, units = "in", res = 300)
p1
dev.off()

# slingshot
pdf("documents/pdf for figures/temporal_annotation_slingshot.pdf",
    width = 6, height = 5)
p2
dev.off()

png("documents/pdf for figures/temporal_annotation_slingshot.png",
    width = 6, height = 5, units = "in", res = 300)
p2
dev.off()

