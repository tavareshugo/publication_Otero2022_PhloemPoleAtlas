# specific plotting functions for our UMAP plots
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
    coord_equal() + theme_void() +
    theme_void()
}
