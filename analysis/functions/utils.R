#' fortify S3 method for BioConductor::DataFrame
#'
#' This allows using DFrame with ggplot2
#'
#' @param data a DataFrame object to plot.
fortify.DFrame <- function(data) as.data.frame(data)

#' helper function to fetch reducedDim data
#' @param x a SingleCellExperiment object.
#' @param type the name of the dimensionality reduction result in x.
#' @param colData logical indicating whether to include colData or not.
#' @param genes a vector of gene names to get logcounts for.
#' @param exprs_values the assay slot to retrive gene expression values from.
#' @param as_dt a logical indicating whether the object retured should be a `data.table`.
#' @param melted whether to return the table in a long format.
#' @param zeros_to_na whether zero expression should be converted to NA.
getReducedDim <- function(x, type = "UMAP", colData = TRUE, genes = NULL,
                          exprs_values = "logcounts",
                          as_dt = TRUE, melted = FALSE, zeros_to_na = TRUE){
  stopifnot(type %in% reducedDimNames(x))

  res <- data.table::as.data.table(reducedDim(x, type = type))

  # fetch colData
  if(colData){
    res <- cbind(res, colData(x))
  }

  # fetch genes
  if(any(genes %in% rownames(x))){
    if(sum(genes %in% rownames(x)) < length(genes)){
      warning("Not all genes are present in the data")
      genes <- genes[which(genes %in% rownames(x))]
    }
    genes_idx <- which(rownames(x) %in% genes)

    if(length(genes_idx) == 1){
      gene_cts <- as.matrix(assay(x, exprs_values)[genes_idx, ])
    } else if(length(genes_idx) > 1){
      gene_cts <- as.matrix(t(assay(x, exprs_values)[genes_idx, ]))
    }
    colnames(gene_cts) <- rownames(x)[genes_idx]

    res <- cbind(res, gene_cts)
  }

  if(melted & any(genes %in% colnames(res))){
    res <- data.table::as.data.table(res)
    res <- data.table::melt(res, measure.vars = genes,
                            variable.name = "id", value.name = "expr")

    if(zeros_to_na){
      res$expr[which(res$expr == 0)] <- NA
    }
  }

  if (as_dt) {
    return(data.table::as.data.table(res))
  } else {
    return(as.data.frame(res))
  }
}


#' Plotting centroids
#'
#' this is a ggplot2 extension that calculates centroids of x and y aesthetics
StatCentroid <- ggplot2::ggproto("StatCentroid", Stat,
                                 compute_group = function(data, scales) {
                                   data.frame(x = mean(data$x), y = mean(data$y))
                                 },

                                 required_aes = c("x", "y")
)

stat_centroid <- function(mapping = NULL, data = NULL, geom = "point",
                          position = "identity", na.rm = FALSE, show.legend = NA,
                          inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatCentroid, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
