
# load packages and read data ---------------------------------------------

library(SingleCellExperiment)
library(float)
library(dplyr)
library(tidyverse)

library(bigSCale)

library(WGCNA)
options(stringsAsFactors = FALSE);

sce <- readRDS("data/processed/SingleCellExperiment/ring_batches_strictfilt.rds")

trajectories_pseudotime <- read.csv("data/external/trajectories_pseudotime.csv")

# build gene network ------------------------------------------------------

expr <- assay(sce,"counts")
model <- compute.network.model(expr)

result <- compute.network(expr, gene.names = rownames(expr), model = model,
                          quantile.p =0.9,
                          polish = FALSE)

# add edge weights as correlation
correlations <- abs(dbl(result$correlations))
network <- result$graph

edgelist <- data.frame(get.edgelist(network),cor=0)
for (i in c(1:nrow(edgelist))) {
  edgelist[i,3] <- correlations[edgelist[i,1],edgelist[i,2]]
}

E(network)$weight <- edgelist[,3]


# module detection --------------------------------------------------------

modules_louvain <- cluster_louvain(network, weights = E(network)$weight)

# order by module size
gene_module <- data.frame(gene=V(network)$name,louvain=modules_louvain$membership)

module_summary <- data.frame(louvain=modules_louvain$membership) %>% 
  group_by(louvain) %>% 
  summarise(number_genes=n()) %>% 
  arrange(desc(number_genes)) %>% 
  mutate(module=paste0("M",1:nrow(.)),m=c(1:nrow(.)))

gene_module <- module_summary %>% 
  left_join(gene_module,.,by="louvain") %>% 
  select(gene,module,m)

V(network)$module <- gene_module$module
V(network)$louvain <- gene_module$m


# sub-modules in Module 1
network_m1 <- induced_subgraph(network,gene_module[gene_module$module=="M1",]$gene)

sub_module_m1 <- cluster_louvain(network_m1, weights = E(network_m1)$weight)$membership

gene_sub_module <- data.frame(gene=V(network_m1)$name,louvain=sub_module_m1)
# order by sub-module size
sub_module_summary <- data.frame(louvain=sub_module_m1) %>% 
  group_by(louvain) %>% 
  summarise(number_genes=n()) %>% 
  arrange(desc(number_genes)) %>% 
  mutate(sub_module=paste0("M1.",1:nrow(.)),subm=c(1:nrow(.)))

gene_module <- sub_module_summary %>% 
  left_join(gene_sub_module,.,by="louvain") %>% 
  select(gene,sub_module,subm) %>% 
  left_join(gene_module,.,by = "gene")

V(network)$sub_module <- gene_module$sub_module
V(network)$sub_louvain <- gene_module$subm


# eigengene ---------------------------------------------------------------

# module eigengenes
MEs_unscaled <- moduleEigengenes(t(as.matrix(assay(sce,"logcounts")[gene_module$gene,])), 
                                         colors = gene_module$m, 
                                         scale = FALSE,
                                         impute = FALSE,
                                         nPC = 10,
                                         align = "along average", 
                                         excludeGrey = TRUE, 
                                         grey = if (is.numeric(colors)) 0 else "grey",
                                         subHubs = FALSE,
                                         trapErrors = FALSE)

expr_eigengene <- MEs_unscaled$eigengenes %>% 
  as_tibble(rownames = "cell")
colnames(expr_eigengene) <- c("cell",paste0("eigen",1:(ncol(expr_eigengene)-1)))


# sub-eigengenes of Module 1
sub_MEs_unscaled <- moduleEigengenes(t(as.matrix(assay(sce,"logcounts")[V(network_m1)$name,])), 
                                 colors = V(network_m1)$sub_louvain, 
                                 scale = FALSE,
                                 impute = FALSE,
                                 nPC = 10,
                                 align = "along average", 
                                 excludeGrey = TRUE, 
                                 grey = if (is.numeric(colors)) 0 else "grey",
                                 subHubs = FALSE,
                                 trapErrors = FALSE)

expr_sub_eigengene <- sub_MEs_unscaled$eigengenes %>% 
  as_tibble(rownames = "cell")
colnames(expr_sub_eigengene) <- c("cell",paste0("eigen1.",1:(ncol(expr_sub_eigengene)-1)))

# prepare tables ----------------------------------------------------------

gene_module_assignment <- gene_module %>% 
  select(gene,module,sub_module)

module_summary <- module_summary %>% 
  select(module,number_genes) %>% 
  mutate(pct_variance=t(MEs_unscaled$varExplained[1,]))

source("analysis/functions/utils.R")
cell_eigengenes <- sce %>%
  getReducedDim(type = "UMAP30_MNN_logvst",
                colData = TRUE) %>%
  mutate(cell = paste0(Sample, "_", Barcode)) %>% 
  select(cell,V1,V2) %>% 
  left_join(trajectories_pseudotime,"cell") %>% 
  full_join(expr_eigengene, "cell") %>% 
  full_join(expr_sub_eigengene,"cell")


# save objects ------------------------------------------------------------

saveRDS(network,"data/processed/network/igraph_object.rds")

write.csv(gene_module_assignment,"data/processed/network/gene_module_assignment.csv",row.names = FALSE)

write.csv(module_summary,"data/processed/network/module_summary.csv",row.names = FALSE)

write.csv(cell_eigengenes,"data/processed/network/cell_eigengenes.csv",row.names = FALSE)
