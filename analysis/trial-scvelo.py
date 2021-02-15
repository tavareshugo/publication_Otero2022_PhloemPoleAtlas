#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:29:20 2020

@author: hugo
"""

import scvelo as scv
import scanpy as sc

scv.set_figure_params()

#### Read data ####

# read data
#adata = scv.read("data/intermediate/velocyto/APL.loom", cache = True)
#adata = scv.read("data/intermediate/velocyto/ring_combined.loom")
#adata = scv.read("data/intermediate/velocyto/ring_combined_1000hvgs.loom")
#adata = scv.read("data/intermediate/velocyto/ring_combined_nofilt.loom")
#adata = scv.read("data/intermediate/velocyto/ring_ppp_clusters_nocycling.loom")  
#adata = scv.read("data/intermediate/velocyto/ring_combined_noouter.loom")
adata = scv.read("data/intermediate/velocyto/ring_strictfilt.loom")


#### Calculate velocities ####

# pre-processing - genes with at least 5 counts and expressed in at least 5 cells
scv.pp.filter_and_normalize(adata, min_shared_counts = 5, min_shared_cells=5, n_top_genes=2000)
scv.pp.moments(adata)

# fit dynamical model
scv.tl.recover_dynamics(adata)

# calculate velocities
scv.tl.velocity(adata)

# get latent time
scv.tl.latent_time(adata)

# calculate velocity graph
scv.tl.velocity_graph(adata)

# UMAP embedding
scv.tl.umap(adata)

# project arrows on embedding
scv.tl.velocity_embedding(adata, basis = 'umap')

# clustering
sc.tl.louvain(adata)

# velocity confidence
scv.tl.velocity_confidence(adata)

# PAGA trajectories
# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='louvain')


### Visualise ####

# umap
#sc.pl.umap(adata)

# check https://github.com/theislab/scvelo/issues/103
adata.obs['cluster_mnn_logvst'] # clusters are stored here
adata.obs['cluster_mnn_logvst'].cat.categories # levels of cluster
adata.uns['cluster_mnn_logvst_colors'] # colours currently assigned

# colours to match those used in other figures
tableau_pallete = ["#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6", "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295"]
tableau_pallete = [tableau_pallete[int(i)-1] for i in adata.obs['cluster_mnn_logvst'].cat.categories.tolist()]
adata.uns['cluster_mnn_logvst_colors'] = tableau_pallete

# figure for paper
scv.pl.velocity_embedding_stream(adata, basis='umap', color = 'cluster_mnn_logvst', title = "")


# project on embedding
scv.pl.velocity_embedding_stream(adata, basis='umap', color = 'louvain')
scv.pl.velocity_embedding_grid(adata, basis='umap', color = 'louvain', arrow_length = 3, arrow_size = 2)
scv.pl.velocity_embedding(adata, basis = 'umap', color = 'louvain', arrow_length = 3, arrow_size = 2)
scv.pl.velocity_embedding(adata, basis = 'umap', color = 'n_counts', arrow_length = 3, arrow_size = 2)

# with trajectory inferred by PAGA
scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5)

# coloured by latent time
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80)

# coloured by velocity
scv.pl.scatter(adata, c=('velocity_length', 'velocity_confidence'), cmap='coolwarm', perc=[5, 95])

# look at top genes
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False)

# look at marker genes
scv.tl.rank_velocity_genes(adata, groupby='louvain', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
scv.pl.scatter(adata, df['0'][:5])

# look at known genes
scv.pl.velocity(adata, ['AT2G22850',  'AT3G14570'], ncols=2) # S17 and CALS8 (PPP)
scv.pl.velocity(adata, ["AT3G12730", "AT5G57350", "AT1G22710"], ncols=2) # CC markers
scv.pl.velocity(adata, ["AT5G14750", "AT1G79580", "AT1G79840"], ncols=2) # outer layers
scv.pl.velocity(adata, ["AT2G37590", "AT5G02460", "AT1G54330", "AT1G05470"], ncols=2) # SE (early)
scv.pl.velocity(adata, ["AT5G17260", "AT1G06490"], ncols=2) # SE (late)
scv.pl.velocity(adata, ["AT4G39810", "AT5G17260"], ncols=2) # SE markers of enucleation
scv.pl.velocity(adata, ["AT2G26760"], ncols=2) # cyclin



#### Export ####

#adata.write("temp.h5ad")
