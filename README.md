# Code for analysis of scRNA-seq data of phloem ring

Assuming full pipeline is run, the file structure is as follows:

- `data`: contains all the data needed for the project
  - `raw`: the raw data files. This should be shared in public repositories upon publication.
  - `intermediate`: intermediate files used in the analysis. These can be re-created by running the workflow and can be removed to save space.
  - `processed`: data files that can also be re-created by running the workflow, but we keep them as they are used for making key figures/analysis. These should also be shared in a public repository upon publication.
  - `external`: contains external data downloaded from web (e.g. reference genome).
- `logs`: contains log files resulting from running different steps of the workflow.
- `workflow`: contains the actual workflow scripts. We're using `snakemake` for our pipelines. The file `pipeline.pdf` contains a diagram of the workflow.
  - `Snakefile`: is the script with the workflow specification
  - `scripts`: contains scripts used in the pipeline
  - `envs`: contain _YAML_ files with the conda environment for each step of the pipeline.

## Data

- 10x data generated in this project is deposited on SRA: [accession GSE181999](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181999).
- The list of public datasets used for the integrated analysis is available [here](public_datasets.csv).
- For a quick look at our data, two files with processed results are attached to [this repository's release](https://github.com/tavareshugo/publication_Otero2022_PhloemPoleAtlas/releases/latest). 
  These are detailed below. 
  
### cell_metadata.csv

A CSV file containing information about our cell annotation.
This can be opened in any spreadsheet or data analysis software.
The columns in the table are:

- `cell_id` is a unique identifier of each cell (a combination of the sorted sample identifier and the 10x barcode).
- `sample` is the sorted sample identifier (refers to the marker used for cell-sorting).
- `barcode` is the 10x barcode of the cell.
- `total_counts` is the total UMI count detected in the cell.
- `detected_genes` is the number of genes with detected expression (at least 1 UMI).
- `cluster` is the cluster number used throughout the paper.
- `annotation` is the cell type annotation inferred in our paper.
- `UMAP1` is the first axis of the UMAP projection used in the paper.
- `UMAP2` is the second axis of the UMAP projection used in the paper.


### SingleCellExperiment_filtered.rds

This is the `SingleCellExperiment` object (for use with R/Bioconductor packages), containing the full analysis of the data.
The object contains:

- Matrices of raw and normalised counts (in the `assays` slots).
- Dimensionality reduction matrices including UMAP, t-SNE and diffusion maps (in the `reducedDims` slots).
- The trajectory analysis results from slingshot, whose pseudotime can be found in the `colData` slot of the object.

Here is some R code to read this object: 

```R
library(SingleCellExperiment)

# read the object
sce <- readRDS("SingleCellExperiment_filtered.rds")

# cluster assignment used in paper
colData(sce)$cluster_mnn_logvst

# assays
assay(sce, "counts")    # raw counts
assay(sce, "logvst")    # counts normalised using sctransform
assay(sce, "logcounts") # counts normalised using scran's deconvolution method
```


## (re-)Running the pipeline

:warning: This pipeline is not being maintained and some of the code may not run as expected.

To re-run the analysis (or if new samples are generated), this command alone should
re-run everything (make sure you `cd` to the project's directory):

`snakemake --use-conda`

Note for SLCU HPC users - to submit the jobs on the cluster, do the following instead:

- run `tmux` (this will launch a terminal that will persist even when you logout)
- run: `snakemake --use-conda --cluster "sbatch -c {threads} --mem-per-cpu=5G -J {rulename} -o logs/slurm/job#%j-{rulename}-{wildcards.sample}.log" --jobs 10`
- you can exit `tmux` by using the keyboard shortcut 'Ctrl + b' followed by 'd'
- when you log back in, to get back your `tmux` session back run `tmux attach -t 0`
