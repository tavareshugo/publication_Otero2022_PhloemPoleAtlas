# Code for analysis of scRNA-seq data of phloem ring

The file structure is as follows:

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
- `notebooks`: contains RMarkdown or other scripts that were used for exploratory analysis of the data. These are meant to be run interactively.


## (re-)Running the pipeline

To re-run the analysis (or if new samples are generated), this command alone should
re-run everything (make sure you `cd` to the project's directory):

`snakemake --use-conda`

Note for SLCU HPC users - to submit the jobs on the cluster, do the following instead:

- run `tmux` (this will launch a terminal that will persist even when you logout)
- run: `snakemake --use-conda --cluster "sbatch -c {threads} --mem-per-cpu=5G -J {rulename} -o logs/slurm/job#%j-{rulename}-{wildcards.sample}.log" --jobs 10`
- you can exit `tmux` by using the keyboard shortcut 'Ctrl + b' followed by 'd'
- when you log back in, to get back your `tmux` session back run `tmux attach -t 0`
