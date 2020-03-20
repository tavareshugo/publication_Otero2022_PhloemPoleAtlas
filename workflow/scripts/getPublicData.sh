#!/bin/bash

# Download data from Denyer et al 2019 - https://doi.org/10.1016/j.devcel.2019.02.022
# GEO GSE123818 - https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP173393
# there are several samples in the repository, we download the ones labelled "WT rep1" and "WT rep2"
# data are provided as 10x barcoded BAMs - https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam

# prepare output directory
mkdir -p data/external/denyer2019/
cd data/external/denyer2019/


#### Download data ####
echo "Starting data download..."

# WT rep 1: run SRR10136142
wget -O denyer1.bam https://sra-pub-src-2.s3.amazonaws.com/SRR10136142/Rep1_col_root_wt_postsorted_genome_bam.bam.1

# WT rep 2: run SRR10136143
wget -O denyer2.bam https://sra-pub-src-2.s3.amazonaws.com/SRR10136143/Rep2_MIR166_root_wt_postsorted_genome_bam.bam.1


#### Convert to FASTQ ####
echo "Starting bam-to-fastq conversion..."

# see https://support.10xgenomics.com/docs/bamtofastq
wget -O bamtofastq http://cf.10xgenomics.com/misc/bamtofastq-1.2.0
chmod 700 bamtofastq

# convert
bamtofastq denyer1.bam rep1
bamtofastq denyer2.bam rep2

