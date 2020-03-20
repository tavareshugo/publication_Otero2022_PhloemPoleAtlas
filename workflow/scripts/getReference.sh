#!/bin/bash

#### Input arguments (in order) ####

RELEASE="$1" # ENSEMBL genome release number
CPUS="$2"    # number of CPUs to use
OUTDIR="$3"  # output directory name


#### prepare directories ####

# create directory
mkdir -p "$OUTDIR"

# get GFP annotations
cp data/raw/GFP/* "$OUTDIR"

# work from this directory
cd "$OUTDIR"


#### Download reference genome ####

# Download genome from ENSEMBL
wget -O temp.fa.gz ftp://ftp.ensemblgenomes.org/pub/plants/release-${RELEASE}/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# Decompress
gunzip temp.fa.gz

# Add GFP proteins to reference
cat temp.fa GFP.fa > genome.fa
rm temp.fa


#### Download and prepare gene annotation ####

# download from ENSEMBL
wget -O temp.gff3.gz ftp://ftp.ensemblgenomes.org/pub/plants/release-${RELEASE}/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.${RELEASE}.gff3.gz

# Decompress
gunzip temp.gff3.gz

# Add GFP annotation
# remove annotations done by ENSEMBL (these are ncRNA, pre_miRNA, tRNA, etc)
# also the chromosomes are included in the annotation, this is also removed
cat temp.gff3 GFP.gff3 |
  grep -v "Ensembl_Plants" | \
  grep -v "The Arabidopsis Information Resource" > \
  annotation.gff3

rm temp.gff3

# Convert annotation to gtf format
# clean the conversion a bit
gffread annotation.gff3 -o - -T | \
  sed 's/transcript_id "transcript:/transcript_id "/' | \
  sed 's/gene_id "gene:/gene_id "/' | \
  sed 's/ gene_name .*//' \
  > annotation.gtf

# Extract transcript sequences from genome (can be used for SALMON index)
gffread -w annotation.fa -g genome.fa annotation.gtf

# create a transcript-to-gene map
cat annotation.gtf | sed 's/.*transcript_id "//' | sed 's/"; gene_id "/\t/' | sed 's/";//' > transcript_to_gene.tsv


#### cellranger indexing ####

cellranger mkref \
  --nthreads="${CPUS}" \
  --genome=cellranger_index \
  --fasta=genome.fa \
  --genes=annotation.gtf


#### SALMON indexing ####

# not sure about building this "gentrome". There's some conflicting instructions between 
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
# and
# https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode 
# The first has some script, which I'd done below - but this resulted in very low mapping rates 
# The latter has a more complex .sh script, which I didn't try running.

# grep "^>" genome.fa | cut -d " " -f 1 > decoys.txt
# sed -i -e 's/^>//g' decoys.txt
# cat annotation.fa genome.fa > gentrome.fa
# salmon index -t gentrome.fa -d decoys.txt -p $CPUS -i salmon_index --gencode

# indexing
# salmon index -t annotation.fa -p $CPUS -i salmon_index --gencode


