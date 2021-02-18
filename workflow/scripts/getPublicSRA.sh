#!/bin/bash
#SBATCH -A LEYSER-SL2-CPU
#SBATCH -D /rds/user/hm533/rds-ol235-leyser-hpc/projects/2020-phloem_scRNAseq
#SBATCH -J public_data
#SBATCH -o logs/getPublicData/getPublicSRA_%a.log
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -p skylake-himem
#SBATCH --mem-per-cpu=12030MB
#SBATCH -a 12,13,16,17,21,65 # number of rows in the CSV file with SRA information

# activate conda environment
source activate sra

# get sample information from the table
INFO=$(head -n "$SLURM_ARRAY_TASK_ID" data/external/public_datasets.csv | tail -n 1)

# split relevant parts
PUBLICATION=$(echo $INFO | cut -d "," -f 1)
SAMPLE=$(echo $INFO | cut -d "," -f 2)
RUN=$(echo $INFO | cut -d "," -f 3)
TYPE=$(echo $INFO | cut -d "," -f 5)
LINK=$(echo $INFO | cut -d "," -f 6)

# move to data folder
cd data/external/


#### download and process ####

# Denyer data are provided as bam files
if [[ $PUBLICATION == "denyer" ]]
then
  OUTDIR="${PUBLICATION}/"

  # create directory structure
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"

  # download file
  NAME="${PUBLICATION}-${SAMPLE}"
  wget -O "${NAME}.bam" "$LINK"

  # see https://support.10xgenomics.com/docs/bamtofastq
  wget -O bamtofastq http://cf.10xgenomics.com/misc/bamtofastq-1.2.0
  chmod 700 bamtofastq

  # convert to fastq
  bamtofastq "${NAME}.bam" "${SAMPLE}"

  # rename the sub-directories created to run1 and run2
  DIRS=$(ls ${SAMPLE})
  i=1
  for DIR in $DIRS
  do
    mv "${SAMPLE}/$DIR" "${SAMPLE}/run${i}"
    i=$((i + 1))
  done

  # remove original file
  #rm "${NAME}.bam"
fi


# Wendrich data are fastq files
if [[ $PUBLICATION == "wendrich" ]]
then
  OUTDIR="${PUBLICATION}/${SAMPLE}/"

  # create directory structure
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"

  # download file
  NAME="${PUBLICATION}-${SAMPLE}-${RUN}"
  #wget -O "${NAME}.sra" "$LINK"

  # convert to fastq
  fastq-dump -O "${RUN}" --split-files "${NAME}.sra"

  # rename files for cellranger
  # 1 = I1
  # 2 = R1
  # 3 = R2
  mv "${RUN}/${NAME}_1.fastq" "${RUN}/fastqdump_S1_L001_I1_001.fastq"
  mv "${RUN}/${NAME}_2.fastq" "${RUN}/fastqdump_S1_L001_R1_001.fastq"
  mv "${RUN}/${NAME}_3.fastq" "${RUN}/fastqdump_S1_L001_R2_001.fastq"

  # compress the files
  gzip -r "${RUN}"

  # remove original file
  #rm "${NAME}.sra"
fi


# Shahan data are fastq files
if [[ $PUBLICATION == "shahan" ]]
then
  OUTDIR="${PUBLICATION}/${SAMPLE}/"

  # create directory structure
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"

  # download file
  NAME="${PUBLICATION}-${SAMPLE}-${RUN}"
  #wget -O "${NAME}.sra" "$LINK"

  # convert to fastq
  fastq-dump -O "${RUN}" --split-files "${NAME}.sra"

  # rename files for cellranger
  # 1 = R1
  # 2 = R2
  # 3 = I1
  mv "${RUN}/${NAME}_1.fastq" "${RUN}/fastqdump_S1_L001_R1_001.fastq"
  mv "${RUN}/${NAME}_2.fastq" "${RUN}/fastqdump_S1_L001_R2_001.fastq"
  mv "${RUN}/${NAME}_3.fastq" "${RUN}/fastqdump_S1_L001_I1_001.fastq"

  # compress the files
  gzip -r "${RUN}"

  # remove original file
  #rm "${NAME}.sra"
fi


# Kim data are fastq files
if [[ $PUBLICATION == "kim" ]]
then
  OUTDIR="${PUBLICATION}/${SAMPLE}/"

  # create directory structure
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"

  # download file
  NAME="${PUBLICATION}-${SAMPLE}-${RUN}"
  #wget -O "${NAME}.sra" "$LINK"

  # convert to fastq
  fastq-dump -O "${RUN}" --split-files "${NAME}.sra"

  # rename files for cellranger
  # 1 = R1
  # 2 = R2
  mv "${RUN}/${NAME}_1.fastq" "${RUN}/fastqdump_S1_L001_R1_001.fastq"
  mv "${RUN}/${NAME}_2.fastq" "${RUN}/fastqdump_S1_L001_R2_001.fastq"

  # compress the files
  gzip -r "${RUN}"

  # remove original file
  #rm "${NAME}.sra"
fi


