#!/bin/bash 
#SBATCH --qos=batch # job queue 
#SBATCH --ntasks=1 # number of nodes 
#SBATCH --cpus-per-task=1 # number of cores 
#SBATCH --mem=1G # memory pool for all cores 
#SBATCH --time=48:00:00 # wall clock time (D-HH:MM)

################################################################################
# Run GBRS on each sample. GRCm39. Ensembl 105.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-09-28
################################################################################

set -e -u -x -o pipefail


##### VARIABLES #####

# Base directory for project.
BASE_DIR=/projects/research_lab_directory/diabetes

# Sample metadata file.
SAMPLE_META_FILE=${BASE_DIR}/data/gbrs_sample_metadata.csv

# Base /flashscratch directory.
FLASHSCRATCH=/flashscratch/dgatti

# Temporary working directory.
TMP_DIR=${FLASHSCRATCH}/tmp

# Temporary output directory.
OUT_DIR=${FLASHSCRATCH}/results

# Final results directory (in backed up space)
DEST_DIR=${BASE_DIR}/results

# GBRS path.
GBRS_GITHUB_PATH=TheJacksonLaboratory/cs-nf-pipelines

# GBRS version on Github.
GBRS_VERSION=v0.4.2

# Singularity cache directory.
#export NXF_SINGULARITY_CACHEDIR=${flashscratch}/singularity_cache

##### MAIN #####


mkdir -p ${OUT_DIR} 
mkdir -p ${DEST_DIR} 
mkdir -p ${TMP_DIR} 
#mkdir -p ${NXF_SINGULARITY_CACHEDIR} 

module load singularity 

cd ${TMP_DIR} 

nextflow run ${GBRS_GITHUB_PATH} \
         -r ${GBRS_VERSION} \
         -profile sumner \
         -w ${TMP_DIR} \
         --workflow gbrs \
         --pubdir ${OUT_DIR} \
         --csv_input ${SAMPLE_META_FILE} \
         -resume

# Copy output files from OUT_DIR to DEST_DIR.
DIRS=`(ls ${OUT_DIR})`

for i in ${DIRS} 
do

  cp ${OUT_DIR}/${i}/stats/* ${DEST_DIR}
  cp ${OUT_DIR}/${i}/gbrs/*_counts ${DEST_DIR}
  cp ${OUT_DIR}/${i}/gbrs/*.tsv ${DEST_DIR}
  cp ${OUT_DIR}/${i}/gbrs/*.pdf ${DEST_DIR}

done

# Clean up tmp directory.
rm -rf ${TMP_DIR}
