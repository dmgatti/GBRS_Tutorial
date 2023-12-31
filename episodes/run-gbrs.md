---
title: "Running GBRS"
teaching: 60
exercises: 60
---

:::::::::::::::::::::::::::::::::::::: questions 

- How do I run GBRS on my FASTQ files at JAX?
- How do I run GBRS on my FASTQ files at another institution?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- At JAX: analyze FASTQ files using GBRS on the High Performance Computing cluster.
- At other institutions: configure Nexflow and the reference files and analyze
FASTQ files using GBRS.

::::::::::::::::::::::::::::::::::::::::::::::::

## At the Jackson Laboratory (JAX)

The GBRS Nextflow pipeline is configured to run on *sumner*, the High Performance 
Computing (HPC) cluster at JAX. 

We will consider two scenarios:

1. You are analyzing expression data from Diversity Outbred mice;
1. You are analyzing expression data from some other cross;

### GBRS for Diversity Outbred Mice

The Next-Generaiong Operations (NGSOps) team has configured the default 
arrguments for the GBRS Nextflow pipeline to use the reference files for GRCm39
and Ensembl version 105. This means that the arguments that point GBRS to the 
locations of the founder genomes, founder transcriptomes, and aligner indices
are already set.

Here is the entire script:

```
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
```


For Diversity Outbred (DO) mice, the reference files are stored in the 
reference data directories `/projects/omics_share` on *sumner*.

JAX uses [slurm](https://jacksonlaboratory.sharepoint.com/sites/ResearchIT/SitePages/Submitting-Basic-Jobs-with-SLURM.aspx)
(NOTE: This is an internal JAX link which requires a JAX login.) to manage 
computing jobs on *sumner*. There are several good tutorials on using *slurm*
on the JAX [Research IT Documentation Library](https://jacksonlaboratory.sharepoint.com/sites/ResearchIT/SitePages/Documentation.aspx).

#### *slurm* Options Block

*slurm* scripts are *bash* scripts, so they start with:

```
#!/bin/bash
```

Next, we use a *slurm* options block to tell *slurm* what computing resources
we would like to use. You could pass *slurm* options in at the command line, but
we recommend using an options block in your script so that you have a record of
the resources required to run your job. 

Each line in the *slurm* options block starts with:

```
#SBATCH
```

This tells *slurm* that what follows on this line is a *slurm* option. There 
are many *slurm* options, which are exhaustively enumerated in the 
[slurm documentation](https://slurm.schedmd.com/sbatch.html). Here we will use
only a few of the options:

```
#SBATCH --qos=batch       # job queue
#SBATCH --ntasks=1        # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
#SBATCH --mem=1G          # memory pool for all cores
#SBATCH --time=48:00:00   # wall clock time (D-HH:MM)
```

The first option specifies the job queue to use for this job. In this case, it
is the "batch" queue on *sumner*.

```
#SBATCH --qos=batch       # job queue
```

The next two options specify the number of nodes and number of CPUs (or cores) 
per node to use. In this case, we will use one node and one CPU. This may seem 
confusing because we would like to use many nodes and cores to reduce the total
run time of our job. But Nextflow handles resource request. This script just
starts the Nextflow job, which does not require many resources. If you request
more nodes and CPUs here, you will only be tying up resources that will sit 
idle while your job run. So **one** node and **one** CPU per node are sufficient.

```
#SBATCH --ntasks=1        # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
```

The next option specifies the memory required for this task. We will request 
one Gigabyte (G). Again, this may seem like a small amount of memory for a 
large computational job, but this is only the memeory required to start and
manage the Nextflow job. Nextflow will request the resources that it needs to 
run GBRS. 

```
#SBATCH --mem=1G          # memory pool for all cores
```

The last option that we will speicify is the total wall-clock time needed to 
complete the job. This needs to be the total amount of time required to complete
the job. This can be difficult to estimate and depends upon the number of 
samples, depth of coverage, and the configuration and utilization of the cluster.

> DMG: TBD: Should we run 50, 100, 200 & 300 samples and include a time/samples plot? 
I have 350 samples that we could use for this. It's not going to be linear with
more sample, but it's something.

```
#SBATCH --time=48:00:00   # wall clock time (D-HH:MM)
```

#### Header Section

```
################################################################################
# Run GBRS on each sample. GRCm39. Ensembl 105.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-09-28
################################################################################
```

This section is not required, but it allows you to add some comments about what
are doing. It is also a consistent place to put your name, contact information,
and the date. Putting your name and contact information helps to keep you 
accountable for the code. Putting the date in the header helps you to remember
when you wrote the script.

#### Error Handling

```
set -e -u -x -o pipefail
```

This commands tells *bash* to exit if any line of the script fails (-e), to 
catch uninitialized variables (-u), to print commands as we run (-x), and to  
output the last exit code (-o).

#### Variable Section

We strongly recommend that you use *bash* variables to build file paths and
specify arguments to the GBRS script. The reason for this is that it is easier
to remember what each varible is when it is given a meaningful name.

Consider this script which runs the "analysis_tool_of_science":

```
analysis_tool_of_science -t 8 \
  -g /path/to/another/file \
  -a /really/long/path/to/yet/another/file \
  -i /path/to/mystery/file/number1 /path/to/mystery/file/number2 \
  -x /path/to/some/directory/on/the/cluster
  -o /path/to/another/directory/with/a/file/prefix
```

You can read the tool's documentation to find out what each file is, but that
requires extra time. This script would be easier to read if the paths had 
meaningful names. In the code block below, we give each file a meaningful 
variable name and provide a comment explaining what each file is. We also use
this opportunity to specify the genome and annotation builds.

```
# FASTA file containing genome sequence for GRCm39.
GENOME_FASTA=/path/to/another/file

# GTF file containing genome annotation for Ensembl 105.
ANNOTATION_GTF=/really/long/path/to/yet/another/file

# Paths to the forward and reverse FASTQ files.
FASTQ_FILE1=/path/to/mystery/file/number1 
FASTQ_FILE2=/path/to/mystery/file/number2

# Path to temporary directory for analysis_tool.
TEMP_DIR=/path/to/some/directory/on/the/cluster

# Path to output file, with a file prefix.
OUTPUT_PATH=/path/to/another/directory/with/a/file/prefix

analysis_tool_of_science -t 8 \
  -g ${GENOME_FASTA} \
  -a ${ANNOTATION_GTF} \
  -i $(FASTQ_FILE1} ${FASTQ_FILE2} \
  -x ${TEMP_DIR}
  -o ${OUTPUT_PATH}
```

This makes the script much easier for yourself (a year from now) and others
to read and understand.

First, we will focus on variables that set the location of your files and 
working directories.

```
# Base directory for project.
BASE_DIR=/projects/research_lab_directory/diabetes
```

We like to set a "base directory", which is the high-level directory for the
project on which you are working. In this case, we are using a diabetes project
as an example. We use the base directory to create other file paths which
are below this directory.

```
# Sample metadata file.
SAMPLE_META_FILE=${BASE_DIR}/data/sodo_gbrs_metadata.csv
```

This is the path to the sample metadata file, which contains information
about your samples. This includies the sample ID, sex, outbreeding generation,
sequencer lane, and paths to the forward and reverse FASTQ files. Each row will
contain information for one sample.

The `sampleID` column should contain values that you use to identify the samples.

The `sex` column should contain the mouse's sex recorded as 'F' for females or 
'M' for males.

The `generation` column should contain the outbreeding generation for the mouse.

The `lane` column should contain the lane number on which the sample was run.
This can sometimes be found in the FASTQ file name. In the example below, it is
recorded as 'L004' in the FASTQ file names and this is what we used in the 
`lane` column.

The `fastq_1` and `fastq_2` columns contain the full path to the forward and 
reverse FASTQ files. In this case, we copied the files to /flashscratch, which
is the temporary sratch storage space for large files.

Below is an example of the top of a sample metadata file. It is comma-delimited
and has column names as the top.

```
sampleID,sex,generation,lane,fastq_1,fastq_2
D2,F,46,L004,/flashscratch/dgatti/fastq/D2_GT22-11729_CTAGGCAT-ACAGAGGT_S224_L004_R1_001.fastq.gz,/flashscratch/dgatti/fastq/D2_GT22-11729_CTAGGCAT-ACAGAGGT_S224_L004_R2_001.fastq.gz
D3,F,46,L001,/flashscratch/dgatti/fastq/D3_GT22-11665_CGGCTAAT-CTCGTTCT_S21_L001_R1_001.fastq.gz,/flashscratch/dgatti/fastq/D3_GT22-11665_CGGCTAAT-CTCGTTCT_S21_L001_R2_001.fastq.gz
D4,F,46,L001,/flashscratch/dgatti/fastq/D4_GT22-11670_TACGCTAC-CGTGTGAT_S65_L001_R1_001.fastq.gz,/flashscratch/dgatti/fastq/D4_GT22-11670_TACGCTAC-CGTGTGAT_S65_L001_R2_001.fastq.gz
D5,F,46,L001,/flashscratch/dgatti/fastq/D5_GT22-11708_GAGCAGTA-TGAGCTGT_S57_L001_R1_001.fastq.gz,/flashscratch/dgatti/fastq/D5_GT22-11708_GAGCAGTA-TGAGCTGT_S57_L001_R2_001.fastq.gz
```

Notice that we have copied the FASTQ files over to /flashscratch. This isn't 
required but we do this to simplify the file paths.

```
# Base /flashscratch directory.
FLASHSCRATCH=/flashscratch/username
```

This is the path to your directory on /flashscratch. You should create this 
before you begin your analysis. And you should replace 'username' with your
username.

```
# Temporary directory.
TMP_DIR=${FLASHSCRATCH}/tmp
```

This is the path to the temporary directory where Nextflow and GBRS can store
temporary working files during the analysis.

```
# Output directory.
OUT_DIR=${FLASHSCRATCH}/results
```

This is the directory where you would like GBRS to write your results. You 
should place this on /flashscratch and then copy the files that you need over
to /projects when the analysis is complete.


```
# Final results directory (in backed up space)
DEST_DIR=${BASE_DIR}/results/gbrs
```

This is the directory to which you will copy your final results. It should be
a location that is backed up (like /projects). Remember, /flashscratch is not
backed up and your files will be deleted without notice after 10 days.

Once you have set all of the paths and created your sample metadata file, you 
can submit the script to the *slurm* queue. Save the script to a file called
'run_gbrs.sh', and then submit it using 'sbatch'.

```
sbatch run_gbrs.sh
```

### GBRS for Other Mouse Crosses

workflow GBRS {
    // Step 0: Download data and concat Fastq files if needed. 
    meta_ch = null

    if (params.download_data){
        FILE_DOWNLOAD(ch_input_sample)

        if (params.read_type == 'PE'){
            FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[2][0], 'R1']}.set{r1}
            FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[2][1], 'R2']}.set{r2}
            read_ch = r1.mix(r2)
        } else if (params.read_type == 'SE'){
            FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[2][0], 'R1']}.set{read_ch}
        }

        FILE_DOWNLOAD.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
    }

    // Step 00: Concat local Fastq files from CSV input if required.
    if (!params.download_data && params.csv_input){
        CONCATENATE_LOCAL_FILES(ch_input_sample)
        
        if (params.read_type == 'PE'){
            CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[2][0], 'R1']}.set{r1}
            CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[2][1], 'R2']}.set{r2}
            read_ch = r1.mix(r2)
        } else if (params.read_type == 'SE'){
            CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[2][0], 'R1']}.set{read_ch}
        }
        CONCATENATE_LOCAL_FILES.out.read_meta_ch.map{it -> [it[0], it[1]]}.set{meta_ch}
    }

    // Step 00: Concatenate Fastq files if required. 
    if (params.concat_lanes && !params.csv_input){
        if (params.read_type == 'PE'){
            CONCATENATE_READS_PE(read_ch)
            temp_read_ch = CONCATENATE_READS_PE.out.concat_fastq
            temp_read_ch.map{it -> [it[0], it[1][0], 'R1']}.set{r1}
            temp_read_ch.map{it -> [it[0], it[1][1], 'R2']}.set{r2}
            read_ch = r1.mix(r2)
        } else if (params.read_type == 'SE'){
            CONCATENATE_READS_SE(read_ch)
            temp_read_ch = CONCATENATE_READS_SE.out.concat_fastq
            temp_read_ch.map{it -> [it[0], it[1], 'R1']}.set{read_ch}
        }
    }

    RUN_EMASE(read_ch)
    // workflow found in: subworkflows/run-emase

    if (meta_ch) {
        RUN_EMASE.out.emase_genes_tpm.join(meta_ch)
        .map{it -> [it[0], it[1], it[2].sex, it[2].generation]}
        .set{grbs_recon_input}
    } else {
        RUN_EMASE.out.emase_genes_tpm
        .map{it -> [it[0], it[1], params.sample_sex, params.sample_generation]}
        .set{grbs_recon_input}              
    }

    GBRS_RECONSTRUCT(grbs_recon_input)

    GBRS_QUANTIFY_GENOTYPES(RUN_EMASE.out.compressed_emase_h5.join(GBRS_RECONSTRUCT.out.genotypes_tsv))

    GBRS_INTERPOLATE(GBRS_RECONSTRUCT.out.genoprobs_npz)

    GBRS_PLOT(GBRS_INTERPOLATE.out.interpolated_genoprobs)

    GBRS_EXPORT(GBRS_INTERPOLATE.out.interpolated_genoprobs)

}

For Diversity Outbred (DO) mice, the reference files are stored in the 
reference data directories `/projects/omics_share` on *sumner*. We will use 
these files as examples ot show you which arguments to provide. 

JAX uses [slurm](https://jacksonlaboratory.sharepoint.com/sites/ResearchIT/SitePages/Submitting-Basic-Jobs-with-SLURM.aspx)
(NOTE: This is an internal JAX link which requires a JAX login.) to manage 
computing jobs on *sumner*. There are several good tutorials on using *slurm*
on the JAX [Research IT Documentation Library](https://jacksonlaboratory.sharepoint.com/sites/ResearchIT/SitePages/Documentation.aspx).

#### *slurm* Options Block

*slurm* scripts are *bash* scripts, so they start with:

```
#!/bin/bash
```

Next, we use a *slurm* options block to tell *slurm* what computing resources
we would like to use. You could pass *slurm* options in at the command line, but
we recommend using an options block in your script so that you have a record of
the resources required to run your job. 

Each line in the *slurm* options block starts with:

```
#SBATCH
```

This tells *slurm* that what follows on this line is a *slurm* option. There 
are many *slurm* options, which are exhaustively enumerated in the 
[slurm documentation](https://slurm.schedmd.com/sbatch.html). Here we will use
only a few of the options:

```
#SBATCH --qos=batch       # job queue
#SBATCH --ntasks=1        # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
#SBATCH --mem=1G          # memory pool for all cores
#SBATCH --time=48:00:00   # wall clock time (D-HH:MM)
```

The first option specifies the job queue to use for this job. In this case, it
is the "batch" queue on *sumner*.

```
#SBATCH --qos=batch       # job queue
```

The next two options specify the number of nodes and number of CPUs (or cores) 
per node to use. In this case, we will use one node and one CPU. This may seem 
confusing because we would like to use many nodes and cores to reduce the total
run time of our job. But Nextflow handles resource request. This script just
starts the Nextflow job, which does not require many resources. If you request
more nodes and CPUs here, you will only be tying up resources that will sit 
idle while your job run. So **one** node and **one** CPU per node are sufficient.

```
#SBATCH --ntasks=1        # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
```

The next option specifies the memory required for this task. We will request 
one Gigabyte (G). Again, this may seem like a small amount of memory for a 
large computational job, but this is only the memeory required to start and
manage the Nextflow job. Nextflow will request the resources that it needs to 
run GBRS. 

```
#SBATCH --mem=1G          # memory pool for all cores
```

The last option that we will speicify is the total wall-clock time needed to 
complete the job. This needs to be the total amount of time required to complete
the job. This can be difficult to estimate and depends upon the number of 
samples, depth of coverage, and the configuration and utilization of the cluster.

> DMG: TBD: Should we run 50, 100, 200 & 300 samples and include a time/samples plot? 
I have 350 samples that we could use for this. It's not going to be linear with
more sample, but it's something.

```
#SBATCH --time=48:00:00   # wall clock time (D-HH:MM)
```

#### Header Section

This section is not required, but it allows you to add some comments about what
are doing. It is also a consistent place to put your name, contact information,
and the date. Putting your name and contact information helps to keep you 
accountable for the code. Putting the date in the header helps you to remember
when you wrote the script.

```
################################################################################
# Run GBRS on each sample.
# GRCm39. Ensembl 105.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-09-28
################################################################################
```

#### Variable Section

We strongly recommend that you use *bash* variables to build file paths and
specify arguments to the GBRS script. The reason for this is that it is easier
to remember what each varible is when it is given a meaningful name.

Consider this script which runs the "analysis_tool_of_science":

```
analysis_tool_of_science -t 8 \
  -g /path/to/another/file \
  -a /really/long/path/to/yet/another/file \
  -i /path/to/mystery/file/number1 /path/to/mystery/file/number2 \
  -x /path/to/some/directory/on/the/cluster
  -o /path/to/another/directory/with/a/file/prefix
```

You can read the tool's documentation to find out what each file is, but that
requires extra time. It would be easier to read if the paths had meaningful
names. In the code block below, we give each file a meaningful variable name
and prrovide a comment explaining what each file is. We also use this
opportunity to specify the genome and annotation builds.

```
# FASTA file containing genome sequence for GRCm39.
GENOME_FASTA=/path/to/another/file

# GTF file containing genome annotation for Ensembl 105.
ANNOTATION_GTF=/really/long/path/to/yet/another/file

# Paths to the forward and reverse FASTQ files.
FASTQ_FILE1=/path/to/mystery/file/number1 
FASTQ_FILE2=/path/to/mystery/file/number2

# Path to temporary directory for analysis_tool.
TEMP_DIR=/path/to/some/directory/on/the/cluster

# Path to output file, with a file prefix.
OUTPUT_PATH=/path/to/another/directory/with/a/file/prefix

analysis_tool_of_science -t 8 \
  -g ${GENOME_FASTA} \
  -a ${ANNOTATION_GTF} \
  -i $(FASTQ_FILE1} ${FASTQ_FILE2} \
  -x ${TEMP_DIR}
  -o ${OUTPUT_PATH}
```

This makes the script much easier for yourself (a year from now) and others
to read and understand.

The GBRS script has a lot of variables. We will explain each one so that you 
know what to use for your samples.

In this first section, we will focus on variables that set the location of
your files and working directories.

```
# Base directory for project.
BASE_DIR=/projects/research_lab_directory/diabetes
```
We like to set a "base directory", which is the high-level directory for the
project which you are working on. In this case, we are using a diabetes project
as an example. We use the base directory to create other file paths which
are below this directory.


```
# Sample metadata file.
SAMPLE_META_FILE=${BASE_DIR}/data/sodo_gbrs_metadata.csv
```

This is the path to the sample metadata file. This file contains information
about your samples, including the sample ID, sex, outbreeding generation,
sequencer lane, and paths to the forward and reverse FASTQ files. Each row will
contain information for one sample.

The `sampleID` should contain values that you use to identify the samples.

The `sex` column should contain the mouse's sex recorded as 'F' for females or 
'M' for males.

The `generation` column should contain the outbreeding generation for the mouse.

The `lane` column should contain the lane number on which the sample was run.
This can sometimes be found in the FASTQ file name. In the example below, it is
recorded as 'L004' in the FASTQ file names and this is what we used in the 
`lane` column.

The `fastq_1` and `fastq_2` columns contain the full path to the forward and 
reverse FASTQ files. In this case, we copied the files to /flashscratch, which
is the temporary sratch storage space for large files.


> DMG: TBD: explain why we should copy the files over to /flashscrath instead
of using them from /projects.

```
sampleID,sex,generation,lane,fastq_1,fastq_2
D2,F,46,L004,/flashscratch/dgatti/fastq/D2_GT22-11729_CTAGGCAT-ACAGAGGT_S224_L004_R1_001.fastq.gz,/flashscratch/dgatti/fastq/D2_GT22-11729_CTAGGCAT-ACAGAGGT_S224_L004_R2_001.fastq.gz
D3,F,46,L001,/flashscratch/dgatti/fastq/D3_GT22-11665_CGGCTAAT-CTCGTTCT_S21_L001_R1_001.fastq.gz,/flashscratch/dgatti/fastq/D3_GT22-11665_CGGCTAAT-CTCGTTCT_S21_L001_R2_001.fastq.gz
D4,F,46,L001,/flashscratch/dgatti/fastq/D4_GT22-11670_TACGCTAC-CGTGTGAT_S65_L001_R1_001.fastq.gz,/flashscratch/dgatti/fastq/D4_GT22-11670_TACGCTAC-CGTGTGAT_S65_L001_R2_001.fastq.gz
D5,F,46,L001,/flashscratch/dgatti/fastq/D5_GT22-11708_GAGCAGTA-TGAGCTGT_S57_L001_R1_001.fastq.gz,/flashscratch/dgatti/fastq/D5_GT22-11708_GAGCAGTA-TGAGCTGT_S57_L001_R2_001.fastq.gz
```

```
# Base /flashscratch directory.
FLASHSCRATCH=/flashscratch/username
```

This is the path to your directory on /flashscratch. You should create this 
before you begin your analysis.

```
# Temporary directory.
TMP_DIR=${FLASHSCRATCH}/tmp
```

This is the path to the temporary directory where Nextflow and GBRS can store
temporary working files during the analysis.

```
# Output directory.
OUT_DIR=${FLASHSCRATCH}/results
```

This is the directory where you would like GBRS to write your results. You 
should place this on /flashscratch and then copy the files that you need over
to /projects when the analysis is complete.


```
# Final results directory (in backed up space)
DEST_DIR=${BASE_DIR}/results/gbrs
```

This is the directory to which you will copy your final results. It should be
a location that is backed up (like /projects). Remember, /flashscratch is not
backed up and your files will be deleted without notice after 10 days.

In the next section, we will focus on the reference files that are required by
GBRS.

```
# GBRS Reference File Directory
GBRS_REF_DIR=/projects/omics_share/mouse/GRCm39/supporting_files/emase_gbrs/rel_2112_v8
```

This is the directory where the GBRS reference files are stored. These are a
series of files which contain information about genes and transcripts in the
eight DO founder strains. The numbers in the directory refere to the 
[Sanger Mouse Genomes Project](https://www.sanger.ac.uk/data/mouse-genomes-project/)
release. "2112" refers to the 12th month of 2021. "v8" refers to version 8 of
their genetic variant release.

There are eight arguments that point to the reference files used by GBRS.

```
# GBRS Transcripts.
GBRS_TRANSCRIPTS=${GBRS_REF_DIR}/emase.fullTranscripts.info
```

This file contains the list of transcripts in Ensembl 105 that GBRS will quantify.

```
# GBRS/EMASE gene to transcript file.
GBRS_GENE2TRANS=${GBRS_REF_DIR}/emase.gene2transcripts.tsv
```

This file contains the list of transcripts which map to each gene. Each row
contains one gene, followed by its associated Ensembl transcripts.

```
# GBRS Full Transcript.
GBRS_FULL_TRANSCRIPTS=${GBRS_REF_DIR}/emase.pooled.fullTranscripts.info
```

This file contains a list of transcript lengths in each of the eight DO founders.

> DMG: Confirm this with Mike.

```
# GBRS Emission Probs.
GBRS_EMIS_PROBS=${GBRS_REF_DIR}/gbrs_emissions_all_tissues.avecs.npz
```

This file contains the allele emission probabilities for the hidden Markov model.
These are used at each gene to estimate the probability that a given allele
is obsered, given that the model is in a certain genotype state.

```
# GBRS Transmission Probs.
GBRS_TRANS_PROBS=${GBRS_REF_DIR}/transition_probabilities
```

This directory contains the transmission probabilities for the hidden Markov model.
There are many files in this directory, one for each DO outbreeding generation.
Each file contains the probability of changing genotype state between each gene
at a given DO outbreeding generation.

```
# Ensembl 105 gene positions.
ENSEMBL_105=${GBRS_REF_DIR}/ref.gene_pos.ordered_ensBuild_105.npz
```

This file contains the Ensembl 105 gene positions in a python format.

```
# GBRS 69K Marker Grid.
MARKER_GRID=${GBRS_REF_DIR}/ref.genome_grid.GRCm39.tsv
```

This file contains the pseudomarker genome grid which is used to report results.
Each tissue expresses different gene, which leads GBRS to estimate genotypes
at different genomic positions in each tissue. In order to standardize results,
GBRS interpolates its founder haplotype estimates to a common grid of 74,165
positions.

```
# Bowtie index for GBRS.
BOWTIE_INDEX=/projects/compsci/omics_share/mouse/GRCm39/transcriptome/indices/imputed/rel_2112_v8/bowtie/bowtie.transcripts
```

This is the path to the bowtie index files. There are multiple files and this
variable should contain the directory and the file prefix for the index files.
In this case, the file prefix is "bowtie.transcripts" and the files are named
"bowtie.transcripts.1.ewbt", "bowtie.transcripts.2.ewbt", etc.

```
# GBRS path.
GBRS_GITHUB_PATH=TheJacksonLaboratory/cs-nf-pipelines
```

This is the Github path to the Computational Science `nextflow` pipelines.

```
# Singularity cache directory.
export NXF_SINGULARITY_CACHEDIR=${flashscratch}/singularity_cache
```

This is the cache directory where Singularity will store containers which it
downloads during the analysis run.

> DMG: confirm with Mike.

#### Entire GBRS Script

Here is the whole script in one piece.

```
#!/bin/bash
#SBATCH --qos=batch       # job queue 
#SBATCH --ntasks=1        # number of nodes
#SBATCH --cpus-per-task=1 # number of cores
#SBATCH --mem=1G          # memory pool for all cores
#SBATCH --time=48:00:00   # wall clock time (D-HH:MM)

################################################################################
# Run GBRS on each sample.
# GRCm39. Ensembl 105.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-09-28
################################################################################


##### VARIABLES #####

# Base directory for project.
BASE_DIR=/projects/bolcun-filas-lab/DO_Superovulation

# Sample metadata file.
SAMPLE_META_FILE=${BASE_DIR}/gbrs_sample_metadata.csv

# Base /flashscratch directory.
FLASHSCRATCH=/flashscratch/dgatti

# Temporary directory.
TMP_DIR=${FLASHSCRATCH}/tmp

# Output directory.
OUT_DIR=${FLASHSCRATCH}/results

# Final results directory (in backed up space)
DEST_DIR=${BASE_DIR}/results/gbrs/grcm39

# GBRS Reference File Directory
GBRS_REF_DIR=/projects/churchill-lab/projects/GBRS_GRCm39

# GBRS Transcripts.
GBRS_TRANSCRIPTS=${GBRS_REF_DIR}/emase.fullTranscripts.info

# GBRS/EMASE gene to transcript file.
GBRS_GENE2TRANS=${GBRS_REF_DIR}/emase.gene2transcripts.tsv

# GBRS Full Transcript.
GBRS_FULL_TRANSCRIPTS=${GBRS_REF_DIR}/emase.pooled.fullTranscripts.info

# GBRS Emission Probs.
GBRS_EMIS_PROBS=${GBRS_REF_DIR}/gbrs_emissions_all_tissues.avecs.npz

# GBRS Transmission Probs.
GBRS_TRANS_PROBS=${GBRS_REF_DIR}/transition_probabilities

# Ensembl 105 gene positions.
ENSEMBL_105=${GBRS_REF_DIR}/ref.gene_pos.ordered_ensBuild_105.npz

# GBRS 69K Marker Grid.
MARKER_GRID=${GBRS_REF_DIR}/ref.genome_grid.GRCm39.tsv

# Bowtie index for GBRS.
BOWTIE_INDEX=/projects/compsci/omics_share/mouse/GRCm39/transcriptome/indices/imputed/rel_2112_v8/bowtie/bowtie.transcripts

# GBRS path.
GBRS_GITHUB_PATH=TheJacksonLaboratory/cs-nf-pipelines

# Singularity cache directory.
export NXF_SINGULARITY_CACHEDIR=${flashscratch}/singularity_cache


##### MAIN #####

mkdir -p ${OUT_DIR}
mkdir -p ${DEST_DIR}
mkdir -p ${TMP_DIR}
mkdir -p ${NXF_SINGULARITY_CACHEDIR}

module load singularity

cd ${TMP_DIR}

nextflow run ${GBRS_GITHUB_PATH} \
         -profile sumner \
         --workflow gbrs \
         --pubdir ${OUT_DIR} \
         -w ${TMP_DIR} \
         --bowtie_index ${BOWTIE_INDEX} \
         --csv_input ${SAMPLE_META_FILE} \
         --transcripts_info ${GBRS_TRANSCRIPTS} \
         --gene2transcript_csv ${GBRS_GENE2TRANS} \
         --full_transcript_info ${GBRS_FULL_TRANSCRIPTS} \
         --emission_prob_avecs ${GBRS_EMIS_PROBS} \
         --trans_prob_dir ${GBRS_TRANS_PROBS} \
         --gene_position_file ${ENSEMBL_105} \
         --genotype_grid ${MARKER_GRID} \
         -resume

# Copy output files from OUT_DIR to DEST_DIR.
DIRS=`ls ${OUT_DIR}`

for i in ${DIRS}
do

  CURR_DEST_DIR=${DEST_DIR}/$i

  mkdir -p ${CURR_DEST_DIR}
  cp ${OUT_DIR}/${i}/stats/* ${CURR_DEST_DIR}
  cp ${OUT_DIR}/${i}/gbrs/*_counts ${CURR_DEST_DIR}
  cp ${OUT_DIR}/${i}/gbrs/*.tsv ${CURR_DEST_DIR}
  cp ${OUT_DIR}/${i}/gbrs/*.pdf ${CURR_DEST_DIR}

done
```

::::::::::::::::::::::::::::::::::::: keypoints 

- Use *bash* variables to build file paths and analysis tool arguments.

::::::::::::::::::::::::::::::::::::::::::::::::

[r-markdown]: https://rmarkdown.rstudio.com/
