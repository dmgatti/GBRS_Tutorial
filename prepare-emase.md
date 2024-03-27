---
title: "Prepare EMASE Reference Files"
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- What does EMASE do?
- What input files are required for the EMASE reference?
- What output files are created by the prepare EMASE process?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Explain what analysis EMASE performs.
- Understand the files needed to prepare an EMASE reference.
- Understand the files produced by the Prepare EMASE process.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

[Expectation Maximization for Allele-Specific Expression](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6022640/)
(EMASE) is software which estimates allele-specific transcript expression in
genetically diverse populations. It was designed to be used with genetically
diverse mice whose genomes are comprised of two or more inbred strains. EMASE
uses the reference genome, strain-specific genomes, and GTF files to create
strain-specific transcriptomes.

Before running EMASE, we must prepare the reference files which EMASE needs. 
Prepare EMASE will take each founder strain genome and its corresponding GTF 
file and will create a [bowtie](https://bowtie-bio.sourceforge.net/index.shtml)
index and a pooled transcript list.

> DMG: Need more detailed explanation above.

On sumner2, the GRCm39 strain-specific genome FASTA files are in: /projects/compsci/omics_share/mouse/GRCm39/genome/sequence/imputed/rel_2112_v8.
These files contain the genome sequence of each strain with the strain-specific 
SNPs and Indels inserted. These were created by the "prepare-pseudo-reference"
process.

The corresponding GTFs, which were also created by the "prepare-pseudo-reference"
process, are in: /projects/compsci/omics_share/mouse/GRCm39/transcriptome/annotation/imputed/rel_2112_v8

We will create two variables for these paths:

```
GENOME_REF_DIR=/projects/compsci/omics_share/mouse/GRCm39/genome/sequence/imputed/rel_2112_v8

GTF_REF_DIR=/projects/compsci/omics_share/mouse/GRCm39/transcriptome/annotation/imputed/rel_2112_v8
```

The Prepare EMASE process requires the following files:

| File Type | Argument Name | File Path | Description | 
|-----------|---------------|-----------|-------------|
| Genome FASTA files | --genome_file_list | `${GENOME_REF_DIR}/A_J.39.fa,${GENOME_REF_DIR}/C57BL_6J.39.fa,${GENOME_REF_DIR}/129S1_SvImJ.39.fa,,${GENOME_REF_DIR}/NOD_ShiLtJ.39.fa,${GENOME_REF_DIR}/NZO_HlLtJ.39.fa,${GENOME_REF_DIR}/CAST_EiJ.39.fa,${GENOME_REF_DIR}/PWK_PhJ.39.fa,${GENOME_REF_DIR}/WSB_EiJ.39.fa` | Comma-separated list of FASTA files for each input strain |
| Transcript GTF files | --gtf_file_list | `${GTF_REF_DIR}/A_J.39_DroppedChromAppended.gtf,${GTF_REF_DIR}/Mus_musculus.GRCm39.105.filtered.gtf,${GTF_REF_DIR}/129S1_SvImJ.39_DroppedChromAppended.gtf,${GTF_REF_DIR}/NOD_ShiLtJ.39_DroppedChromAppended.gtf,${GTF_REF_DIR}/NZO_HlLtJ.39_DroppedChromAppended.gtf,${GTF_REF_DIR}/CAST_EiJ.39_DroppedChromAppended.gtf,${GTF_REF_DIR}/PWK_PhJ.39_DroppedChromAppended.gtf,${GTF_REF_DIR}/WSB_EiJ.39_DroppedChromAppended.gtf` | Comma-separated list of GTF files for each input strain |

We also need to provide a set of letters that EMASE will use as short 
abbreviations for each strain. Note that the order of the strains above will
match the order of the strains below.

| Argument Name | Value | Description |
|---------------|-------|-------------|
| --haplotype_list | A,B,C,D,E,F,G,H | Comma-separated list of letters to be used as an abbreviation for each strain | 


The EMASE software has the following arguments:

```
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────╮
│ --genome-file      -G      FILE     Genome files, can seperate files by "," or have multiple -G      |
|                                     [default: None] [required]                                       │
│ --haplotype-char   -s      TEXT     haplotype, either one per -h option, i.e. -h A -h B -h C, or a   |
|                                     shortcut -h A,B,C [default: None]                                │
│ --gtf-file         -g      FILE     Gene Annotation File files, can seperate files by "," or have    |
|                                     multiple -G [default: None]                                      │
│ --out_dir          -o      TEXT     Output folder to store results (default: the current working     |
|                                     directory) [default: None]                                       │
│ --save-g2tmap      -m               saves gene id to transcript id list in a tab-delimited text file │
│ --save-dbs         -d               save dbs                                                         │
│ --no-bowtie-index  -m               skips building bowtie index                                      │
│ --verbose          -v      INTEGER  specify multiple times for more verbose output [default: 0]      │
│ --help                              Show this message and exit.                                      │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

Below is an example of the Prepare EMASE command using these arguments:

```
emase prepare --genome-file ${genome_file_list} \
              --gtf-file ${gtf_file_list} \
              --haplotype-char ${haplotype_list} \
              --save-g2tmap 
              --out_dir ./ \
              --save-g2tmap \
              --no-bowtie-index
```

```
workflow PREPARE_EMASE {
    // Prepare emase reference, given list of genomes and gtf files. 
    EMASE_PREPARE_EMASE()
    BOWTIE_BUILD(EMASE_PREPARE_EMASE.out.pooled_transcript_fasta, 'bowtie.transcripts')
    // clean transcript lists to add transcripts absent from certain haplotypes.
    CLEAN_TRANSCRIPT_LISTS(EMASE_PREPARE_EMASE.out.pooled_transcript_info)
}
```

```
#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=gbrs_mouse
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=1G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow

# RUN PIPELINE
nextflow ../main.nf \
-profile sumner \
--workflow prepare_emase \
--pubdir "/flashscratch/${USER}/outputDir" \
-w /flashscratch/${USER}/outputDir/work \
--genome_file_list "/path/to/genome/A.fa,/path/to/genome/B.fa,..." \
--gtf_file_list "/path/to/gtf/A.gtf,/path/to/gtf/B.gtf,..." \
--haplotype_list "A,B,..." \
--comment "This script will run prepare_emase to generate multiway references based on default parameters"
```



[r-markdown]: https://rmarkdown.rstudio.com/
