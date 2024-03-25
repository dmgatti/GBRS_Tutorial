---
title: "Running EMASE to Quantify Allele-Specific Transcript Expression"
teaching: 60
exercises: 20
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is Expectation Maximization for Allele-Specific Expression (EMASE)?
- What analyses does EMASE provide?
- What are the steps in an EMASE analysis?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Explain what EMASE does.
- Perofrm an EMASE analysis using Nextflow.

::::::::::::::::::::::::::::::::::::::::::::::::

workflow EMASE {
    // Step 0: Download data and concat Fastq files if needed. 
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
    // workflow run as subworkflow due to re-use in GBRS workflow. 

}

## Introduction

Expectation Maximization for Allele-Specific Expression (EMASE) is 
a method for estimating allele-specific transcript expression in organisms 
with diploid genomes. It was developed with The [Diversity Outbred](https://www.jax.org/strain/009376) 
mice in mind, but also works in other organisms. 

Briefly, in a mouse cross comprised of many different inbred founder strains,
RNASeq reads may map to multiple genes, mulitple isoforms, and/or to multiple
founder alleles (termed "multi-mapping reads"). Previous transcript estimation methods 
either discarded multi-mapping reads or treated different types of multi-mapping 
reads equivalently. EMASE uses a hierarchical approach to allocate 
read counts by first allocating reads among genes, then among isoforms, and finally
between alleles. This method improves transcript and isoform abundance estimates
over methods which handle multi-mapping reads differently. For full details,
please refer to 
[the publication](https://academic.oup.com/bioinformatics/article/34/13/2177/4850941).

## At the Jackson Laboratory (JAX)

The GBRS Nextflow pipeline is configured to run on sumner2, the High Performance Computing (HPC) cluster at JAX.

We will consider two scenarios:

    You are analyzing expression data from Diversity Outbred mice;
    You are analyzing expression data from some other cross;

### GBRS for Diversity Outbred Mice

The Next-Generation Operations (NGSOps) team has configured the default arrguments for the GBRS Nextflow pipeline to use the reference files for GRCm39 and Ensembl version 105. This means that the arguments that point GBRS to the locations of the founder genomes, founder transcriptomes, and aligner indices are already set.

Here is the entire script:

```
#!/bin/bash
#SBATCH --mail-user=first.last@jax.org
#SBATCH --job-name=emase_mouse
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
--workflow emase \
--pubdir "/flashscratch/${USER}/outputDir" \
-w /flashscratch/${USER}/outputDir/work \
--sample_folder <PATH_TO_YOUR_SEQUENCES> \
--comment "This script will run emase analysis on mouse samples"
```

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: instructor

Inline instructor notes can help inform instructors of timing challenges
associated with the lessons. They appear in the "Instructor View"

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1: Can you do it?

What is the output of this command?

```r
paste("This", "new", "lesson", "looks", "good")
```

:::::::::::::::::::::::: solution 

## Output
 
```output
[1] "This new lesson looks good"
```

:::::::::::::::::::::::::::::::::


## Challenge 2: how do you nest solutions within challenge blocks?

:::::::::::::::::::::::: solution 

You can add a line with at least three colons and a `solution` tag.

:::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::

## Figures

You can use standard markdown for static figures with the following syntax:

`![optional caption that appears below the figure](figure url){alt='alt text for
accessibility purposes'}`

![You belong in The Carpentries!](https://raw.githubusercontent.com/carpentries/logo/master/Badge_Carpentries.svg){alt='Blue Carpentries hex person logo with no text.'}

::::::::::::::::::::::::::::::::::::: callout

Callout sections can highlight information.

They are sometimes used to emphasise particularly important points
but are also used in some lessons to present "asides": 
content that is not central to the narrative of the lesson,
e.g. by providing the answer to a commonly-asked question.

::::::::::::::::::::::::::::::::::::::::::::::::


## Math

One of our episodes contains $\LaTeX$ equations when describing how to create
dynamic reports with {knitr}, so we now use mathjax to describe this:

`$\alpha = \dfrac{1}{(1 - \beta)^2}$` becomes: $\alpha = \dfrac{1}{(1 - \beta)^2}$

Cool, right?

::::::::::::::::::::::::::::::::::::: keypoints 

- Use `.md` files for episodes when you want static content
- Use `.Rmd` files for episodes when you need to generate output
- Run `sandpaper::check_lesson()` to identify any issues with your lesson
- Run `sandpaper::build_lesson()` to preview your lesson locally

::::::::::::::::::::::::::::::::::::::::::::::::

[r-markdown]: https://rmarkdown.rstudio.com/
