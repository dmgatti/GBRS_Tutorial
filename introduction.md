---
title: "Genotyping by RNA Sequencing"
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is Genotyping by RNA Sequencing (GBRS)?
- What analyses does GBRS provide?
- What are the steps in a GBRS analysis?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Explain what GBRS does.
- Understand when and why you would use GBRS for transcript quantification.
- List the steps in a GBRS analysis.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

Genotyping by RNA Sequencing (**GBRS**) is a suite of tools which estimates 
allele-specific transcript expression and performs haplotype reconstruction
in multi-founder mouse crosses. GBRS can also by used to analyze transcript 
expression in other organisms, but this tutorial will focus on mouse crosses.
We have configured GBRS to use the [Nextflow](https://www.nextflow.io/) 
workflow software and [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)
containers.

## Prerequisites

This tutorial assumes that you have a background in genetics and are familiar
with bash command-line programming. 

In the area of genetics, we assume that:

1. you understand the folloiwng terms: inbred strain, recombinant inbred
strain, multi-founder advanced generation intercross;
1. you are familiar with the 
[Diversity Outbred](https://www.jax.org/strain/009376) (J:DO) mouse population;
1. you understand the concepts of haplotype reconstruction and allele-
specific expression;
1. you understand what RNASeq analysis is used for and what type of data it
produces.

In the area of bash programming, we assume that:

1. you can navigate the UNIX/Linux file system;
1. you understand the difference between absolute and relative file paths;
1. you can create and use bash variables;
1. you understand what a software container is and know how to use one.

## Background

There are many tool chains that can align RNASeq reads and quantify transcript 
expression. The alignment tools (i.e. [bowtie](https://bowtie-bio.sourceforge.net/index.shtml)
or [STAR](https://github.com/alexdobin/STAR)) generally align the reads to the
reference genome or transcriptome, allowing for gaps and mismatches. These 
methods work well for reads derived from mice with genomes that are similar to 
the reference, but may not provide accurate alignments for mouse genome that 
diverge from the reference genome. This may be particularly important when 
users need allele-specific expression estimates. Allele-specific expression 
provides an estimate of the expression of each of the two alleles in a diploid 
sample.

GBRS provides a suite of tools which aligns reads to pseudogenomes and
then quantifies allele-specific expression in multi-founder mouse crosses. A 
"**pseudo-genome**" is a reference genome with SNPs and indels inserted from 
the genome of a non-reference mouse strain. By analogy, we can also create a 
"**pseudo-transcriptome**" for a non-reference inbred strain, which involves 
inserting SNPs and indels into the reference transcriptome. We then align RNASeq
reads to multiple pseudo-transriptomes and use 
[expectation maximization](https://pubmed.ncbi.nlm.nih.gov/29444201/)
to estimate allele-specific transcript counts and total gene counts. In the 
process, we also reconstruct the haplotypes of each mouse in terms of the 
founder strain haplotypes.

The GBRS suite consists of three tools:

1. g2gtools: given a reference genome and a VCF containing SNPs and indels, 
creates pseudo-genomes by inserting the SNPs and indels into the reference
genome.
1. [EMASE](https://pubmed.ncbi.nlm.nih.gov/29444201/): given a set of 
pseudogenomes and a transcriptome, create a multi-way transcriptome. EMASE also
produces estimates of allele-specific transcript counts.
1. GBRS: Genotyping-by-RNA-Sequencing uses the pesudo-genomes created by 
g2gtools and the RNASeq read alignments from EMASE to reconstruct the 
haplotypes of each sample in terms of the founder strain haplotypes.

GBRS is designed to work well with bulk RNASeq off with sufficient read depth
to confidently call genetic variants. It does not work well with single-cell
RNASeq because of the low sequencing depth in each cell.


## Organization of this Tutorial

This tutorial is intended to cover all aspects of GBRS, from building reference
genomes and transcriptomes through haplotype reconstruction and estimation of
allele-specific transcript expression. All users may not need to perform all 
steps.

We envision four types of users:

1. Users who are analyzing [J:DO](https://www.jax.org/strain/009376) data
  - inside of The Jackson Laboratory
  This is the simplest case. The Next-Generation Sequencing Operations group
  has created the reference files and a pipeline to execute this workflow.
  This case is covered in the "Running GBRS" lesson.
  - outside of The Jackson Laboratory
  In this case, users will need to download and store the DO-related 
  reference files and have Nextflow installed on their computing cluster. 
  This case is covered in the "Running GBRS" lesson.
2. Users who are analyzing other mouse crosses
  - inside of The Jackson Laboratory
  Users will need to create reference files for the founder strains in their
  cross and point to them in their GBRS scripts. This is described in the 
  "Preparing GBRS Reference Genomes" and "Preparing GBRS Reference Transcriptomes"
  lessons.
  - outside of The Jackson Laboratory
  Users will need to create reference files for the founder strains in their
  cross and point to them in their GBRS scripts. This is described in the 
  "Preparing GBRS Reference Genomes" and "Preparing GBRS Reference Transcriptomes"
  lessons.
  

<!-- Reference files for DO:  https://zenodo.org/record/8186981  -->


::::::::::::::::::::::::::::::::::::: keypoints 

- GBRS is a suite of tools which includes pseudo-genome creation, allele-
specific transcript estimation, and haplotype reconstruction.
- Users who have J:DO data and work at The Jackson Laboratory can use
reference data stored on the high-performance computing cluster.
- Users outside of The Jackson Laboratory will need to create or download
the reference files.

::::::::::::::::::::::::::::::::::::::::::::::::

[r-markdown]: https://rmarkdown.rstudio.com/
