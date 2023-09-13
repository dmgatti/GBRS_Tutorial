# Genotyping-by-RNA-Sequencing Tutorial

This is a tutorial for the Genotyping-by-RNA-Sequencing (GBRS) suite of tools.
These tools perform several tasks to assist with aligning RNA-Seq reads to 
genomes which differ from the reference. More specifically, GBRS is designed 
to align reads to multi-founder crosses between model organisms. This tutorial
focuses on mouse genomes, which are well annotated.

Broadly, these tasks are:

1. Insert SNPs and Indels into the reference genome for each strain of 
interest (i.e. founder strains of a cross) to create imputed genomes;
1. Create strain-specific transcriptomes from the imputed genomes;
1. Align RNS-Seq reads to all of the founder strain genomes;
1. Perform haplotype reconstruction in the cross using genetic variants
in the reads;
1. Estimate allele-specific gene and transcript abundance.

[workbench]: https://carpentries.github.io/sandpaper-docs/
