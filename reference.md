---
title: 'Reference'
---

## Glossary


BAM file: A binary alignment/map file which contains the positions in the 
reference genome where RNA-Seq reads align. This is a compressed version of the
SAM file.

base pair (bp): a pair of nucleotides on opposite strands of a chromosoe at a
specific genomic position. Often used in the context of genome postions
e.g. "chromosome 1 at 3,000 base pairs" would indicate a position 3,000 from
the start of chromosome 1.

FASTQ file: A text file containing the reads from an RNA-Seq experiment, along 
with base-level quality data for each read.

FastQC: ([FASTQ file Quality Control](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
is a program which performs a variety of quality control checks on FASTQ files.

IGV: The ([Integrative Genome Viewer](https://software.broadinstitute.org/software/igv/))
is a program which allows you to visualize the read alignments to the reference
genome. It is used after you have aligned your reads to a reference genome and
you have BAM files.

MultiQC: ([MultiQC]()) is a program which summarizes output from multiple 
FastQC files to make it easier to digest the results of many samples.

read: A DNA sequence which is derived from an RNA transcript in a biological
sample. These reads are often between 25 and 150 nucleotides long. 

RSEM: ([RNA-Seq by Expectation Maximization](https://github.com/deweylab/RSEM))
is a method of quantifying gene and transcript isoform abundance which handles
reads which align to multiple locations in the genome.

SAM file: A [Sequence Alignment/Map](https://samtools.github.io/hts-specs/SAMv1.pdf) 
file is a text file which contains the positions in the reference genome where
RNA-Seq reads align. One read may align to multiple positions. It also contains
the sequence of the read and its base-level quality scores.

STAR: ([Spliced Transcripts Alignment to a Reference](https://github.com/alexdobin/STAR)) 
A program which aligns RNA derived from transcripts to a reference genome in
a manner which accounts for reads which cross exon splice junctions. This is
referred to as a "splice-aware aligner."
 

