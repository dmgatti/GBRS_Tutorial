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
allele-specific transcript expression in multi-founder mouse crosses. GBRS can
also by used to analyze transcript expression in other organisms, but this 
tutorial will focus on mouse crosses.

## Prerequisites

This tutorial assumes that you have a background in genetics and are familiar
with bash command-line programming. 

In the area of genetics, we assume that:

1. you understand the folloiwng terms: inbred strain, recombinant inbred
strain, multi-founder advanced generation intercross;
1. you are familiar with the 
[Diversity Outbred](https://www.jax.org/strain/009376) mouse population;
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
users need allele-specific expression estimates.

GBRS provides a suite of tools which aligns reads to pseudogenomes and
then quantifies allele-specific expression in multi-foundr mouse crosses. A 
"**pseudo-genome**"" is a reference genome with SNPs and indels inserted from 
the genome of a non-reference mouse strain. By analogy, we can also create a 
"**pseudo-transcriptome**" for a non-reference inbred strain, which involves 
inserting SNPs and indels into the reference transcriptome. We then align RNASeq
reads to multiple pseudo-transriptomes and use 
[expectation maximization](https://pubmed.ncbi.nlm.nih.gov/29444201/)
to estimate allele-specific transcript and total gene counts. In the process,
we also reconsruct the haplotypes of each mouse in terms of the founder strain
haplotypes.

The GBRS suite consists of three tools:

1. g2gtools: given a reference genome and a VCF containing SNPs and indels, 
creates pseudo-genomes by inserting the SNPs and indels into the reference
genome.
1. [EMASE](https://pubmed.ncbi.nlm.nih.gov/29444201/): given a set of 
pseudogenomes and a transcriptome, create a multi-way transcriptome. EMASE also
produces estimates of allele-specific transcript counts.
1. GBRS: 

GBRS is designed to work well with bulk RNASeq off with sufficient read depth
to confidently call genetic variants. It does not work well with single-cell
RNASeq because of the low sequencing depth in each cell.

> DMG: What else do we need? Is this enough detail? Is it correct? 

## Organization of this Tutorial

This tutorial is intended to cover all aspects of GBRS, from building reference
genomes and transcriptomes through haplotype reconstruction and estimation of
allele-specific transcript expression. All users may not need to perform all 
steps.

We envision four types of users:

1. Users who are analyzing J:DO data
  - inside of The Jackson Laboratory
  - outside of The Jackson Laboratory
2. Users who are analyzing other mouse crosses
  - inside of The Jackson Laboratory
  - outside of The Jackson Laboratory





This is a lesson created via The Carpentries Workbench. It is written in
[Pandoc-flavored Markdown](https://pandoc.org/MANUAL.txt) for static files and
[R Markdown][r-markdown] for dynamic files that can render code into output. 
Please refer to the [Introduction to The Carpentries 
Workbench](https://carpentries.github.io/sandpaper-docs/) for full documentation.

What you need to know is that there are three sections required for a valid
Carpentries lesson:

 1. `questions` are displayed at the beginning of the episode to prime the
    learner for the content.
 2. `objectives` are the learning objectives for an episode displayed with
    the questions.
 3. `keypoints` are displayed at the end of the episode to reinforce the
    objectives.

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
