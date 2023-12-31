---
title: "Running EMASE to Quantify Allele-Specific Transcript Expression"
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
- List the steps in a GBRS analysis.

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
