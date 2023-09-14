---
title: "Preparing GBRS Pseudo-References"
teaching: 60
exercises: 15
---

:::::::::::::::::::::::::::::::::::::: questions 

- Why do we need to align to a non-reference genome?
- What is a pseudo-reference genome?
- How do we create a pseudo-reference genome?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Construct a pseudo-reference genome for Diversity Outbred mice.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

[Genome-to-Genome Tools](https://github.com/churchill-lab/g2gtools), (g2gtools) 
is a suite of tools which incorporates SNPS and indels from mouse strains into
the reference genome. Inserting these variants may lead to changes in the length
of the genome and g2gtools is able to map positions between the reference
genome and the modified genome.


## Creating Pseudo-reference Genomes and Transcriptomes

The first step in creating the reference files for GBRS is to create a set of
genomes and transcriptomes with SNPs and Indels from other strains inserted
into the reference genome. 

### Setup

We will use the g2gtools container stored in the public reference area on sumner. 

```
G2GTOOLS=/projects/omics_share/meta/containers/quay.io-jaxcompsci-g2gtools-74926ad.img
```

We will also use a tool called [samtools](http://www.htslib.org/doc/samtools.html) 
to query and modify sequence and annotation files.

```
SAMTOOLS=/projects/omics_share/meta/containers/quay.io-biocontainers-samtools-1.14--hb421002_0.img
```

We need two input files and we will produce one output file per strain.

The first file is the reference genome in FASTA format. 

```
REF_FASTA=/projects/omics_share/mouse/GRCm39/genome/sequence/ensembl/v105/Mus_musculus.GRCm39.dna.primary_assembly.fa
```

The second file is a VCF containing the SNPs of one or more strains. This file
should contain the strain(s) that you have in your mouse cross.

```
VCF_INDELS=/projects/omics_share/mouse/GRCm39/genome/annotation/snps_indels/rel_2112_v8/mgp_REL2021_indels.vcf.gz
VCF_SNPS=/projects/omics_share/mouse/GRCm39/genome/annotation/snps_indels/rel_2112_v8/mgp_REL2021_snps.vcf.gz
```

Next, we need to pass in the name of the strain for which we are creating the
VCI file. In this case, we will use DBA/2J.

```
STRAIN=DBA_2J
```

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1: Viewing Strain Names in a VCF

There are 50 strains in the VCF SNP file. What are their names? The strain
names are stored in the VCF header. The `samtools` container also includes a
tool called `tabix`, which queries and indexes tab-delimited files. Run the 
`samtools` container and look at the tabix help page. Find a command 
that will print out the SNP VCF header and find the strain names.

```
singularity run ${SAMTOOLS} tabix --help
```

:::::::::::::::::::::::: solution 

## Output
  
Call `tabix` with the '--only-header' option and pass in the SNP VCF. In this
course, we are using long option names. You could also use '-H' to achieve the
same result.

```
$ singularity run ${SAMTOOLS} tabix --only-header ${VCF_SNPS}
```

```output
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.13+htslib-1.13
##bcftoolsCommand=mpileup -f Mus_musculus.GRCm39.dna.toplevel.fa.gz -b samples -g 10 -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/AD -E -Q 0 -pm 3 -F 0.25 -d 500
##reference=Mus_musculus.GRCm39.dna.toplevel.fa.gz
##contig=<ID=1,length=195154279>
##contig=<ID=2,length=181755017>
##contig=<ID=3,length=159745316>
##contig=<ID=4,length=156860686>
##contig=<ID=5,length=151758149>
##contig=<ID=6,length=149588044>
##contig=<ID=7,length=144995196>
##contig=<ID=8,length=130127694>
##contig=<ID=9,length=124359700>
##contig=<ID=10,length=130530862>
##contig=<ID=11,length=121973369>
##contig=<ID=12,length=120092757>
##contig=<ID=13,length=120883175>
##contig=<ID=14,length=125139656>
##contig=<ID=15,length=104073951>
##contig=<ID=16,length=98008968>
##contig=<ID=17,length=95294699>
##contig=<ID=18,length=90720763>
##contig=<ID=19,length=61420004>
##contig=<ID=X,length=169476592>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Total allelic depths (high-quality bases)">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=MinDP,Number=1,Type=Integer,Description="Minimum per-sample depth in this gVCF block">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred-scaled Genotype Quality">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##bcftools_callCommand=call -mAv -f GQ,GP -p 0.99; Date=Wed Aug 11 21:20:03 2021
##bcftools_normCommand=norm --fasta-ref Mus_musculus.GRCm39.dna.toplevel.fa.gz -m +indels; Date=Fri Aug 13 11:11:49 2021
##FORMAT=<ID=FI,Number=1,Type=Integer,Description="High confidence (1) or low confidence (0) based on soft filtering values">
##FILTER=<ID=LowQual,Description="Low quality variants">
##VEP="v104" time="2021-08-30 23:27:00" cache="mus_musculus/104_GRCm39" ensembl-funcgen=104.59ae779 ensembl-variation=104.6154f8b ensembl=104.1af1dce ensembl-io=104.1d3bb6e assembly="GRCm39" dbSNP="150" gencode="GENCODE M27" regbuild="1.0" sift="sift"
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|SIFT|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS">
##bcftools_viewVersion=1.13+htslib-1.13
##bcftools_viewCommand=view -i 'FORMAT/FI[*] = 1' mgp_REL2021_snps.vcf.gz; Date=Sat Dec 18 19:08:09 2021
##bcftools_annotateVersion=1.13+htslib-1.13
##bcftools_annotateCommand=annotate -x INFO/VDB,INFO/SGB,INFO/RPBZ,INFO/MQBZ,INFO/MQBZ,INFO/MQSBZ,INFO/BQBZ,INFO/SCBZ,INFO/FS,INFO/MQOF,INFO/AC,INFO/AN,FORMAT/SP,FORMAT/ADF,FORMAT/ADR,FORMAT/GP; Date=Sat Dec 18 19:08:09 2021
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##bcftools_viewCommand=view -a -Oz -o final_mgp_REL2021_snps.vcf.gz; Date=Sat Dec 18 19:08:09 2021
##bcftools_annotateCommand=annotate -x INFO/MQ0F -Oz -o final_mgp_REL2021_snps.vcf.gz mgp_REL2021_snps.vcf.gz; Date=Mon Dec 20 07:12:23 2021
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  129P2_OlaHsd    129S1_SvImJ    129S5SvEvBrd     A_J     AKR_J   B10.RIII        BALB_cByJ       BALB_cJ BTBR_T+_Itpr3tf_J       BUB_BnJC3H_HeH  C3H_HeJ C57BL_10J       C57BL_10SnJ     C57BL_6NJ       C57BR_cdJ       C57L_J  C58_J   CAST_EiJCBA_J   CE_J    CZECHII_EiJ     DBA_1J  DBA_2J  FVB_NJ  I_LnJ   JF1_MsJ KK_HiJ  LEWES_EiJ       LG_J   LP_J     MAMy_J  MOLF_EiJ        NOD_ShiLtJ      NON_LtJ NZB_B1NJ        NZO_HlLtJ       NZW_LacJ       PL_J     PWK_PhJ QSi3    QSi5    RF_J    RIIIS_J SEA_GnJ SJL_J   SM_J    SPRET_EiJ       ST_bJ   SWR_J  WSB_EiJ  ZALENDE_EiJ
```
  
The VCF header contains information about the file format version, chromosomes,
metadata, and tools used to create the VCF. The column names are on the last 
line of the VCF header and this line contains the strain names. If you had a 
cross comprised of other strains, you would query the VCF to find their names.

:::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::
  
We will work on the sumner /fastscratch area. Create a directory with your 
user ID. Then create a directory called 'gbrs'. For example:

```
mkdir -p /fastscratch/dgatti/gbrs/${STRAIN}
```

Then change into the directory that you just created.

```
cd /fastscratch/dgatti/gbrs
```

Load in the Singularity software.

```
module load singularity
```

### Use `g2gtools vcf2vci` to Gather Strain-Specific SNPs and Indels

In the first step, we insert indels from a specific strain into the reference
genome using the `vcf2vci` command. This is written out to a text file in a format called "VCI". 

The arguments are:
  
```
$ g2gtools vcf2vci --help
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --vcf       -i      FILE     VCF files can seperate files by "," or have multiple -i [default: None]  │
│                                 [required]                                                               │
│ *  --fasta     -f      FILE     Fasta file matching VCF information [default: None] [required]           │
│ *  --strain    -s      TEXT     Name of strain/sample (column in VCF file) [default: None] [required]    │
│    --output    -o      FILE     Name of output file [default: None]                                      │
│    --diploid   -d               Create diploid VCI file                                                  │
│    --keep                       Keep track of VCF lines that could not be converted to VCI file          │
│    --pass                       Use only VCF lines that have a PASS for the filter value                 │
│    --quality                    Filter on quality, FI=PASS                                               │
│    --no-bgzip                   DO NOT compress and index output                                         │
│    --verbose   -v      INTEGER  specify multiple times for more verbose output [default: 0]              │
│    --help                       Show this message and exit.                                              │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```
  
We will use a subset of the arguments, passing in the reference FASTA file, 
the SNPD and indel VCFs, the strain name to search for in the indel VCF, and the
output file path.
  
Inputs:
--fasta: path to the GRCm39 reference FASTA file
--vcf: path(s) to the Sanger VCFs for SNPs and indels. Both can be passed in.
  
Output:
--output: A VCI file, which is a custom file format, akin to a BED format, 
which lists the positions and sequences of indels. Note that, while we list
the output file extension as "vci", `g2gtools` will gzip and index the file, 
so we will need to add the '.gz' extension when we use the file in downstream
commands.
  
```
STRAIN_VCI=${STRAIN}/REF-to-${STRAIN}.vci
singularity run ${G2GTOOLS} g2gtools vcf2vci \
                            --fasta ${REF_FASTA} \
                            --vcf ${VCF_INDELS} \
                            --vcf ${VCF_SNPS} \
                            --strain ${STRAIN}  \
                            --output ${STRAIN_VCI}
```

### Use `g2gtools patch` to Insert Strain-Specific SNPs into Reference

Next, we insert SNPs from a specific strain into the reference genome using the
'patch' command.

The arguments are:
  
```
$ g2gtools patch --help
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────╮
│ *  --input    -i      FILE     Fasta file to extract from [default: None] [required]             │
│ *  --vci      -c      FILE     VCI File to use [default: None] [required]                        │
│    --output   -o      FILE     Name of output file [default: None]                               │
│    --bed      -b      FILE     BED file name [default: None]                                     │
│    --region   -r      TEXT     Region to extract in chromosome:start-end format [default: None]  │
│    --reverse                   Reverse the direction of VCI file                                 │
│    --bgzip                     compress and index output                                         │
│    --verbose  -v      INTEGER  specify multiple times for more verbose output [default: 0]       │
│    --help                      Show this message and exit.                                       │
╰─────────────────────────────────────────────────────────────────────────────────────────────────╯
```
  
We will use the following arguments:
  
Inputs:
--input: path to the GRCm39 reference FASTA file
--vci: path to the **gzipped** VCI file created by `vcf2vci`.
  
Output:
--output: path to patched FASTA file with SNPs inserted into the
reference genome sequence.
  
```
PATCHED_FASTA=${STRAIN}/${STRAIN}.patched.fa
singularity run ${G2GTOOLS} g2gtools patch \
                            --input ${REF_FASTA} \
                            --vci ${STRAIN_VCI}.gz \
                            --output ${PATCHED_FASTA}
```
  
### Use `g2gtools transform`  to Insert Strain-Specific Indels into Strain FASTA

The 'transform' function takes the strain-specific SNP-patched FASTA file,
inserts indels from the VCI file, and outputs a FASTA file. This FASTA file 
contains the inferred sequence of a specific strain.

Next, we insert SNPs from a specific strain into the reference genome using the
'patch' command.

The arguments are:
   
```
$ g2gtools transform --help
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────╮
│ *  --input    -i      FILE     Fasta file to extract from [default: None] [required]             │
│ *  --vci      -c      FILE     VCI File to use [default: None] [required]                        │
│    --output   -o      FILE     Name of output file [default: None]                               │
│    --bed      -b      FILE     BED file name [default: None]                                     │
│    --region   -r      TEXT     Region to extract in chromosome:start-end format [default: None]  │
│    --reverse                   Reverse the direction of VCI file                                 │
│    --bgzip                     compress and index output                                         │
│    --verbose  -v      INTEGER  specify multiple times for more verbose output [default: 0]       │
│    --help                      Show this message and exit.                                       │
╰─────────────────────────────────────────────────────────────────────────────────────────────────╯
```
  
We will use the following arguments:
  
Inputs:
--input: path to the strain-specific patched FASTA file created by `patch`.
--vci: path to the **gzipped** VCI file created by `vcf2vci`.
  
Output:
--output: path to transformed FASTA file with SNPs and indels inserted into
the reference genome sequence.
  
```
STRAIN_FASTA=${STRAIN}/${STRAIN}.${GENOME_VERSION}.fa
singularity run ${G2GTOOLS} g2gtools transform \
                            --input ${PATCHED_FASTA} \
                            --vci ${STRAIN_VCI}.gz \
                            --output ${STRAIN_FASTA}
```

We now have a pseudo-reference genome in the strain-specific FASTA file.

### Use `samtools faidx` to Index the Strain-Specific FASTA file

This is a step which uses the [samtools](https://www.htslib.org/doc/samtools.html) 
suite of tools to index the FASTA file.

Remember that we created a variable for the `samtools` Singularity container 
at the top of the lesson.



Then we will use `samtools` to index the FASTA file.

```
singularity run ${SAMTOOLS} samtools faidx ${STRAIN_FASTA}
```

### Use `g2gtools convert` to Create a Strain-Specific GTF File

 The `convert` function takes the strain-specific FASTA file, the VCI file, 
 and the reference GTF and creates a strain-specific annotation file in GTF.

The arguments are:
  
```
$ g2gtools convert --help
╭─ Options ────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *  --input-file   -i      FILE                   Input file to convert to new coordinates [default: None] [required]  │
│ *  --vci-file     -c      FILE                   VCI file [default: None] [required]                                  │
│ *  --file-format  -f      [BAM|SAM|GFF|GTF|BED]  Input file format [default: None] [required]                         │
│    --output       -o      FILE                   Name of output file [default: None]                                  │
│    --reverse      -r                             Reverse the direction of the conversion                              │
│    --verbose      -v      INTEGER                specify multiple times for more verbose output [default: 0]          │
│    --help                                        Show this message and exit.                                          │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```
  
We will use the following arguments:
  
Inputs:
--input-file: path to the reference GTF.
--vci-file: path to the **gzipped** VCI file created by `vcf2vci`.
--file-format: string indicating a GTF file. Lower case since it will be used
as a file name extension.
  
Output:
--output: path to converted GTF file with SNPs and indels inserted into
the reference transcript and coordinates shifted to the strain-specific 
coordinates.
  
```
singularity run ${G2GTOOLS} 
```

Next, we need the path to the reference GTF. We have selected the Ensembl 105
annotation for this lesson.

```
REF_GTF=/projects/omics_share/mouse/GRCm39/transcriptome/annotation/ensembl/v105/Mus_musculus.GRCm39.105.chr.gtf
```

We also need to specify the annotation file format in lower case because it
will be used as a file name extension.

```
ANNOT_FORMAT=gtf
```

```
STRAIN_GTF=${STRAIN}/${STRAIN}.${GENOME_VERSION}.${ANNOT_FORMAT}
singularity run ${G2GTOOLS} g2gtools convert \
                            --input-file ${REF_GTF} \
                            --vci-file ${STRAIN_VCI}.gz \
                            --file-format ${ANNOT_FORMAT} \
                            --output ${STRAIN_GTF}
```

We now have two key files that are used by GBRS:

- the strain-specific FASTA file, which contains the strain's inferred sequence,
and
- the strain-specific GTF file, which contains the strain's inferred annotation.


::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1: Create a set of reference files for another strain in the VCF.

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
::::::::::::::::::::::::::::::::::::::::::::::::


```
g2gtools transform -i ${STRAIN}/${STRAIN}.patched.fa -c ${STRAIN}/REF-to-${STRAIN}.chain -o ${STRAIN}/${STRAIN}.fa
g2gtools convert   -c ${STRAIN}/REF-to-${STRAIN}.chain -i ${GTF} -f gtf -o ${STRAIN}/${STRAIN}.gtf

g2gtools gtf2db                -i ${STRAIN}/${STRAIN}.gtf -o ${STRAIN}/${STRAIN}.gtf.db
g2gtools extract --transcripts -i ${STRAIN}/${STRAIN}.fa -db ${STRAIN}/${STRAIN}.gtf.db > ${STRAIN}/${STRAIN}.transcripts.fa
g2gtools extract --genes       -i ${STRAIN}/${STRAIN}.fa -db ${STRAIN}/${STRAIN}.gtf.db > ${STRAIN}/${STRAIN}.genes.fa
g2gtools extract --exons       -i ${STRAIN}/${STRAIN}.fa -db ${STRAIN}/${STRAIN}.gtf.db > ${STRAIN}/${STRAIN}.exons.fa
```

```
 g2gtools --help

 Usage: g2gtools [OPTIONS] COMMAND [ARGS]...

 g2gtools

╭─ Options ─────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --version                                                                                                 │
│ --install-completion        [bash|zsh|fish|powershell|pwsh]  Install completion for the specified shell.  │
│                                                              [default: None]                              │
│ --show-completion           [bash|zsh|fish|powershell|pwsh]  Show completion for the specified shell, to  │
│                                                              copy it or customize the installation.       │
│                                                              [default: None]                              │
│ --help                                                       Show this message and exit.                  │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ────────────────────────────────────────────────────────────────────────────────────────────────╮
│ convert           Convert coordinates of BAM|SAM|GTF|GFF|BED files                                        │
│ extract           Extract subsequence from a fasta file given a region                                    │
│ fasta-format      Reformat a Fasta file                                                                   │
│ gtf2db            Convert a GTF file to a G2G DB file                                                     │
│ parse-region      Parse a region from the command line.  Useful for verifying regions.                    │
│ patch             Patch SNPs onto the reference sequence                                                  │
│ transform         Incorporate indels onto the input sequence                                              │
│ vcf2vci           Create VCI file from VCF file(s)                                                        │
│ vci-query         Query a VCI file.                                                                       │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

For each strain the following is produced: 

```
A_J.39.exons.fa
A_J.39.fa
A_J.39.fa.fai
A_J.39.genes.fa
A_J.39.gtf
A_J.39.gtf.db
A_J.39.gtf.unmapped
A_J.39.patched.fa
A_J.39.patched.fa.fai
A_J.39.transcripts.fa
A_J.39.vci.gz
A_J.39.vci.gz.tbi
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
