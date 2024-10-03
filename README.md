# 🧬 Lab: Next-Generation Sequencing Analysis 1 & 2 🧬

In this lab, we will guide you through the essential steps involved in processing high-quality reads from Next-Generation Sequencing (NGS) and mapping them to a reference genome. You will learn how to perform quality control on short-read NGS data, map high-quality reads, and understand the post-mapping processes necessary for effective variant calling and annotation. 

By the end of this tutorial, you will have a comprehensive understanding of the entire workflow, including:
- Data preparation
- Quality control
- Read mapping
- Generation of variant databases

## Learning Outcomes

### NGS Analysis 1
- Review the structure and origin of sequencing reads.
- Outline the file formats in which reads and mapped reads are stored.
- Describe the basic quality control steps necessary for processing read data.
- Describe the process of read mapping.

### NGS Analysis 2
- Understand the motivation to identify variants from NGS data in a clinical setting.
- Describe the GATK best practices for variant identification.
- Describe the file formats in which variants are stored.
- Outline the primary measures of variant quality control.

---

# 🛠 Before You Begin

1. **Create a Galaxy user account** using your RCSI email: [Sign up here](https://usegalaxy.org/login/start?redirect=None).
2. **Download the "data_files" folder** from [this link](https://rcsicampus-my.sharepoint.com/:f:/g/personal/laurawhelan_rcsi_com/EkI2pyMKZNxOjeDJOtqnB9EB3L5pV0j_TEIivBL5suTB7A?e=x4cKab).
3. **Create a new "history"** within Galaxy named "ATT_NGS_LAB".

![Create New History](https://github.com/user-attachments/assets/4d6c3652-22f2-4612-ad62-ae78b5c13c4b)

> **Note:** I have "pre-made" all the files for you. You're going to perform all the steps to make these files, but some of these steps take a long time computationally. That's why we have pre-made files ready for you, similar to a cooking show — *here’s one we made earlier!*

---

# 🧬 Lab: Next-Generation Sequencing Analysis 1 🧬
 

# Exome sequencing data analysis for diagnosing a genetic disease


> Overview
> 
> *   How do you identify genetic variants in samples based on exome sequencing data?
>     
> *   How do you, among the set of detected variants, identify candidate causative variants for a given phenotype/disease?
>     
> 
> **Objectives:**
> 
> *   Jointly call variants and genotypes for a family trio from whole-exome sequencing data
>     
> *   Use variant annotation and the observed inheritance pattern of a phenotype to identify candidate causative variants and to prioritize them
>     
>

Exome sequencing is a method that enables the selective sequencing of the exonic regions of a genome - that is the transcribed parts of the genome present in mature mRNA, including protein-coding sequences, but also untranslated regions (UTRs).

In humans, there are about 180,000 exons with a combined length of ~ 30 million base pairs (30 Mb). Thus, the exome represents only 1% of the human genome, but has been estimated to harbor up to 85% of all disease-causing variants ([Choi et al., 2009](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2768590/)).

Exome sequencing, thus, offers an affordable alternative to whole-genome sequencing in the diagnosis of genetic disease, while still covering far more potential disease-causing variant sites than genotyping arrays. This is of special relevance in the case of rare genetic diseases, for which the causative variants may occur at too low a frequency in the human population to be included on genotyping arrays.

Of note, a recent study focusing on the area of clinical pediatric neurology indicates that the costs of exome sequencing may actually not be higher even today than the costs of conventional genetic testing ([Vissers et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5589982/)).

> Details: Exome sequencing _vs_ whole-genome sequencing
> 
> In principle, the steps illustrated in this tutorial are suitable also for the analysis of whole-genome sequencing (WGS) data. At comparable mean coverage, however, WGS datasets will be much larger than exome sequencing ones and their analysis will take correspondingly more time.
> 
> The obvious benefit of WGS compared to exome-sequencing, of course, is that it will allow variant detection in even more regions of the genome. As a less apparent advantage, the more complete information of WGS data can make it easier to detect _copy number variation (CNV)_ and _structural variants_ such as translocations and inversions (although such detection will require more sophisticated analysis steps, which are not covered by this tutorial).
> 
> Very generally, one could argue that exome-sequencing captures most of the information that can be analyzed with standard bioinformatical tools today at reasonable costs. WGS, on the other hand, captures as much information as today’s sequencing technology can provide, and it may be possible to reanalyze such data with more powerful bioinformatical software in the future to exploit aspects of the information that were not amenable to analysis at the time of data acquisition.

The identification of causative variants underlying any particular genetic disease is, as we will see in this tutorial, not just dependent on the successful detection of variants in the genome of the patient, but also on variant comparison between the patient and selected relatives. Most often family trio data, consisting of the genome sequences of the patient and their parents, is used for this purpose. With multisample data like this it becomes possible to search for variants following any kind of Mendelian inheritance scheme compatible with the observed inheritance pattern of the disease, or to detect possibly causative _de-novo_ mutations or _loss-of-heterozygosity_ (LOH) events.



> Agenda
> 
> In this tutorial, we will cover:
> 
> 1.  [Data Preparation](#data-preparation)
>     1.  [Get data](#get-data)
> 2.  [Quality control](#quality-control)
> 3.  [Read mapping](#read-mapping)
> 4.  [Mapped reads postprocessing](#mapped-reads-postprocessing)
>     1.  [Filtering on mapped reads properties](#filtering-on-mapped-reads-properties)
>     2.  [Removing duplicate reads](#removing-duplicate-reads)
> 5.  [Variant calling](#variant-calling)
>     1.  [Generating FreeBayes calls](#generating-freebayes-calls)
>     2.  [Post-processing FreeBayes calls](#post-processing-freebayes-calls)
> 6.  [Variant annotation and reporting](#variant-annotation-and-reporting)
>     1.  [Get data](#get-data-1)
>     2.  [Variant annotation with functional genomic effects](#variant-annotation-with-functional-genomic-effects)
>     3.  [Generating a GEMINI database of variants for further annotation and efficient variant queries](#generating-a-gemini-database-of-variants-for-further-annotation-and-efficient-variant-queries)
>     4.  [Candidate variant detection](#candidate-variant-detection)
> 7.  [Conclusion](#conclusion)
>

# Let's Get Started
> **There will be a video to follow on screen!**

# Data Preparation

In this tutorial, we are going to analyze exome sequencing data from a family trio, in which the boy child is affected by the disease [osteopetrosis](https://ghr.nlm.nih.gov/condition/osteopetrosis/), while both parents, who happen to be consanguineous, are unaffected. Our goal is to identify the genetic variation that is responsible for the disease.

## Get data

tip This tutorial offers **two alternative entry points** allowing you to

*   conduct a full analysis starting from original sequenced reads in `fastq` format or to
*   start your analysis with premapped reads in `bam` format that are (almost) ready for variant calling.

The following hands-on section will guide you through obtaining the right data for either analysis.

> Hands-on: Data upload
> 
> 1.  Create a new history for this tutorial and give it a meaningful name
>     
>     > Tip: Creating a new history
>     > 
>     > To create a new history simply click the new-history icon at the top of the history panel:
>     > 
>     > ![UI for creating new history](/training-material/shared/images/history_create_new.svg)
>     
>     > Tip: Renaming a history
>     > 
>     > 1.  Click on galaxy-pencil (**Edit**) next to the history name (which by default is “Unnamed history”)
>     > 2.  Type the new name
>     > 3.  Click on **Save**
>     > 4.  To cancel renaming, click the galaxy-undo “Cancel” button
>     > 
>     > If you do not have the galaxy-pencil (**Edit**) next to the history name (which can be the case if you are using an older version of Galaxy) do the following:
>     > 
>     > 1.  Click on **Unnamed history** (or the current name of the history) (**Click to rename history**) at the top of your history panel
>     > 2.  Type the new name
>     > 3.  Press Enter
>     
> 2.  Obtain the raw sequencing data
>     
> 3.  Check that the newly created datasets in your history have their datatypes assigned correctly to `fastqsanger.gz`, and fix any missing or wrong datatype assignment
>     
>     > Tip: Changing the datatype
>     > 
>     > *   Click on the galaxy-pencil **pencil icon** for the dataset to edit its attributes
>     > *   In the central panel, click galaxy-chart-select-data **Datatypes** tab on the top
>     > *   In the galaxy-chart-select-data **Assign Datatype**, select `fastqsanger.gz` from “_New type_” dropdown
>     >     *   Tip: you can start typing the datatype into the field to filter the dropdown menu
>     > *   Click the **Save** button
>     
>    Congratulations for obtaining the datasets required for an analysis including reads mapping. You should now **proceed with Step 7** below.
>     
> 4.  Obtain the premapped sequencing data
>     
> 5.  Check that the newly created datasets in your history have their datatypes assigned correctly to `bam`, and fix any missing or wrong datatype assignment
>     
>     > Tip: Changing the datatype
>     > 
>     > *   Click on the galaxy-pencil **pencil icon** for the dataset to edit its attributes
>     > *   In the central panel, click galaxy-chart-select-data **Datatypes** tab on the top
>     > *   In the galaxy-chart-select-data **Assign Datatype**, select `bam` from “_New type_” dropdown
>     >     *   Tip: you can start typing the datatype into the field to filter the dropdown menu
>     > *   Click the **Save** button
>     
> 6.  Specify the genome version that was used for mapping
>     
>     Change the database/build (dbkey) for each of your bam datasets to `hg19`.
>     
>     > Tip: Changing database/build (dbkey)
>     > 
>     > *   Click the desired dataset’s name to expand it.
>     > *   Click on the “?” next to database indicator:
>     >     
>     >     ![UI for changing dbkey](/training-material/shared/images/datasets_dbkey.svg)
>     >     
>     > *   In the central panel, change the **Database/Build** field
>     > *   Select your desired database key from the dropdown list: `Human Feb. 2009 (GRCh37/hg19) (hg19)`
>     > *   Click the **Save** button

>     
>     Congratulations for obtaining the premapped sequencing datasets. Now, **follow the remaining steps** to set everything up for a successful analysis.
>     
> 7.  Rename the datasets [NOT NEEDED IN THIS TUTORIAL!]
>     
>     For datasets that you upload via a link, Galaxy will pick the link address as the dataset name, which you will likely want to shorten to just the file names.
>     
>     > Tip: Renaming a dataset
>     > 
>     > *   Click on the galaxy-pencil **pencil icon** for the dataset to edit its attributes
>     > *   In the central panel, change the **Name** field
>     > *   Click the **Save** button
>     
> 8.  Add #father/#mother/#child tags to the datasets
>     
>     Parts of the analysis in this tutorial will consist of identical steps performed on the data of each family member.
>     
>     To make it easier to keep track of which dataset represents which step in the analysis of which sample, Galaxy supports dataset tags. In particular, if you attach a tag starting with `#` to any dataset, that tag will automatically propagate to any new dataset derived from the tagged dataset.
>     
>     Tags are supposed to help you identify the origin of datasets quickly, but you can choose them as you like.
>     
>     > Tip: Adding a tag
>     > 
>     > Datasets can be tagged. This simplifies the tracking of datasets across the Galaxy interface. Tags can contain any combination of letters or numbers but cannot contain spaces.
>     > 
>     > **To tag a dataset**:
>     > 
>     > 1.  Click on the dataset to expand it
>     > 2.  Click on **Add Tags** galaxy-tags
>     > 3.  Add tag text. Tags starting with `#` will be automatically propagated to the outputs of tools using this dataset (see below).
>     > 4.  Press Enter
>     > 5.  Check that the tag appears below the dataset name
>     
> 9.  Obtain the reference genome
>     
>     
>     Import the `hg19` version of the human chromosome 8 sequence:
>     
>     ```
>     https://zenodo.org/record/3243160/files/hg19_chr8.fa.gz
>     ```
>     
>     In the upload dialog, make sure you specify:
>     
>     *   **Type**: `fasta`
>     *   **Genome**: `Human Feb. 2009 (GRCh37/hg19) (hg19)`
>     
>     Alternatively, load the dataset from a shared data library.
>     
> 10.  Rename the reference genome
>     
>     The reference genome you have imported above came as a compressed file, but got unpacked by Galaxy to plain `fasta` format according to your datatype selection. At a minimum, you may now wish to remove the `.gz` suffix from the dataset name to avoid confusion.
>     

**Congratulations!** You are all set for starting the analysis now.

If you have chosen to follow the complete analysis from the original sequenced data, just proceed with the next section.

If, on the other hand, you have prepared to start from the premapped data, skip the sections on _Quality control_ and _Read mapping_, and conitnue with **Mapped reads postprocessing**.

# Quality control

This step serves the purpose of identifying possible issues with the raw sequenced reads input data before embarking on any “real” analysis steps.

Some of the typical problems with NGS data can be mitigated by preprocessing affected sequencing reads before trying to map them to the reference genome. Detecting some other, more severe problems early on may at least save you a lot of time spent on analyzing low-quality data that is not worth the effort.


> Hands-on: Quality control of the input datasets
> 
> 1.  Run **FastQC** ( Galaxy version 0.74+galaxy0) on each of your six fastq datasets
>     
>     *   param-files _“Short read data from your current history”_: all 6 FASTQ datasets selected with **Multiple datasets**
>     
>     > Tip: Select multiple datasets
>     > 
>     > 1.  Click on param-files **Multiple datasets**
>     > 2.  Select several files by keeping the Ctrl (or COMMAND) key pressed and clicking on the files of interest
>     
>     When you start this job, twelve new datasets (one with the calculated raw data, another one with an html report of the findings for each input dataset) will get added to your history.
>     
> 2.  Use **MultiQC** ( Galaxy version 1.11+galaxy1) to aggregate the raw **FastQC** data of all input datasets into one report
>     *   In _“Results”_
>         *   _“Which tool was used generate logs?”_: `FastQC`
>         *   In _“FastQC output”_
>             *   _“Type of FastQC output?”_: `Raw data`
>             *   param-files _“FastQC output”_: all six _RawData_ outputs of **FastQC** tool)
> 3.  Inspect the _Webpage_ output produced by the tool
>     
>     > Question
>     > 
>     > 1.  Based on the report, do you think preprocessing of the reads (trimming and/or filtering) will be necessary before mapping?
>     > 2.  Why do all samples show a non-normal GC content distribution, and should you be worried?


# Read mapping

Now that you confirmed that the quality of the input data is good enough to warrant further analysis, it is time to map the sequenced reads to the reference genome.

> Hands-on: Read Mapping
> 
> 1.  **Map with BWA-MEM** ( Galaxy version 0.7.17.2) to map the reads from the **father** sample to the reference genome
>     
>     *   _“Will you select a reference genome from your history or use a built-in index?”_: `Use a built-in genome index`
>         
>         *   _“Using reference genome”_: `Human: hg19` (or a similarly named option)
>         
>         > Comment: Using the imported \`hg19\` sequence
>         > 
>         > If you have imported the `hg19` chr8 sequence as a fasta dataset into your history instead:
>         > 
>         > *   _“Will you select a reference genome from your history or use a built-in index?”_: `Use a genome from history and build index`
>         >     *   param-file _“Use the following dataset as the reference sequence”_: your imported `hg19` fasta dataset.
>         
>     *   _“Single or Paired-end reads”_: `Paired`
>         
>         *   param-file _“Select first set of reads”_: the forward reads (R1) dataset of the **father** sample
>         *   param-file _“Select second set of reads”_: the reverse reads (R2) dataset of the **father** sample
>         
>         > Tip: No FASTQ datasets selectable?
>         > 
>         > Please confirm that the problematic datasets declare _format_: `fastqsanger.gz`.
>         > 
>         > Most Galaxy tools that accept FASTQ input expect the data to be formatted as _FASTQ with Sanger-scaled quality values_, the most widely spread version of the FASTQ format. To make this requirement explicit (instead of generating possibly wrong results) these tools require you to set the dataset type to `fastqsanger` (`fastqsanger.gz` for data compressed with gzip). You can do so either on data upload or later from the _Edit dataset attributes_ view (which you can reach by clicking on the galaxy-pencil pencil icon.
>         
>     *   _“Set read groups information?”_: `Set read groups (SAM/BAM specification)`
>         *   _“Auto-assign”_: `No`
>             *   _“Read group identifier (ID)”_: `000`
>         *   _“Auto-assign”_: `No`
>             *   _“Read group sample name (SM)”_: `father`
>     
>     > Warning: Read group IDs and sample names - choose, but choose wisely
>     > 
>     > In general, you are free to choose ID and SM values to your liking, but …
>     > 
>     > The **ID** should **unambiguously identify** the sequencing run that produced the reads. At the very least, no two input datasets in any given analysis should define the same ID twice, or tools like _FreeBayes_, which we are going to use in the next step, will refuse to work with the data.
>     > 
>     > The **SM** value, on the other hand, should identify the biological sample represented by the data and is used by many tools (like _GEMINI_ which we will use later) to let you refer to one specifc sample in a multisample analysis. Choose descriptive, but short and easy to remember sample names since you will have to type them in again!
>     
> 2.  **Map with BWA-MEM** ( Galaxy version 0.7.17.2) to map the reads from the **mother** sample to the reference genome **using the same parameters as before** except
>     
>     *   _“Single or Paired-end reads”_: `Paired`
>         *   param-file _“Select first set of reads”_: the forward reads (R1) dataset of the **mother** sample
>         *   param-file _“Select second set of reads”_: the reverse reads (R2) dataset of the **mother** sample
>     *   _“Set read groups information?”_: `Set read groups (SAM/BAM specification)`
>         *   _“Auto-assign”_: `No`
>             *   _“Read group identifier (ID)”_: `001`
>         *   _“Auto-assign”_: `No`
>             *   _“Read group sample name (SM)”_: `mother`
> 3.  **Map with BWA-MEM** ( Galaxy version 0.7.17.2) to map the reads from the **child** sample to the reference genome **using the same parameters as before** except
>     
>     *   _“Single or Paired-end reads”_: `Paired`
>         *   param-file _“Select first set of reads”_: the forward reads (R1) dataset of the **child** sample
>         *   param-file _“Select second set of reads”_: the reverse reads (R2) dataset of the **child** sample
>     *   _“Set read groups information?”_: `Set read groups (SAM/BAM specification)`
>         *   _“Auto-assign”_: `No`
>             *   _“Read group identifier (ID)”_: `002`
>         *   _“Auto-assign”_: `No`
>             *   _“Read group sample name (SM)”_: `proband`

# Mapped reads postprocessing

At this point in the analysis you should have obtained three mapped reads datasets in `bam` format. Each of these datasets should:

*   have its _database_ set to the key `hg19`
    
    Please correct any missing (`?`) or wrong keys now!
    
    > Tip: Changing database/build (dbkey)
    > 
    > *   Click the desired dataset’s name to expand it.
    > *   Click on the “?” next to database indicator:
    >     
    >     ![UI for changing dbkey](/training-material/shared/images/datasets_dbkey.svg)
    >     
    > *   In the central panel, change the **Database/Build** field
    > *   Select your desired database key from the dropdown list: `Human Feb. 2009 (GRCh37/hg19) (hg19)`
    > *   Click the **Save** button
    
*   ideally, carry one of the `#father`, `#mother` or `#child` tags for quick identification of the samples they provide data for.
    

In principle, you could use these datasets directly for variant calling, and in many cases, including this one, this would be sufficient to identify the sought-after variants.

To obtain an accurate picture of the variant spectrum found in your samples it is good practice though to perform various postprocessing steps on the mapped reads before passing them to a variant caller.

The optimal set of postprocessing steps required depends on the variant calling software used at the next step. The **FreeBayes** variant caller that we are going to use in this tutorial is particularly well suited for use with minimal mapped reads postprocessing pipelines, so all we are going to do here is:

*   filter the paired-end reads of all samples to retain only those read pairs, for which both the forward and the reverse read have been mapped to the reference successfully
    
    For such pairs of reads, we can be extra confident that they don’t come from some non-human contaminant DNA or represent a sequencing artefact of some sort.
    
*   deduplicate reads
    
    Duplicate reads, which typically arise from PCR-overamplification of genomic fragments during sequencing library preparation, can, to some extent, lead to wrong genotype assignments at variant sites (if, for example, a sample is heterozygous for a variant, but fragments with one of the two alleles get amplified more efficiently than the others).
    

## Filtering on mapped reads properties

To produce new filtered BAM datasets with only mapped reads the mate of which is also mapped:

> Hands-on: Filtering for read pair mapping status
> 
> 1.  **Samtools view** ( Galaxy version 1.15.1+galaxy0) with the following parameters (leaving non-mentioned ones at their defaults):
>     
>     *   param-files _“SAM/BAM/CRAM data set”_: all 3 mapped reads datasets of the family trio, outputs of **Map with BWA-MEM** tool
>     *   _“What would you like to look at?”_: `A filtered/subsampled selection of reads`
>     
>     *   In _“Configure filters”_:
>         *   _“Exclude reads with any of the following flags set”_: `Read is unmapped` **and** `Mate is unmapped`

This will result in three new datasets, one for each sample in the analysis.

> Details: More than one way to filter
> 
> Instead of the above filter conditions we could also have exploited the _Read is mapped in a proper pair_ flag bit.
> 
> For a read to be flagged as being mapped in a proper pair its mate needs to be mapped, but the mapped pair also needs to meet additional, aligner-specific criteria. These may include (and do so for **BWA-MEM**):
> 
> *   both reads need to map to the same chomosome
> *   the two reads need to map to the reference in the orientation expected by the aligner
> *   the two read pairs need to map to the reference within an aligner-determined distance
> 
> Thus, filtering based on the flag has two consequences:
> 
> *   filtering will be stricter than with just the _Read is unmapped_ and _Mate is unmapped_ flags above
> *   you will eliminate read pairs that could be informative with regard to chromosomal rearrangements and insertion/deletion events
>     
>     \=> Do not filter for properly paired reads if you plan to detect such structural variants!
>     
> 
> In addition, the _proper pair_ flag is considered undefined if the read itself is unmapped, so a _proper pair_ filter should eliminate unmapped reads explicitly to be on the safe side.
> 
> Thus, if you would like to use proper pair filtering (we have no intention to detect structural variants in this tutorial) instead of just filtering for mapped reads with a mapped mate, you could run the alternative:
> 
> > Hands-on
> > 
> > 1.  **Samtools view** ( Galaxy version 1.15.1+galaxy0):
> >     
> >     *   param-files _“SAM/BAM/CRAM data set”_: all 3 mapped reads datasets of the family trio, outputs of **Map with BWA-MEM** tool
> >     *   _“What would you like to look at?”_: `A filtered/subsampled selection of reads`
> >     
> >     *   In _“Configure filters”_:
> >         *   _“Require that these flags are set”_: `Read is mapped in a proper pair`
> >         *   _“Exclude reads with any of the following flags set”_: `Read is unmapped`

## Removing duplicate reads

> Hands-on: Remove duplicates
> 
> 1.  **RmDup** ( Galaxy version 2.0.1) with the following parameters:
>     
>     *   param-files _“BAM file”_: all 3 filtered reads datasets; the outputs of **Samtools view**
>     *   _“Is this paired-end or single end data”_: `BAM is paired-end`
>     
>     *   _“Treat as single-end”_: `No`

Again, this will produce three new datasets, one for each member of the family trio.

# Variant calling

With the sequenced reads of all samples mapped and postprocessed, we can start looking for evidence of sequence deviations, _i.e._ variants, between the sequenced genomic samples and the reference genome.

This task has been automated and optimized continuously over the last decade, and modern variant calling software hides much of the complexity involved in it. At least a basic understanding of the underlying concepts is still highly recommended though and, if you are new to variant calling, the tutorial on [Calling variants in diploid systems](../dip/tutorial.html) may be a good starting point for you.

## Generating FreeBayes calls

We will use **FreeBayes** to call our variants. **FreeBayes** is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.

> Hands-on: Generating FreeBayes calls
> 
> 1.  Run **FreeBayes** ( Galaxy version 1.3.6+galaxy0):
>     *   _“Choose the source for the reference genome”_: `Locally cached`
>         
>         *   _“Run in batch mode?”_: `Merge output VCFs`
>             *   param-files _“BAM dataset(s)”_: all three mapped reads datasets of the family trio; the outputs of **RmDup**
>         *   _“Using reference genome”_: `Human: hg19` (or a similarly named option)
>         
>         > Comment: Using the imported \`hg19\` sequence
>         > 
>         > If you have imported the `hg19` chr8 sequence as a fasta dataset into your history instead:
>         > 
>         > *   _“Choose the source for the reference genome”_: `History`
>         >     *   _“Run in batch mode?”_: `Merge output VCFs`
>         >         *   param-files _“BAM or CRAM dataset(s)”_: all three mapped and fully post-processed reads datasets of the family trio; outputs of **RmDup**
>         >     *   param-file _“Use the following dataset as the reference sequence”_: your imported `hg19` fasta dataset.
>         
>     *   _“Limit variant calling to a set of regions?”_: `Do not limit`
>     *   _“Choose parameter selection level”_: `1. Simple diploid calling`

**Congratulations!** You have created you first multisample VCF file, one of the most complicated file formats in bioinformatics. For every variant detected in at least one of your samples, this tab-separated format uses a single line to store all information about the variant. This includes but is not limited to:

*   the position of the variant in the genome (with respect to the reference genome used for the analysis)
*   the nature of the variant (the actual sequence change associated with it)
*   the detected genotype of every sample at the variant position
*   measures of the reliability of the variant call and of all individual genotype calls

> Hands-on: Optional hands-on: Inspect the VCF output produced by FreeBayes
> 
> 1.  Display the VCF dataset:
>     
>     *   Click the galaxy-eye icon next to the VCF dataset generated by FreeBayes to display its contents.
>         
>         VCF is a tabular plain text format though its information density makes it complicated to understand.
>         
>     
>     > Question
>     > 
>     > Can you locate at least some of the above-listed information in the dataset?
>     > 
>     > Hints:
>     > 
>     > *   Lines starting with `##` are comment lines explaining the content of the file.
>     > *   Diploid genotypes at biallelic sites are encoded using `0/0`, `0/1` and `1/1` to represent homozygous reference, heterozygous and homozygous variant states, respectively.

## Post-processing FreeBayes calls

Before starting to analyze the detected variants, we need to post-process the VCF dataset to fix a few incompatibilities between Freebayes and downstream analysis tools. In particular, we want to:

*   Split multiallelic variant records, _i.e._, records that list more than one alternate allele at a given genomic position, into separate record lines.
    
    This will allow us to annotate each record with information about the impact of one specific variant allele further on.
    
*   Make sure that indels are represented in left-aligned and normalized form because this is how known indels are stored in public annotation databases.
    

A tool that can do this and also ensures that a VCF dataset conforms to standards in some other, less important respects is **bcftools norm**.

> Hands-on: Post-processing FreeBayes calls
> 
> 1.  **bcftools norm** ( Galaxy version 1.15.1+galaxy3) with the following parameters:
>     *   _“VCF/BCF Data”_: the VCF output of **FreeBayes** tool
>     *   _“Choose the source for the reference genome”_: `Use a built-in genome`
>         
>         *   _“Reference genome”_: `Human: hg19` (or a similarly named option)
>         
>         > Comment: Using the imported \`hg19\` sequence
>         > 
>         > If you have imported the `hg19` chr8 sequence as a fasta dataset into your history instead:
>         > 
>         > *   _“Choose the source for the reference genome”_: `Use a genome from the history`
>         >     *   _“Reference genome”_: your imported `hg19` fasta dataset
>         
>     *   _“When any REF allele does not match the reference genome base”_: `ignore the problem (-w)`
>     *   _“Atomize”_: `No`
>     *   _“Left-align and normalize indels?”_: `Yes`
>     *   _“Perform deduplication for the folowing types of variant records”_: `do not deduplicate any records`
>         
>         Freebayes is not producing any duplicate calls.
>         
>     *   _“~multiallelics”_: `split multiallelic sites into biallelic records (-)`
>         *   _“split the following variant types”_: `both`
>             
>             We want to split both, multiallelic SNP and indel records.
>             
>     *   _“output\_type”_: `uncompressed VCF`
>         
>         We would like to keep the results human-readable. VCF is also what tools like _SnpEff_ and _GEMINI_ expect as input. Compressed, binary BCF is interesting for space-efficient long-term storage of large lists of variants.
>         

You could try to look for the differences between the original and the normalized VCF dataset, but for convenience **bcftools norm** reports a brief summary of the actions it performed. Expand the dataset in the history (by clicking on its name) to see this output listing the total number of variant lines processed, along with the number of split, realigned and skipped records.

# Variant annotation and reporting

A list of variants detected in a set of samples is a start, but to discover biologically or clinically relevant information in it is almost impossible without some additional tools and data. In particular, the variants in the list need to be:

*   **prioritized** with respect to their potential relevance for the biological / clinical phenotype that is studied
    
    Even with exome sequencing, only a fraction of the detected variants will have a clear impact on the function of a protein (many variants will introduce silent mutations, or reside in intronic regions still covered by the exome-enriched sequencing data). Of these, many will have been observed before in healthy individuals arguing against them playing an important role in an adverse phenotype.
    
*   **filtered** based on the inheritance pattern expected for a causative variant
    
    A multisample VCF file records the most likely genotypes of all samples at every variant site. Knowing which individuals (samples) are affected by a phenotype we can exclude variants with inheritance patterns that are incompatible with the observed inheritance of the phenotype.
    
*   **reported** in a more human-friendly form
    
    While the VCF format can be used to encode all relevant information about any variant, it is hard for humans to parse that information.
    
    Ideally, one would like to generate simpler reports for any set of filtered and prioritized variants.
    

## Get data

Our workhorse for annotating and reporting variants and the genes affected by them will be the **GEMINI** framework. GEMINI comes bundled with a wealth of annotation data for human variants from many different sources. These can be used to annotate any list of human variants conveniently, without the need for separate downloads and conversion between different annotation data formats.


The only additional annotation tool we need, for the purpose of this tutorial, is the tool **SnpEff**, which can annotate variants with their calculated effects on known genomic features. Because SnpEff is a generic tool that can be used on variants found in the genome of any organism we need to provide it with a so-called **SnpEff genome file** that holds the annotated features (genes, transcripts, translated regions, _etc._) for our genome of interest.

While annotated variants are all we need to _prioritize_ them as described above, _filtering based on inheritance patterns_ requires a way to inform GEMINI about the relationship between our samples and their observed phenotypes. This is done through a so-called **pedigree file** in PED format, which is rather simple to generate manually.

> Hands-on: Obtain SnpEff genome and GEMINI pedigree files
> 
> 1.  Download **SnpEff functional genomic annotations**
>     
>     > Comment: Shortcut
>     > 
>     > You can skip this step if the Galaxy server you are working on offers `Homo sapiens: hg19` as a locally installed snpEff database. You can check the **Genome source** select list of the **SnpEff eff** ( Galaxy version 4.3+T.galaxy2) tool to see if this is the case.
>     
>     Use **SnpEff Download** ( Galaxy version 4.3+T.galaxy2) to download genome annotation database `hg19`.
>     
> 2.  Create a PED-formatted pedigree dataset describing our single-family sample trio:
>     
>     ```
>     #family_id    name     paternal_id    maternal_id    sex    phenotype
>     FAM           father   0              0              1      1
>     FAM           mother   0              0              2      1
>     FAM           proband  father         mother         1      2
>     ```
>     
>     and set its datatype to `tabular`.
>     
>     > Tip: Creating a new file
>     > 
>     > *   Click galaxy-upload **Upload Data** at the top of the tool panel
>     > *   Select galaxy-wf-edit **Paste/Fetch Data** at the bottom
>     > *   Paste the file contents into the text field
>     > *   Change **Type** from “Auto-detect” to `tabular`\* Press **Start** and **Close** the window
>     
>     

## Variant annotation with functional genomic effects

We need to start annotating our variants with SnpEff simply because Gemini knows how to parse SnpEff-annotated VCFs, while GEMINI output cannot be used with SnpEff.

> Hands-on: Adding annotations with SnpEff
> 
> 1.  **SnpEff eff** ( Galaxy version 4.3+T.galaxy2)
>     *   param-file _“Sequence changes (SNPs, MNPs, InDels)”_: the output of **bcftools norm** tool
>     *   _“Input format”_: `VCF`
>     *   _“Output format”_: `VCF (only if input is VCF)`
>     *   _“Genome source”_: `Locally installed reference genome`
>         
>         *   _“Genome”_: `Homo sapiens: hg19` (or a similarly named option)
>         
>         > Comment: Using the imported \`hg19\` SnpEff genome database
>         > 
>         > If you have imported the `hg19` SnpEff genome database into your history instead:
>         > 
>         > *   _“Genome source”_: `Downloaded snpEff database in your history`
>         >     *   param-file _“SnpEff4.3 Genome Data”_: your imported `hg19` SnpEff dataset.
>         
>     *   _“Produce Summary Stats”_: `Yes`

Running the above job will produce two datasets. One is a _Summary Stats_ HTML report, which contains some interesting general metrics such as a distribution of variants across gene features. The other one is the main annotation result - a VCF like the input, but with annotations of variant effects added to the INFO column.

> Hands-on: Optional hands-on: Inspect the Summary Stats output produced by SnpEff
> 
> 1.  Display the dataset:
>     
>     *   Click the galaxy-eye icon next to the HTML dataset generated by SnpEff to display its contents.
>     
>     > Question
>     > 
>     > One section in the report is **Number of effects by type and region**. Given that you are analyzing exome data, what is the most surprising aspect in this section? Do you have an idea how to explain it?
>     >      

## Generating a GEMINI database of variants for further annotation and efficient variant queries

Next, we are going to use the SnpEff-annotated VCF as the basis for more exhaustive annotation with GEMINI. Unlike SnpEff, GEMINI does not just add annotations to a list of variants in VCF format. Instead the framework extracts the variants from the VCF input and stores them, together with newly added annotations, in an SQL database. It then lets you formulate queries for retrieving and reporting subsets of variants. The combined variant extraction/annotation/storage step is performed by the **GEMINI load** tool. In addition, that same tool can be used to incorporate sample pedigree info into the database.

> Hands-on: Creating a GEMINI database from a variants dataset
> 
> 1.  **GEMINI load** ( Galaxy version 0.20.1+galaxy2) with
>     *   param-file _“VCF dataset to be loaded in the GEMINI database”_: the output of **SnpEff eff** tool
>     *   _“The variants in this input are”_: `annotated with snpEff`
>     *   _“This input comes with genotype calls for its samples”_: `Yes`
>         
>         Sample genotypes were called by Freebayes for us.
>         
>     *   _“Choose a gemini annotation source”_: select the latest available annotations snapshot (most likely, there will be only one)
>     *   _“Sample and family information in PED format”_: the pedigree file prepared above
>     *   _“Load the following optional content into the database”_
>         
>         *   param-check _“GERP scores”_
>         *   param-check _“CADD scores”_
>         *   param-check _“Gene tables”_
>         *   param-check _“Sample genotypes”_
>         *   param-check _“variant INFO field”_
>         
>         Leave **unchecked** the following:
>         
>         *   _“Genotype likelihoods (sample PLs)”_
>             
>             Freebayes does not generate these values
>             
>         *   _“only variants that passed all filters”_
>             
>             This setting is irrelevant for our input because Freebayes did not apply any variant filters.
>             

Running this job generates a GEMINI-specific database dataset, which can only be processed with other GEMINI tools. The benefit, however, is that we now have variants, rich annotations and pedigree info stored in a format that enables flexible and highly efficient queries, which will greatly simplify our actual task to identify the variant responsible for the child’s disease!


## Candidate variant detection

Let us now try to identify variants that have the potential to explain the boy child’s osteopetrosis phenotype. Remember that the parents are consanguineous, but both of them do not suffer from the disease. This information allows us to make some assumptions about the inheritance pattern of the causative variant.

> Question
> 
> Which inheritance patterns are in line with the phenotypic observations for the family trio?
> 
> Hint: GEMINI easily lets you search for variants fitting any of the following inheritance patterns:
> 
> *   Autosomal recessive
> *   Autosomal dominant
> *   X-linked recessive
> *   X-linked dominant
> *   Autosomal de-novo
> *   X-linked de-novo
> *   Compound heterozygous
> *   Loss of heterozygosity (LOH) events
> 
> Think about which of these might apply to the causative variant.
> 
> > Solution
> > 
> > *   Since both parents are unaffected the variant cannot be dominant and inherited.
> > *   A de-novo acquisition of a dominant (or an X-linked recessive) mutation is, of course, possible.
> > *   A recessive variant is a possibility, and a more likely one given the parents’ consanguinity.
> > *   For both the de-novo and the inherited recessive case, the variant could reside on an autosome or on the X chromosome. Given that we provided you with only the subset of sequencing reads mapping to chr8, an X-linked variant would not be exactly instructive though ;)
> > *   A compound heterozygous combination of variant alleles affecting the same gene is possible, but less likely given the consanguinity of the parents (as this would require two deleterious variant alleles in the gene circulating in the same family).
> > *   A loss of heterozygosity (LOH) turning a heterozygous recessive variant into a homozygous one could be caused by uniparental disomy or by an LOH event early in embryonic development, but both these possibilities have an exceedingly small probability.
> > 
> > Based on these considerations it makes sense to start looking for inherited autosomal recessive variants first. Then, if there is no convincing candidate mutation among them, you could extend the search to de-novo variants, compund heterozygous variant pairs and LOH events - probably in that order.

Since our GEMINI database holds the variant and genotype calls for the family trio and the relationship between the family members, we can make use of **GEMINI inheritance pattern** tool to report all variants fitting any specific inheritance model with ease.

Below is how you can perform the query for inherited autosomal recessive variants. Feel free to run analogous queries for other types of variants that you think could plausibly be causative for the child’s disease.

> Hands-on: Finding and reporting plausible causative variants
> 
> 1.  **GEMINI inheritance pattern** ( Galaxy version 0.20.1)
>     
>     *   _“GEMINI database”_: the GEMINI database of annotated variants; output of **GEMINI load** tool
>     *   _“Your assumption about the inheritance pattern of the phenotype of interest”_: `Autosomal recessive`
>         *   param-repeat _“Additional constraints on variants”_
>             *   _“Additional constraints expressed in SQL syntax”_: `impact_severity != 'LOW'`
>                 
>                 This is a simple way to prioritize variants based on their functional genomic impact. Variants with _low impact severity_ would be those with no obvious impact on protein function (_i.e._, silent mutations and variants outside coding regions)
>                 
>         *   _“Include hits with less convincing inheritance patterns”_: `No`
>             
>             This option is only meaningful with larger family trees to account for errors in phenotype assessment.
>             
>         *   _“Report candidates shared by unaffected samples”_: `No`
>             
>             This option is only meaningful with larger family trees to account for alleles with partial phenotypic penetrance.
>
>             *   _“Family-wise criteria for variant selection”_: keep default settings
>         
>             This section is not useful when you have data from just one family.
>
>         *   In _“Output - included information”_
> 
>         *   _“Set of columns to include in the variant report table”_: `Custom (report user-specified columns)`
>             *   _“Choose columns to include in the report”_:
>                 *   param-check _“alternative allele frequency (max\_aaf\_all)”_
>             *   _“Additional columns (comma-separated)”_: `chrom, start, ref, alt, impact, gene, clinvar_sig, clinvar_disease_name, clinvar_gene_phenotype, rs_ids`
>     


> Question
> 
> From the GEMINI reports you generated, can you identify the most likely candidate variant responsible for the child’s disease?

# Conclusion

It was not hard to find the most likely causative mutation for the child’s disease (you did find it, right?).

Even though it will not always provide as strong support for just one specific causative variant, analysis of whole-exome sequencing data of family trios (or other related samples) can often narrow down the search for the cause of a genetic disease to just a very small, manageable set of candidate variants, the relevance of which can then be addressed through standard methods.

> Key points
> 
> *   Exome sequencing is an efficient way to identify disease-relevant genetic variants.
>     
> *   Freebayes is a good variant and genotype caller for the joint analysis of multiple samples. It is straightforward to use and requires only minimal processing of mapped reads.
>     
> *   Variant annotation and being able to exploit genotype information across family members is key to identifying candidate disease variants. SnpEff and GEMINI, in particular, are powerful tools offered by Galaxy for that purpose.
>     



###Resources used for this tutorial:

1.  Wolfgang Maier, Bérénice Batut, Torsten Houwaart, Anika Erxleben, Björn Grüning, **Exome sequencing data analysis for diagnosing a genetic disease (Galaxy Training Materials)**. [https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html) Online; accessed TODAY
2.  Hiltemann, Saskia, Rasche, Helena et al., 2023 **Galaxy Training: A Powerful Framework for Teaching!** PLOS Computational Biology [10.1371/journal.pcbi.1010752](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010752)
3.  Batut et al., 2018 **Community-Driven Data Analysis Training for Biology** Cell Systems [10.1016/j.cels.2018.05.012](https://doi.org/10.1016%2Fj.cels.2018.05.012)

