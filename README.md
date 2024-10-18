# 🧬 Lab: Next-Generation Sequencing Analysis 1 & 2 🧬

## Learning Outcomes

### NGS Analysis 1
- Review the process of processing high quality reads and mapping to a reference genome.	
- Perform quaity control of short-read NGS data.	 
- Carrying out mapping of high-quality reads to a reference genome.		

### NGS Analysis 2
- Review the file structure which stores discovered variants.
- Outline the quality control parameters that are used to filter detected variants in NGS analysis.
- Perform filtering of NGS genotype data to identify a variant of interest.

---

# 🛠 Before You Begin

1. **Create a Galaxy user account** using your RCSI email: [Sign up here](https://usegalaxy.org/login/start?redirect=None).
2. **Download the "data_files" folder** from [this link](https://rcsicampus-my.sharepoint.com/:f:/g/personal/laurawhelan_rcsi_com/EkI2pyMKZNxOjeDJOtqnB9EB3L5pV0j_TEIivBL5suTB7A?e=x4cKab).

> **Note:** I have "pre-made" all the files for you. You're going to perform all the steps to make these files, but some of these steps take a long time computationally. That's why we have pre-made files ready for you, similar to a cooking show — *here’s one we made earlier!*

---

# Exome Sequencing Data Analysis for Diagnosing a Genetic Disease

> **Overview**  
> 
> - How do you identify genetic variants in samples based on exome sequencing data?
> - How do you, among the set of detected variants, identify candidate causative variants for a given phenotype/disease?

> **Objectives:**
> 
> - Jointly call variants and genotypes for a family trio from whole-exome sequencing data.
> - Use variant annotation and the observed inheritance pattern of a phenotype to identify candidate causative variants and prioritize them.

Exome sequencing is a method that enables the selective sequencing of the exonic regions of a genome — that is, the transcribed parts of the genome present in mature mRNA, including protein-coding sequences and untranslated regions (UTRs).

In humans, there are about 180,000 exons with a combined length of ~30 million base pairs (30 Mb). 

Exome sequencing offers an affordable alternative to whole-genome sequencing in the diagnosis of genetic diseases. 

> **Details: Exome Sequencing _vs_ Whole-Genome Sequencing**  
> 
> The steps in this tutorial are suitable for analyzing whole-genome sequencing (WGS) data. However, at comparable coverage, WGS datasets are much larger than exome datasets, and their analysis takes more time.  
> 
> While WGS allows for variant detection in more regions of the genome, it also enables detection of _copy number variation (CNV)_ and _structural variants_ like translocations and inversions (requiring more sophisticated analyses, not covered here).  
> 
> Generally, exome sequencing captures most information analyzable with standard bioinformatic tools at reasonable costs. WGS, however, captures the maximum information that current sequencing technology can provide, with the potential for future reanalysis using advanced tools.

Identifying causative variants for a genetic disease requires variant comparison between the patient and selected relatives. Family trio data, consisting of the genome sequences of the patient and their parents, is often used to detect variants following Mendelian inheritance, _de novo_ mutations, or _loss-of-heterozygosity_ (LOH) events.

> **Agenda**  
> 
> Across these 2 labs, we will cover:
> 
> 1. [Data Preparation](#data-preparation)
>    - [Get data](#get-data)
> 2. [Quality Control](#quality-control)
> 3. [Read Mapping](#read-mapping)
> 4. [Mapped Reads Postprocessing](#mapped-reads-postprocessing)
>    - [Filtering on Mapped Reads Properties](#filtering-on-mapped-reads-properties)
>    - [Removing Duplicate Reads](#removing-duplicate-reads)
> 5. [Variant Calling](#variant-calling)
>    - [Generating FreeBayes Calls](#generating-freebayes-calls)
>    - [Post-processing FreeBayes Calls](#post-processing-freebayes-calls)
> 6. [Variant Annotation and Reporting](#variant-annotation-and-reporting)
>    - [Get Data](#get-data-1)
>    - [Variant Annotation with Functional Genomic Effects](#variant-annotation-with-functional-genomic-effects)
>    - [Generating a GEMINI Database](#generating-a-gemini-database)
>    - [Candidate Variant Detection](#candidate-variant-detection)
> 7. [Conclusion](#conclusion)

---

# 🧬 Lab: Next-Generation Sequencing Analysis 1 🧬

## Let's Get Started

> **There will be a video to follow on screen!**

---

# Data Preparation

In this tutorial, we will analyze exome sequencing data from a family trio, where the boy is affected by [osteopetrosis](https://ghr.nlm.nih.gov/condition/osteopetrosis/), and both consanguineous parents are unaffected. Our goal is to identify the genetic variation responsible for the disease.

## Get Data

### Hands-on: Data Upload

1. **Create a new history** for this tutorial and called ATT_NGS_ANALYSIS
   > **Tip: Creating a New History**
   > - Click the new-history icon at the top of the history panel
   > - **Tip: Renaming a History**  
   >   1. Click the galaxy-pencil (**Edit**) next to the history name (default: “Unnamed history”)  
   >   2. Type the new name  
   >   3. Click **Save**  
   >   4. To cancel renaming, click the “Cancel” button.
   >   5. Upload the data you downloaded earlier (fastq.gz)


2. **Check** that the newly created datasets have their datatypes correctly assigned to `fastqsanger.gz`. Fix any missing or incorrect datatype assignment.

   > **Tip: Changing the Datatype**  
   > - Click the galaxy-pencil **pencil icon** for the dataset.  
   > - In the central panel, click the **Datatypes** tab.  
   > - Select `fastqsanger.gz` from the dropdown list and click **Save**.


3. **Add tags** (#father, #mother, #child) to the datasets.

   > **Tip: Adding a Tag**  
   > - Click on the dataset to expand it.  
   > - Click **Add Tags** and enter the desired text (e.g., `#father`).  
   > - Press Enter and verify the tag appears below the dataset name.


**Congratulations!** You are all set for starting the analysis now.


# Quality Control

This step serves the purpose of identifying possible issues with the raw sequenced reads before embarking on any “real” analysis steps.

Some typical problems with NGS data can be mitigated by preprocessing affected sequencing reads before mapping them to the reference genome. Detecting other, more severe problems early on can at least save you a lot of time analyzing low-quality data.

> **Hands-on: Quality control of the input datasets**
> 
> 1.  Run **FastQC** (Galaxy version 0.74+galaxy0) on each of your six fastq datasets:
>     * _“Short read data from your current history”_: Select all 6 FASTQ datasets with **Multiple datasets**.
>     
>     > **Tip: Select multiple datasets**
>     > 1. Click on **Multiple datasets**.
>     > 2. Select multiple files by holding the Ctrl (or COMMAND) key and clicking on the files.
>     
>     This will add twelve new datasets (one with raw data, another with an HTML report for each input dataset) to your history.
> 
> 2. Use **MultiQC** (Galaxy version 1.11+galaxy1) to aggregate the raw **FastQC** data of all input datasets into one report:
>     * In _“Results”_
>         * _“Which tool was used to generate logs?”_: `FastQC`
>         * In _“FastQC output”_
>             * _“Type of FastQC output?”_: `Raw data`
>             * _“FastQC output”_: Select all six _RawData_ outputs from **FastQC**.
> 
> 3. Inspect the _Webpage_ output produced by the tool.
> 
>     > **Questions:**
>     > 1. Based on the report, will preprocessing of the reads (trimming/filtering) be necessary before mapping?
>     > 2. Why do all samples show a non-normal GC content distribution, and should you be worried?

---

# Read Mapping

Now that you’ve confirmed that the quality of the input data is good enough for further analysis, it’s time to map the reads to the reference genome.

> **Hands-on: Read Mapping**
> 
> 1. **Map with BWA-MEM** (Galaxy version 0.7.17.2) to map the reads from the **father** sample to the reference genome:
>     * _“Will you select a reference genome from your history or use a built-in index?”_: `Use a built-in genome index`.
>         * _“Using reference genome”_: `Human: hg19` (or similar option).
>         
>         > **Comment:** Using the imported `hg19` sequence
>         > If you’ve imported the `hg19` chr8 sequence as a fasta dataset into your history:
>         > * _“Will you select a reference genome from your history or use a built-in index?”_: `Use a genome from history and build index`.
>         > * _“Use the following dataset as the reference sequence”_: Your imported `hg19` fasta dataset.
>     
>     * _“Single or Paired-end reads”_: `Paired`.
>         * _“Select first set of reads”_: Forward reads (R1) of the **father** sample.
>         * _“Select second set of reads”_: Reverse reads (R2) of the **father** sample.
>         
>         > **Tip: No FASTQ datasets selectable?**
>         > Ensure the dataset format is `fastqsanger.gz`.
>         
>     * _“Set read groups information?”_: `Set read groups (SAM/BAM specification)`.
>         * _“Auto-assign”_: `No`.
>             * _“Read group identifier (ID)”_: `000`.
>             * _“Read group sample name (SM)”_: `father`.
>     
> 2. Repeat the process for the **mother** and **child** samples, adjusting the following parameters:
>     - **Mother** sample:  
>         * _“Read group identifier (ID)”_: `001`.  
>         * _“Read group sample name (SM)”_: `mother`.
>     - **Child** sample:  
>         * _“Read group identifier (ID)”_: `002`.  
>         * _“Read group sample name (SM)”_: `proband`.

---

# Mapped Reads Postprocessing

Mapping takes a considerable amount of time, so we'll proceed with the 3 "pre-made" mapped reads datasets in `bam` format, each of which:

- Has its _database_ set to the key `hg19`.
    > **Tip: Changing database/build (dbkey)**
    > 1. Click the dataset’s name to expand it.
    > 2. Click on the “?” next to the database indicator:  
    > 3. In the central panel, change the **Database/Build** field to `Human Feb. 2009 (GRCh37/hg19) (hg19)`.
    > 4. Click **Save**.

- Ideally, carries the tags `#father`, `#mother`, or `#child` for easy identification.

You could use these datasets for variant calling directly, but it’s best practice to perform postprocessing steps on mapped reads before variant calling. In this tutorial, we’ll:

- Filter the paired-end reads to retain only those where both forward and reverse reads are mapped.
- Deduplicate reads to avoid PCR overamplification errors.

## Filtering on Mapped Reads Properties

To produce filtered BAM datasets with only mapped reads where both mates are mapped:

> **Hands-on: Filtering for Read Pair Mapping Status**
> 
> 1. **Samtools view** (Galaxy version 1.15.1+galaxy0) with these parameters:
>     * _“SAM/BAM/CRAM data set”_: All 3 mapped reads datasets.
>     * _“What would you like to look at?”_: `A filtered/subsampled selection of reads`.
>     * _“Exclude reads with any of the following flags set”_: `Read is unmapped` **and** `Mate is unmapped`.

This will result in three new datasets, one for each sample.

## Removing Duplicate Reads

> **Hands-on: Remove Duplicates**
> 
> 1. **RmDup** (Galaxy version 2.0.1) with the following parameters:
>     * _“BAM file”_: All 3 filtered reads datasets.
>     * _“Is this paired-end or single-end data”_: `BAM is paired-end`.
>     * _“Treat as single-end”_: `No`.

This will generate three new datasets for the family trio.


# Variant Calling

With the sequenced reads of all samples mapped and post-processed, we can start looking for evidence of sequence deviations, _i.e._, variants, between the sequenced genomic samples and the reference genome.


## Generating FreeBayes Calls

We will use **FreeBayes** to call variants. **FreeBayes** is a Bayesian genetic variant detector designed to find small polymorphisms, including SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.

> **Hands-on: Generating FreeBayes calls**
> 
> 1. Run **FreeBayes** (Galaxy version 1.3.6+galaxy0):
>     * _“Choose the source for the reference genome”_: `Locally cached`
>         * _“Run in batch mode?”_: `Merge output VCFs`
>             * _“BAM dataset(s)”_: All three mapped reads datasets of the family trio (outputs of **RmDup**).
>         * _“Using reference genome”_: `Human: hg19` (or a similarly named option).
>         
>         > **Comment:** Using the imported `hg19` sequence
>         > If you have imported the `hg19` chr8 sequence as a fasta dataset:
>         > * _“Choose the source for the reference genome”_: `History`
>         > * _“Run in batch mode?”_: `Merge output VCFs`
>         >     * _“BAM or CRAM dataset(s)”_: All three mapped reads datasets of the family trio (outputs of **RmDup**).
>         >     * _“Use the following dataset as the reference sequence”_: Your imported `hg19` fasta dataset.
>     
>     * _“Limit variant calling to a set of regions?”_: `Do not limit`
>     * _“Choose parameter selection level”_: `1. Simple diploid calling`

**Congratulations!** You’ve created your first multisample VCF file, which stores detailed information about detected variants such as:

- The position of the variant in the genome (relative to the reference genome)
- The nature of the variant (the actual sequence change)
- The detected genotype of every sample at the variant position
- Measures of the reliability of the variant and genotype calls

> **Hands-on: Optional - Inspect the VCF output produced by FreeBayes**
> 
> 1. Display the VCF dataset:
>     * Click the galaxy-eye icon next to the VCF dataset generated by FreeBayes.
>     
>     > **Question:**
>     > Can you locate some of the information listed above in the VCF file?
>     > 
>     > * Lines starting with `##` are comment lines explaining the content of the file.
>     > * Diploid genotypes at biallelic sites are encoded as `0/0`, `0/1`, and `1/1` representing homozygous reference, heterozygous, and homozygous variant states.

---

## Post-Processing FreeBayes Calls

Before analyzing the detected variants, we need to post-process the VCF file to resolve certain incompatibilities between FreeBayes and downstream analysis tools. Specifically, we will:

- Split multiallelic variant records into separate lines for each alternate allele.
- Ensure indels are left-aligned and normalized, matching the format used by public annotation databases.

We will use **bcftools norm** to perform these steps.

> **Hands-on: Post-processing FreeBayes calls**
> 
> 1. **bcftools norm** (Galaxy version 1.15.1+galaxy3) with the following parameters:
>     * _“VCF/BCF Data”_: The VCF output from **FreeBayes**.
>     * _“Choose the source for the reference genome”_: `Use a built-in genome`
>         * _“Reference genome”_: `Human: hg19` (or a similarly named option).
>         
>         > **Comment:** Using the imported `hg19` sequence
>         > If you have imported the `hg19` chr8 sequence:
>         > * _“Choose the source for the reference genome”_: `Use a genome from the history`
>         >     * _“Reference genome”_: Your imported `hg19` fasta dataset.
>     
>     * _“When any REF allele does not match the reference genome base”_: `ignore the problem (-w)`
>     * _“Left-align and normalize indels?”_: `Yes`
>     * _“Split multiallelic sites into biallelic records (-)”_: `Yes`
>         * _“Split the following variant types”_: `both`
>     * _“Output type”_: `uncompressed VCF`

Running this job will generate a normalized VCF file. You can expand the dataset in your history to view a summary of the actions performed by **bcftools norm**, such as the number of split, realigned, and skipped records.

---

# Variant Annotation and Reporting

Detecting variants is just the beginning. To discover biologically or clinically relevant information, we need to:

- **Prioritize** variants based on their relevance to the phenotype of interest.
- **Filter** variants based on inheritance patterns, especially when working with multisample data.
- **Report** variants in a more human-readable format than VCF.

We will use the **GEMINI** framework for variant annotation and reporting. Additionally, we need **SnpEff** to annotate variants with their functional genomic effects.

## Get Data

The **SnpEff** tool will annotate the functional effects of variants, while **GEMINI** will handle further annotation and reporting.

> **Hands-on: Obtain SnpEff genome and GEMINI pedigree files**
> 
> 1. **Download SnpEff functional genomic annotations**
>     > **Comment:** Shortcut
>     > If your Galaxy server has `Homo sapiens: hg19` as a locally installed SnpEff database, you can skip this step. Check the **Genome source** list in the **SnpEff eff** tool.
>     
>     Use **SnpEff Download** (Galaxy version 4.3+T.galaxy2) to download the genome annotation database `hg19`.
>     
> 2. Create a PED-formatted pedigree dataset for the family trio:
>     
>     ```
>     #family_id    name     paternal_id    maternal_id    sex    phenotype
>     FAM           father   0              0              1      1
>     FAM           mother   0              0              2      1
>     FAM           proband  father         mother         1      2
>     ```
>     
>     Set its datatype to `tabular`.
>     
>     > **Tip: Creating a new file**
>     > * Click galaxy-upload **Upload Data**.
>     > * Select galaxy-wf-edit **Paste/Fetch Data** and paste the file contents.
>     > * Set **Type** to `tabular` and press **Start**.

---

## Variant Annotation with Functional Genomic Effects

We will start by annotating the variants with **SnpEff**, which will add functional information about the impact of the variants.

> **Hands-on: Adding annotations with SnpEff**
> 
> 1. **SnpEff eff** (Galaxy version 4.3+T.galaxy2):
>     * _“Sequence changes (SNPs, MNPs, InDels)”_: The output from **bcftools norm**.
>     * _“Input format”_: `VCF`
>     * _“Output format”_: `VCF`
>     * _“Genome source”_: `Locally installed reference genome`
>         * _“Genome”_: `Homo sapiens: hg19` (or similar).
>         
>         > **Comment:** Using the imported `hg19` SnpEff genome database
>         > If you have imported the `hg19` SnpEff genome into your history:
>         > * _“Genome source”_: `Downloaded snpEff database in your history`
>         >     * _“SnpEff4.3 Genome Data”_: Your imported `hg19` SnpEff dataset.
>     
>     * _“Produce Summary Stats”_: `Yes`

The result will include a _Summary Stats_ HTML report and the annotated VCF file.

> **Hands-on: Optional - Inspect the Summary Stats output from SnpEff**
> 
> 1. Display the dataset:
>     * Click the galaxy-eye icon next to the HTML dataset generated by SnpEff.
>     
>     > **Question:**
>     > In the **Number of effects by type and region** section, what is surprising, given that you are analyzing exome data?

---

## Generating a GEMINI Database for Further Annotation

Next, we will use **GEMINI** to annotate the variants further and store them in a queryable SQL database.

> **Hands-on: Creating a GEMINI database from a variants dataset**
> 
> 1. **GEMINI load** (Galaxy version 0.20.1+galaxy2):
>     * _“VCF dataset to be loaded in the GEMINI database”_: The output from **SnpEff eff**.
>     * _“The variants in this input are”_: `annotated with snpEff`
>     * _“This input comes with genotype calls for its samples”_: `Yes`
>     * _“Sample and family information in PED format”_: The pedigree file created earlier.
>     * _“Load the following optional content into the database”_: Check `GERP scores`, `CADD scores`, `Gene tables`, `Sample genotypes`, and `variant INFO field`.

Running this job will generate a GEMINI-specific database.

---

## Candidate Variant Detection

Let’s now search for variants that could explain the boy’s osteopetrosis phenotype. The parents are consanguineous, but unaffected, so we can make some assumptions about the inheritance pattern.

> **Question:**
> 
> Which inheritance patterns fit the family trio's phenotypic observations?  
> _Hint: GEMINI allows you to search for patterns such as autosomal recessive, de novo, compound heterozygous, and loss of heterozygosity (LOH)._

> **Solution:**
> 
> - The variant cannot be dominant and inherited since both parents are unaffected.
> - A de novo dominant or X-linked recessive mutation is possible.
> - The most likely scenario, given the consanguinity, is a recessive variant.

We will start by looking for inherited autosomal recessive variants.

> **Hands-on: Finding and reporting plausible causative variants**
> 
> 1. **GEMINI inheritance pattern** (Galaxy version 0.20.1):
>     * _“GEMINI database”_: The GEMINI database of annotated variants (output from **GEMINI load**).
>     * _“Your assumption about the inheritance pattern”_: `Autosomal recessive`.
>     * _“Additional constraints expressed in SQL syntax”_: `impact_severity != 'LOW'` (to prioritize functionally impactful variants).
>     * _“Include hits with less convincing inheritance patterns”_: `No`.
>     * _“Report candidates shared by unaffected samples”_: `No`.
>     * _“Output - included information”_:
>         * _“Set of columns to include”_: `Custom (report user-specified columns)`.
>             * Add `alternative allele frequency (max_aaf_all)` and include: `chrom, start, ref, alt, impact, gene, clinvar_sig, clinvar_disease_name, clinvar_gene_phenotype, rs_ids`.


# Conclusion

It wasn’t difficult to find the most likely causative mutation for the child’s disease (you did find it, right?).

While whole-exome sequencing of family trios may not always point to just one causative variant, it can often narrow down the search to a small, manageable set of candidate variants. These can then be further investigated through standard methods.

> **Key points:**
> 
> - **Exome sequencing** is an efficient method for identifying disease-relevant genetic variants.
> 
> - **FreeBayes** is a reliable variant and genotype caller for the joint analysis of multiple samples. It is easy to use and requires minimal processing of mapped reads.
> 
> - **Variant annotation** and the ability to leverage genotype information across family members are key to identifying candidate disease variants. **SnpEff** and **GEMINI** are particularly powerful tools available in Galaxy for this purpose.

---

### Resources used for this tutorial:

1. Wolfgang Maier, Bérénice Batut, Torsten Houwaart, Anika Erxleben, Björn Grüning, **Exome sequencing data analysis for diagnosing a genetic disease (Galaxy Training Materials)**. [https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html) Online; accessed TODAY.

2. Hiltemann, Saskia, Rasche, Helena et al., 2023 **Galaxy Training: A Powerful Framework for Teaching!** PLOS Computational Biology. [10.1371/journal.pcbi.1010752](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010752).

3. Batut et al., 2018 **Community-Driven Data Analysis Training for Biology** Cell Systems. [10.1016/j.cels.2018.05.012](https://doi.org/10.1016%2Fj.cels.2018.05.012).

### Questions? Contact laurawhelan@rcsi.ie

