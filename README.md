# 🧬 Lab: Next-Generation Sequencing Analysis 1 & 2 🧬

In this lab, we will guide you through the essential steps involved in processing high-quality reads from Next-Generation Sequencing (NGS) and mapping them to a reference genome. You will learn how to perform quality control on short-read NGS data, map high-quality reads, and understand the post-mapping processes necessary for effective variant calling and annotation. 

By the end of this tutorial, you will have a comprehensive understanding of the entire workflow, including:
- 🗃️ Data preparation
- ✅ Quality control
- 🧬 Read mapping
- 🧾 Generation of variant databases

---

## 🎯 Learning Outcomes

### 🔬 NGS Analysis 1
- Review the structure and origin of sequencing reads.
- Outline the file formats in which reads and mapped reads are stored.
- Describe the basic quality control steps necessary for processing read data.
- Explain the process of read mapping.

### 🧬 NGS Analysis 2
- Understand the motivation for identifying variants from NGS data in a clinical setting.
- Describe the GATK best practices for variant identification.
- Explain the file formats in which variants are stored.
- Outline the primary measures of variant quality control.

---

# 🛠 Before You Begin

1. **Create a Galaxy user account** using your RCSI email: [Sign up here](https://usegalaxy.org/login/start?redirect=None).
2. **Download the "data_files" folder** from [this link](https://rcsicampus-my.sharepoint.com/:f:/g/personal/laurawhelan_rcsi_com/EkI2pyMKZNxOjeDJOtqnB9EB3L5pV0j_TEIivBL5suTB7A?e=x4cKab).
3. **Create a new "history"** within Galaxy named "ATT_NGS_LAB".

![Create New History](https://github.com/user-attachments/assets/4d6c3652-22f2-4612-ad62-ae78b5c13c4b)

> **Note:** I have "pre-made" all the files for you. You're going to perform all the steps to make these files, but some of these steps take a long time computationally. That’s why we have pre-made files ready for you, similar to a cooking show — *here’s one we made earlier!*

---

# 🧬 Lab: Next-Generation Sequencing Analysis 1 🧬

## 🚀 Let's Get Started
> **There will be a video to follow on screen!**

### 🗃️ Data Preparation

#### Steps
1. **Upload your data** to your "ATT_NGS_LAB" history.
2. **Check dataset types**: Ensure the datasets have their datatypes assigned correctly to `fastqsanger.gz`. Fix any missing or incorrect datatype assignments.
3. **Tag datasets** as #father, #mother, or #child for the `.gz` files.

    - To tag a dataset:
      - Click on the dataset to expand it.
      - Click on "Add Tags".
      - Add the appropriate tag (tags starting with `#` will propagate to tool outputs).
      - Press Enter.
      - Check that the tag appears below the dataset name.

4. **Update file type and genome build** for all `.gz` files:
   - **Type**: `fasta`
   - **Genome**: Human Feb. 2009 (GRCh37/hg19) (hg19)

---

### 🔬 Quality Control

This step is to identify any possible issues with the raw sequencing read data before starting the analysis.

#### Steps
1. **Run FastQC** on each of your six `.fastq` datasets.
2. **Use MultiQC** to aggregate the raw FastQC data into one comprehensive report.
3. **Inspect the MultiQC output** for potential issues.

---

### 🧬 Read Mapping

After confirming that the quality of the input data is acceptable, it’s time to map the sequencing reads to the reference genome.

#### Steps

1. **Map the father’s reads** using BWA-MEM (Galaxy version 0.7.17.2):
   - **Reference Genome**: Use a built-in genome index (Human: hg19).
   - **Read Type**: Paired-end.
   - **First Set of Reads**: Forward reads (R1) of the father.
   - **Second Set of Reads**: Reverse reads (R2) of the father.
   - **Read Group Information**: Set read groups.
     - **Read Group ID**: `000`
     - **Sample Name**: `father`

2. **Map the mother’s reads** using BWA-MEM with the same parameters, except:
   - **Read Group ID**: `001`
   - **Sample Name**: `mother`

3. **Map the child’s reads** using BWA-MEM with the same parameters, except:
   - **Read Group ID**: `002`
   - **Sample Name**: `proband`

> **Note:** Read mapping is a step that can take a considerable amount of time, so for the next steps, we'll use pre-made `.bam` files.

---

### 🧾 Mapped Reads Postprocessing

Earlier, you should have uploaded three `.bam` files — one for the father, one for the mother, and one for the proband. These are exactly what would be produced by BWA-MEM.

#### Steps

1. **Tag datasets** as `#father`, `#mother`, or `#child` for the `.bam` files.

2. **Update database/build** for all `.bam` files:
   - Change the **Database/Build** field to `Human Feb. 2009 (GRCh37/hg19) (hg19)`.

3. **Filtering on mapped reads properties**:
   Run **Samtools view** to exclude unmapped reads from all three mapped reads datasets.

4. **Remove duplicate reads**:
   Run **RmDup** on the filtered reads datasets. Ensure paired-end BAM data is handled appropriately.

---

### 🔬 Variant Calling

Now comes the fun part! We're actually going to find variants in our patient vs the human reference genome — sort of like spot the difference.

#### Steps

Run **FreeBayes** to generate variant calls with these parameters:
- **Reference Genome**: Human: hg19
- **BAM datasets**: The filtered `.bam` datasets for the family trio.
- **Run in batch mode?**: Yes, merge output VCFs.
  
You have now created your first multi-sample VCF file! 🎉

---

# ✅ That's the end of our first NGS lab!

---

# 🧬 Lab: Next-Generation Sequencing Analysis 2 🧬

At the end of our last lab, we created a multi-sample VCF file. Now it's time to do some post-processing and annotation to make it easier to interpret. 

---

### 🧬 Post-processing FreeBayes Calls

#### Candidate Variant Detection

Let’s identify variants that might explain the boy child's osteopetrosis phenotype. Since the parents are consanguineous and unaffected, we can hypothesize the inheritance pattern of the causative variant.

---

### ❓ Question

Which inheritance patterns are in line with the family trio's phenotypic observations?

Hint: GEMINI allows searches for variants fitting the following inheritance patterns:
- Autosomal recessive
- Autosomal dominant
- X-linked recessive
- X-linked dominant
- Autosomal de-novo
- X-linked de-novo
- Compound heterozygous
- Loss of heterozygosity (LOH) events

---

### 💡 Solution

We can use the **GEMINI inheritance pattern tool** to easily report all variants fitting specific inheritance models.

---

### 🧠 Hands-on: Finding and Reporting Plausible Causative Variants

Run **GEMINI inheritance pattern** (Galaxy version 0.20.1) with the following settings:

- **GEMINI database**: Your GEMINI database of annotated variants.
- **Inheritance pattern**: **Autosomal recessive**.
- **Additional constraints**: `impact_severity != 'LOW'` to exclude variants with low impact.

In **Output - included information**:
- **Columns to include**: `chrom, start, ref, alt, impact, gene, clinvar_sig, clinvar_disease_name, clinvar_gene_phenotype, rs_ids`

---

### ❓ Question

Can you identify the most likely candidate variant responsible for the child’s disease?

---
