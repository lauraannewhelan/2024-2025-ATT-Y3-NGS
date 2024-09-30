# 🧬 Lab: Next-Generation Sequencing Analysis 1 & 2 🧬

In this tutorial, we will guide you through the essential steps involved in processing high-quality reads from Next-Generation Sequencing (NGS) and mapping them to a reference genome. You will learn how to perform quality control on short-read NGS data, map high-quality reads, and understand the post-mapping processes necessary for effective variant calling and annotation. 

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

# Let's Get Started
> **There will be a video to follow on screen!**

## Data Preparation

### Steps
1. **Upload your data** to your "ATT_NGS_LAB" history.
2. **Check dataset types**: Ensure the datasets have their datatypes assigned correctly to `fastqsanger.gz`. Fix any missing or incorrect datatype assignments.
3. **Tag datasets** as #father, #mother, or #child for the ".gz" files.

    - To tag a dataset:
      - Click on the dataset to expand it.
      - Click on "Add Tags".
      - Add the appropriate tag (tags starting with `#` will propagate to tool outputs).
      - Press Enter.
      - Check that the tag appears below the dataset name.

4. **Update file type and genome build** for all ".gz" files:
   - **Type**: fasta
   - **Genome**: Human Feb. 2009 (GRCh37/hg19) (hg19)

---

## Quality Control

This step is to identify any possible issues with the raw sequencing read data before starting the analysis.

### Steps
1. **Run FastQC** on each of your six `.fastq` datasets.
2. **Use MultiQC** to aggregate the raw FastQC data into one comprehensive report.
3. **Inspect the MultiQC output** for potential issues.

---

## Read Mapping

After confirming that the quality of the input data is acceptable, it’s time to map the sequencing reads to the reference genome.

### Steps

1. **Map the father’s reads** using BWA-MEM (Galaxy version 0.7.17.2):
   - **Reference Genome**: Use a built-in genome index (Human: hg19).
   - **Read Type**: Paired-end.
   - **First Set of Reads**: Forward reads (R1) of the father.
   - **Second Set of Reads**: Reverse reads (R2) of the father.
   - **Read Group Information**: Set read groups.
     - **Read Group ID**: 000
     - **Sample Name**: father

2. **Map the mother’s reads** using BWA-MEM with the same parameters, except:
   - **Read Group ID**: 001
   - **Sample Name**: mother

3. **Map the child’s reads** using BWA-MEM with the same parameters, except:
   - **Read Group ID**: 002
   - **Sample Name**: proband

> **Note:** Read mapping is step that can take a considerable amount of time, so for the next steps we'll use pre made bam files for the next steps. 

---

## Mapped Reads Postprocessing

Earlier, you should have uploaded three `.bam` files — one for the father, one for the mother, and one for the proband. These are exactly what would be produced by BWA-MEM.

### Steps

1. **Tag datasets** as `#father`, `#mother`, or `#child` for the `.bam` files:
   
   - Click on the dataset to expand it.
   - Click on "Add Tags".
   - Add the appropriate tag (tags starting with `#` will propagate to tool outputs).
   - Press Enter.
   - Check that the tag appears below the dataset name.

2. **Update database/build** for all `.bam` files:
   
   - Click the desired dataset’s name to expand it.
   - Click on the “?” next to the database indicator:
   
     ![Screenshot](https://github.com/user-attachments/assets/3e10afc7-6148-4433-b393-41f945126ada)
   
   - In the central panel, change the **Database/Build** field.
   - Select your desired database key from the dropdown list: `Human Feb. 2009 (GRCh37/hg19) (hg19)`.
   - Click the **Save** button.

3. **Filtering** on mapped reads properties:

   Run **Samtools view** (Galaxy version 1.15.1+galaxy0) with the following parameters (leave non-mentioned ones at their defaults):

   - **SAM/BAM/CRAM data set**: Select all three mapped reads datasets of the family trio (outputs of the BWA-MEM tool).
   - **What would you like to look at?**: A filtered/subsampled selection of reads
   - **Configure filters**: “Exclude reads with any of the following flags set”: Read is unmapped and Mate is unmapped
