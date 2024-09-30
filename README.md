# 🧬 Lab: Next generation sequencing analysis 1&2 🧬

In this tutorial, we will guide you through the essential steps involved in processing high-quality reads from Next-Generation Sequencing (NGS) and mapping them to a reference genome. You will learn how to perform quality control on short-read NGS data, carry out the mapping of high-quality reads, and understand the post-mapping processes necessary for effective variant calling and annotation. By the end of this tutorial, you will have a comprehensive understanding of the entire workflow, including data preparation, quality control, read mapping, and the generation of variant databases. 

## Learning outcomes
### NGS analysis 1
- Review the structure and origin of sequencing reads.
- Outline the file formats in which reads and mapped reads are stored in.
- Describe the basic quality control steps that are necessary for processing of read data.
- Describe the process of read mapping.
  
### NGS analysis 2
- Understand the motivation to identify variants from NGS data in a clinical setting.
- Describe the GATK best practises for variant identification.
- Describe the file format that variants are stored in.
- Outline the primary measures of variant quality control.

# Before you begin
1. Create a Galaxy user account with your RCSI email: [Sign up here](https://usegalaxy.org/login/start?redirect=None) 
2. Download the "data_files" folder [here](https://rcsicampus-my.sharepoint.com/:f:/g/personal/laurawhelan_rcsi_com/EkI2pyMKZNxOjeDJOtqnB9EB3L5pV0j_TEIivBL5suTB7A?e=x4cKab).
3. Create a new "history" within Galaxy called "ATT_NGS_LAB"

![history_create_new](https://github.com/user-attachments/assets/4d6c3652-22f2-4612-ad62-ae78b5c13c4b)

Note: I have "pre-made" all the files for you. You're going to perform all the steps to make these files too, but computationally some of these steps take a long time. Because of this, we have pre-made files waiting for you. Think of it like a cooking show - here's one we made earlier!

# Let's get started
##There will be a video to follow on screen!

#Data Preparation:
###Steps
1. Upload your data to your ATT_NGS_LAB history
   
2. Check that the newly created datasets in your history have their datatypes assigned correctly to fastqsanger.gz, and fix any missing or wrong datatype assignment

3. Add #father/#mother/#child tags to the ".gz" files 
   
   To tag a dataset:
      - Click on the dataset to expand it
      - Click on Add Tags galaxy-tags
      - Add tag text. Tags starting with # will be automatically propagated to the outputs of tools using this dataset (see below).
      - Press Enter
      - Check that the tag appears below the dataset name

4. Update the type of file and genome build for all ".gz" files as follows:
   
      Type: fasta
   
      Genome: Human Feb. 2009 (GRCh37/hg19) (hg19)

#Quality control

This step serves the purpose of identifying possible issues with the raw sequenced reads input data before embarking on any “real” analysis steps.

###Steps
1. Run  FastQC on each of your six fastq datasets
2. Use  MultiQC to aggregate the raw FastQC data of all input datasets into one report
3. Inspect the Webpage output produced by the tool

#Read mapping

Now that you confirmed that the quality of the input data is good enough to warrant further analysis, it is time to map the sequenced reads to the reference genome.
###Steps
1. Map with BWA-MEM ( Galaxy version 0.7.17.2) to map the reads from the father sample to the reference genome
   
- “Will you select a reference genome from your history or use a built-in index?”: Use a built-in genome index
- “Using reference genome”: Human: hg19 (or a similarly named option)
- “Single or Paired-end reads”: Paired
- “Select first set of reads”: the forward reads (R1) dataset of the father sample
- “Select second set of reads”: the reverse reads (R2) dataset of the father sample
- “Set read groups information?”: Set read groups (SAM/BAM specification)
- “Auto-assign”: No
- “Read group identifier (ID)”: 000
- “Auto-assign”: No
- “Read group sample name (SM)”: father

 2. Map with BWA-MEM ( Galaxy version 0.7.17.2) to map the reads from the mother sample to the reference genome using the same parameters as before except

- “Single or Paired-end reads”: Paired
- “Select first set of reads”: the forward reads (R1) dataset of the mother sample
- “Select second set of reads”: the reverse reads (R2) dataset of the mother sample
- “Set read groups information?”: Set read groups (SAM/BAM specification)
- “Auto-assign”: No
- “Read group identifier (ID)”: 001
- “Auto-assign”: No
- “Read group sample name (SM)”: mother

3. Map with BWA-MEM ( Galaxy version 0.7.17.2) to map the reads from the child sample to the reference genome using the same parameters as before except

- “Single or Paired-end reads”: Paired
- “Select first set of reads”: the forward reads (R1) dataset of the child sample
- “Select second set of reads”: the reverse reads (R2) dataset of the child sample
- “Set read groups information?”: Set read groups (SAM/BAM specification)
- “Auto-assign”: No
- “Read group identifier (ID)”: 002
- “Auto-assign”: No
- “Read group sample name (SM)”: proband
