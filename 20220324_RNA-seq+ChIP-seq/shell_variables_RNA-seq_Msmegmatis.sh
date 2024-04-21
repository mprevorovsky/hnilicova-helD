#!/bin/bash

# CPU threads
CPU=7
# output directory for raw read QC
QC_dir="./QC_RNA-seq_Msmegmatis/"
# output directory for quality-trimmed read QC
QC_trim_dir="./QC_trim_RNA-seq_Msmegmatis/"
# FASTQ file directory
fastq_dir="./FASTQ_RNA-seq_Msmegmatis/"
# quality-trimmed FASTQ file directory
fastq_trim_dir="./FASTQ_trim_RNA-seq_Msmegmatis/"
# FASTQ file extension
fastq_file_ext="\.fastq\.gz$"
# genome sequence and annotation folder
genome_dir="./genome_Msmegmatis/"
# file containing reference genome sequence
genome="${genome_dir}M_smegmatisMC2_155.fasta"
# genome annotation file
genome_annot="${genome_dir}GCF_000015005.1_ASM1500v1_genomic.customized.gff"
# BAM file directory
bam_dir="./BAM_RNA-seq_Msmegmatis/"
# how to perform binning of genome coverage
bin_size=1
# images directory
image_dir="./images_Msmegmatis/"
# output of multiBamSummary
bam_summary_file="multiBamSummary.npz"
# directory for genome coverage data
coverage_dir="./coverage_RNA-seq_Msmegmatis/"
# directory where Trimmomatic is installed
trimmomatic_dir="/opt/Trimmomatic-0.39/"
