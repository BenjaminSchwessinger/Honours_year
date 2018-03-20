#!/bin/bash
#PBS -P sd34
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=120GB
#PBS -l ncpus=16
#PBS -l jobfs=200GB

set -vx

##this is a script to do methylation-calling on contig_19 of reference genome only.

#define some variables at the start
#give this job a name
name=test_contig_019
short=/short/sd34/ap5514
###

#INPUTS
albacore_output.fastq=$short/methylation_calling/input/pcontig_019_aln.fastq
fast5_files=$short/methylation_calling/input/pcontig_019_aln_fast5.tar.gz
reference.fasta=$short/methylation_calling/input/ref_pcontig_019.fasta

OUTPUT=$short/methylation_calling/output/test_contig_019/nanopolish

threads=16
mem_size='120G'


#make the output folder
mkdir -p ${OUTPUT}


#now move everything to the node so we can get started
cd $PBS_JOBFS

# make an index file linking fastq to fast5
nanopolish index -d fast5_files/ albacore_output.fastq
#We get the following files: albacore_output.fastq.index, albacore_output.fastq.index.fai, albacore_output.fastq.index.gzi, and albacore_output.fastq.index.readdb 

# map reads to reference assembly and save as bam file

######## add minimap filepath, add BAM files
minimap2 -a -x map-ont reference.fasta albacore_output.fastq | samtools sort -T tmp -o albacore_output.sorted.bam
samtools index albacore_output.sorted.bam

# call methylation and save to tsv file with log likelihood ratio
nanopolish call-methylation -t 8 -r albacore_output.fastq -b albacore_output.sorted.bam -g reference.fasta -w "chr20:5,000,000-10,000,000" > methylation_calls.tsv

# helper script to make tsv file showing how often each reference position was methylated
/short/sd34/ap5514/myapps/nanopolish/0.9.0/bin/nanopolish/scripts/calculate_methylation_frequency.py -i methylation_calls.tsv > methylation_frequency.tsv


