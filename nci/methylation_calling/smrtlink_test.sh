#!/bin/bash
#PBS -P sd34
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=120GB
#PBS -l ncpus=16
#PBS -l jobfs=200GB

set -vx

##this is a script to do methylation-calling on PacBio sequence data on contig_019 of the Pst-104E genome.


#define some variables at the start
#give this job a name
name=contig_019
short=/short/sd34/ap5514

#INPUTS
input=$short/methylation_calling/input/${name}

albacore_output_fastq=$input/pcontig_019_aln.fastq
fast5_dir=$input/pcontig_019_aln_fast5.tar.gz
reference_fasta=$input/ref_pcontig_019.fasta
index=$short/methylation_calling/contig_019_mc/nanopolish/index/*

PBS_JOBFS_fastq=$PBS_JOBFS/fastq/pcontig_019_aln.fastq
PBS_JOBFS_fast5=$PBS_JOBFS/fast5/pcontig_019_aln_fast5.tar.gz
PBS_JOBFS_ref=$PBS_JOBFS/ref/ref_pcontig_019.fasta

OUTPUT=$short/methylation_calling/${name}_mc/smrtlink

threads=16
mem_size='120G'

#OUTPUT
#make the output folder
mkdir -p ${OUTPUT}


module load smrtlink

cd $PBS_JOBFS

#Index raw data
samtools index pathToSubreadBams/*.bam 
pbindex pathToSubreadBams/*.bam

#Create the dataset
dataset create --type SubreadSet subreadset.xml pathToSubreadBams/*.bam 

#Index reference genome
samtools faidx genome.fasta

#Create ReferenceSet
dataset create --type ReferenceSet genome.referenceset.xml genome.fasta

#Run the resequencing pipeline
pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.ds_modification_detection \
-e eid_subread:subreadset.xml -e eid_ref_dataset:genome.referenceset.xml \
--preset-xml preset.xml --preset-xml preset_basemod.xml -o $OUTPUT


