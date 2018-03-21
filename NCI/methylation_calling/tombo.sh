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

##### fix this for data pre-processing.
# albacore_output.fastq=$short/methylation_calling/input/pcontig_019_aln.fastq

fast5_dir=$short/methylation_calling/input/pcontig_019_aln_fast5.tar.gz
reference_fasta=$short/methylation_calling/input/ref_pcontig_019.fasta
sample_alt_model= (????)

OUTPUT=$short/methylation_calling/output/test_contig_019/tombo

threads=16
mem_size='120G'


#make the output folder
mkdir -p ${OUTPUT}


#now move everything to the node so we can get started
cd $PBS_JOBFS

#Comparing to an alternative 5mC and 6mA model (recommended method)
tombo test_significance --fast5-basedirs fast5_dir \
    --alternate-bases 5mC 6mA --statistics-file-basename sample_alt_model


# output the genomic sequence surrounding locations with the largest fraction of modified reads
tombo write_most_significant_fasta --statistics-filename sample_alt_model.6mA.tombo.stats \
    --genome-fasta reference_fasta
