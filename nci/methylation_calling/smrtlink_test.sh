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

reference_fasta=$input/ref_pcontig_019.fasta
bax_files=$short/raw_data/pacbio/pbrun1/SCH1743_lib20160318_0.14_F01/1.bax.h5
preset=$short/myapps/smrtlink/5.1.0/smrtcmds/bin/preset.xml
basemod_preset=$short/myapps/smrtlink/5.1.0/smrtcmds/bin/preset_basemod.xml

OUTPUT=$short/methylation_calling/${name}_mc/smrtlink

threads=16
mem_size='120G'

#OUTPUT
#make the output folder
mkdir -p ${OUTPUT}

module load smrtlink

cd $PBS_JOBFS
subread_dir=${PBS_JOBFS}/subreadset

#Probably don't need this for RSII data
#Index raw data
#samtools index pathToSubreadBams/*.bam 
#pbindex pathToSubreadBams/*.bam

#Create the dataset for RSII bax.h5 files
dataset create --type HdfSubreadSet hdfsubreadset_test.xml $bax_files
pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.sa3_hdfsubread_to_subread --preset-xml preset.xml \
-e eid_hdfsubread:hdfsubreadset_test.xml -o $subread_dir

#make a symbolic link to the suubreadset file
subreadset=${subread_dir}/

#Index reference genome
samtools faidx $reference_fasta

#Create ReferenceSet
dataset create --type ReferenceSet contig_019.referenceset.xml $reference_fasta

#Run the resequencing pipeline
pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.ds_modification_detection \
-e eid_subread:subreadset_test.xml -e eid_ref_dataset:contig_019.referenceset.xml \
--preset-xml $preset --preset-xml $basemod_preset -o $OUTPUT
