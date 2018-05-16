#!/bin/bash
#PBS -P sd34
#PBS -q normal
#PBS -l walltime=3:00:00
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
bax_files=$short/raw_data/pacbio/contig_019_bax/*
preset=$short/myapps/smrtlink/5.1.0/smrtcmds/bin/preset.xml
basemod_preset=$short/myapps/smrtlink/5.1.0/smrtcmds/bin/preset_basemod.xml

OUTPUT=$short/methylation_calling/${name}_mc/smrtlink/pbsmrtpipe

threads=16
mem_size='420G'

#OUTPUT
#make the output folder
mkdir -p ${OUTPUT}

module load smrtlink

cd $PBS_JOBFS
subread_dir=${PBS_JOBFS}/subreadset
mkdir $subread_dir

#Create the dataset for RSII bax.h5 files
time dataset create --type HdfSubreadSet hdfsubreadset.xml $bax_files
time pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.sa3_hdfsubread_to_subread --preset-xml $preset \
-e eid_hdfsubread:hdfsubreadset.xml -o $subread_dir

#make a symbolic link to the suubreadset and referenceset files
subreadset=${subread_dir}/tasks/pbcoretools.tasks.gather_subreadset-1/file.subreadset.xml
referenceset=${PBS_JOBFS}/contig_019.referenceset.xml

#Index reference genome
cp $reference_fasta $PBS_JOBFS
ref=${PBS_JOBFS}/ref_pcontig_019.fasta
samtools faidx $ref

#Create ReferenceSet
dataset create --type ReferenceSet $referenceset $ref

#Run the resequencing pipeline
time pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.ds_modification_detection \
-e eid_subread:$subreadset -e eid_ref_dataset:$referenceset \
--preset-xml $preset --preset-xml $basemod_preset -o $OUTPUT

#copy miscellaneous (hdfsubseatset, subreadset, referenceset) back to nci
cp $referenceset $OUTPUT
cp hdfsubreadset.xml $OUTPUT
