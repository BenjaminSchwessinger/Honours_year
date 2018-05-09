#!/bin/bash
#PBS -P sd34
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=120GB
#PBS -l ncpus=16
#PBS -l jobfs=200GB

set -vx

##this is a script to do methylation-calling on contig_019 of the Pst-104E genome.

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

OUTPUT=$short/methylation_calling/${name}_mc/mcaller

threads=16
mem_size='120G'

#OUTPUT
#make the output folder
mkdir -p ${OUTPUT}


#now move everything to the node so we can get started
cd $PBS_JOBFS
mkdir fastq fast5 ref 

cp $albacore_output_fastq $PBS_JOBFS_fastq
cp $fast5_dir $PBS_JOBFS_fast5
cp $reference_fasta $PBS_JOBFS_ref
cp -r $index $PBS_JOBFS/fastq/. #index files from nanopolish were copied over, so first step of mCaller is not needed

#

#go ahead with unzipping
cd $PBS_JOBFS/fast5
for x in *.tar.gz
do
len=${#x}
folder=${x::len-7}
mkdir ${folder}
mv ${x} ${folder}/.
cd ${folder}
tar -xopf ${x}&
cd $PBS_JOBFS/fast5
done
wait
rm */*.tar.gz

#rename this value for the methylation-calling
PBS_JOBFS_fast5=$PBS_JOBFS/fast5/pcontig_019_aln_fast5

#

#load modules
cd $PBS_JOBFS
module load mcaller 
module load samtools

#

# Extract template strand reads from fast5 files using a method that saves the file path in the fastq header
#time nanopolish extract -q -t template $PBS_JOBFS_fast5 -o ${name}.fastq 
#time nanopolish index -d $PBS_JOBFS_fast5 -q ${name}.fastq	# needed for using basecalls from albacore v.2.0.0 onwards

/short/sd34/ap5514/myapps/nanopolish/0.9.0/bin/nanopolish/nanopolish index -d $PBS_JOBFS_fast5 -q $PBS_JOBFS_fastq
#

# Align fastq reads to reference assembly
out_minimap2=${PBS_JOBFS}/"minimap2"
mkdir $out_minimap2
cd $out_minimap2
echo "Mapping with minimap2"
date
time ~/myapps/minimap2/v2.7/minimap2/minimap2 -t $threads -a -x map-ont $PBS_JOBFS_ref $PBS_JOBFS_fastq | samtools view -Sb - | samtools sort -T /tmp/${name}.sorted -o ${name}.sorted.bam
echo "Done Mapping with minimap2"
date
samtools index ${name}.sorted.bam

#copy to short
cp -r $out_minimap2 ${OUTPUT}/.

#

# Run nanopolish with the following command to save a tsv file and the event values scaled towards the model

time /short/sd34/ap5514/myapps/nanopolish/0.9.0/bin/nanopolish/nanopolish eventalign -t $threads --scale-events -n -r $PBS_JOBFS_fastq -b ${out_minimap2}/${name}.sorted.bam -g $PBS_JOBFS_ref > ${name}.eventalign.tsv

#

# Run mCaller to detect m6A
# This returns a tabbed file with chromosome, read name, genomic position, position k-mer context, features, strand, and label.

time mCaller.py -m GATC -r $PBS_JOBFS_ref -d r95_twobase_model_NN_6_m6A.pkl -e ${name}.eventalign.tsv -f $PBS_JOBFS_fastq -b A 

#

# Run summary script to generate a bed file of methylated positions

time make_bed.py -f ${name}.eventalign.diffs.6 -d 15 -m 0.5

#

#Copy back to NCI
cp -r $PBS_JOBFS/* $OUTPUT

#

#delete folders in jobfs
rm -rf ${PBS_JOBFS}/*




################ Maybe keep these from tombo script ######################

#Remember to copy index file over to NCI.
cd $PBS_JOBFS/fast5
cp .${folder}.RawGenomeCorrected_000.tombo.index $OUTPUT

#Copy back to NCI
cp -r $PBS_JOBFS/alt_model $OUTPUT
rm -rf ${PBS_JOBFS}/alt_model

#
