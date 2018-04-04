#!/bin/bash
#PBS -P sd34
#PBS -q express 
#PBS -l walltime=24:00:00
#PBS -l mem=120GB
#PBS -l ncpus=16
#PBS -l jobfs=200GB

set -vx

##this is a script to do methylation-calling on contig_19 of reference genome only.

#define some variables at the start
#give this job a name
name=contig_019
short=/short/sd34/ap5514
###

#INPUTS
input=$short/methylation_calling/input
contig=pcontig_019

PBS_JOBFS_fastq=$PBS_JOBFS/fastq/${contig}_reads.fastq
PBS_JOBFS_fast5=$PBS_JOBFS/fast5/${contig}_aln_fast5.tar.gz
PBS_JOBFS_ref=$PBS_JOBFS/ref/ref_${contig}.fasta

albacore_output_fastq=$input/pcontig_019_reads.fastq
fast5_files=$input/pcontig_019_aln_fast5.tar.gz
reference_fasta=$input/ref_pcontig_019.fasta

OUTPUT=$input/${name}_output/nanopolish

threads=16
mem_size='120G'

#make the output folder
mkdir -p ${OUTPUT}

#

#now move everything to the node so we can get started
cd $PBS_JOBFS
mkdir fastq fast5 ref

cp $albacore_output_fastq $PBS_JOBFS_fastq
cp $fast5_files $PBS_JOBFS_fast5
cp $reference_fasta $PBS_JOBFS_ref

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
PBS_JOBFS_fast5=$PBS_JOBFS/fast5/${contig}_aln_fast5
#

#load modules
module load nanopolish
module load samtools
module load minimap2

#

### Data pre-processing
# make an index file linking fastq to fast5
out_index=$PBS_JOBFS/"index/"
mkdir $out_index
cd $out_index
/short/sd34/ap5514/myapps/nanopolish/0.9.0/bin/nanopolish/nanopolish index -d $PBS_JOBFS_fast5 $PBS_JOBFS_fastq 
#We get the following files: albacore_output.fastq.index, albacore_output.fastq.index.fai, albacore_output.fastq.index.gzi, and albacore_output.fastq.index.readdb 

#move output files to index folder, as they were automatically sent to the fastq folder
cp ../fastq/!(${contig}_reads.fastq) . 
#copy to short
cp -r $out_index ${OUTPUT}/.

#

#map reads to reference assembly and save as bam file
out_minimap2=${PBS_JOBFS}/"minimap2"
mkdir $out_minimap2
cd $out_minimap2
echo "Mapping with minimap2"
date
time ~/myapps/minimap2/v2.7/minimap2/minimap2 -t $threads -a -x map-ont $PBS_JOBFS_ref $PBS_JOBFS_fastq | samtools sort -@ $threads -O BAM -o ${name}.minimap2.sorted.bam
echo "Done Mapping with minimap2"
date
samtools index ${name}.minimap2.sorted.bam
#copy to short
cp -r $out_minimap2 ${OUTPUT}/.

#

# map with ngmlr
out_ngmlr=$PBS_JOBFS/"ngmlr/"
mkdir $out_ngmlr
cd $out_ngmlr
echo "Mapping with ngmlr"
date
time /home/106/ap5514/myapps/ngmlr/bin/ngmlr-0.2.6/ngmlr -t ${threads} -r $PBS_JOBFS_ref -q $PBS_JOBFS_fastq  -x ont | samtools sort -@ $threads -O BAM -o ${name}_pass.ngmlr.sorted.bam
echo "Done Mapping with ngmlr"
date

#copy to short
cp -r $out_ngmlr ${OUTPUT}/.

#

#call methylation and save to tsv file with log likelihood ratio
out_methyl=${PBS_JOBFS}/"methyl"
mkdir $out_methyl
cd $out_methyl
time /short/sd34/ap5514/myapps/nanopolish/0.9.0/bin/nanopolish/nanopolish call-methylation -t $threads -r $PBS_JOBFS_fastq -b $out_minimap2/${name}.minimap2.sorted.bam -g $PBS_JOBFS_ref > ${name}_methylation_calls_minimap2.tsv
time /short/sd34/ap5514/myapps/nanopolish/0.9.0/bin/nanopolish/nanopolish call-methylation -t $threads -r $PBS_JOBFS_fastq -b $out_ngmlr/${name}.ngmlr.sorted.bam -g $PBS_JOBFS_ref > ${name}_methylation_calls_ngmlr.tsv

#copy to short
cp -r $out_methyl ${OUTPUT}/.

#

#helper script to make tsv file showing how often each reference position was methylated
out_mfreq=${PBS_JOBFS}/"mfreq"
mkdir $out_mfreq
cd $out_mfreq
time /short/sd34/ap5514/myapps/nanopolish/0.9.0/bin/nanopolish/scripts/calculate_methylation_frequency.py -i ${out_methyl}/${name}_methylation_calls_minimap2.tsv -s > ${name}_methylation_frequency_minimap2.tsv
time /short/sd34/ap5514/myapps/nanopolish/0.9.0/bin/nanopolish/scripts/calculate_methylation_frequency.py -i ${out_methyl}/${name}_methylation_calls_ngmlr.tsv -s > ${name}_methylation_frequency_ngmlr.tsv

#copy to short
cp -r $out_mfreq ${OUTPUT}/.

#

#delete folders in jobfs 
rm -rf ${PBS_JOBFS}/* 
