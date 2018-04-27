#!/bin/bash
#PBS -P sd34
#PBS -q normal 
#PBS -l walltime=22:00:00
#PBS -l mem=128GB
#PBS -l ncpus=16
#PBS -l jobfs=420GB

set -vx

##this is a script to do methylation-calling on the whole Pst_104E reference genome.

#define some variables at the start
#give this job a name
name=ngmlr_genome
short=/short/sd34/ap5514
###

#INPUTS
input=$short/methylation_calling/input/Pst_104E/

PBS_JOBFS_fastq=$PBS_JOBFS/fastq/${name}_aln.fastq
PBS_JOBFS_fast5=$PBS_JOBFS/fast5/${name}_aln_fast5.tar.gz
PBS_JOBFS_ref=$PBS_JOBFS/ref/${name}_v13_ph_ctg.fa

albacore_output_fastq=$input/genome_aln.fastq
fast5_files=$input/genome_aln_fast5.tar.gz
reference_fasta=$short/Pst_104_v13_assembly/Pst_104E_v13_ph_ctg.fa

index1=$short/methylation_calling/input/whole_genome/Pst_104E_mc/nanopolish/index/Pst_104E_aln.fastq.index
index2=$short/methylation_calling/input/whole_genome/Pst_104E_mc/nanopolish/index/Pst_104E_aln.fastq.index.fai
index3=$short/methylation_calling/input/whole_genome/Pst_104E_mc/nanopolish/index/Pst_104E_aln.fastq.index.gzi
index4=$short/methylation_calling/input/whole_genome/Pst_104E_mc/nanopolish/index/Pst_104E_aln.fastq.index.readdb

OUTPUT=$short/methylation_calling/${name}_mc/nanopolish

threads=16
mem_size='120G'

#make the output folder
mkdir -p ${OUTPUT}

#

#now move everything to the node so we can get started
cd $PBS_JOBFS
mkdir fastq fast5 ref index

cp $albacore_output_fastq $PBS_JOBFS_fastq
cp $fast5_files $PBS_JOBFS_fast5
cp $reference_fasta $PBS_JOBFS_ref

cd $PBS_JOBFS/index/
cp $index1 .
cp $index2 .
cp $index3 .
cp $index4 .

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
PBS_JOBFS_fast5=$PBS_JOBFS/fast5/${name}_aln_fast5
#

#load modules
module load nanopolish
module load samtools

#

### Data pre-processing
# make an index file linking fastq to fast5
out_index=$PBS_JOBFS/"index/"
cd $out_index
#We get the following 4 4 4 4 files: albacore_output.fastq.index, albacore_output.fastq.index.fai, albacore_output.fastq.index.gzi, and albacore_output.fastq.index.readdb 

#copy index files to the fastq folder for later steps
cp -r ./* ../fastq/ 
#copy to short
cp -r $out_index ${OUTPUT}/.

#

#map reads to reference assembly and save as bam file

# map with ngmlr
out_ngmlr=$PBS_JOBFS/"ngmlr"
mkdir $out_ngmlr
cd $out_ngmlr
echo "Mapping with ngmlr"
date
time /home/106/ap5514/myapps/ngmlr/bin/ngmlr-0.2.6/ngmlr -t ${threads} -r $PBS_JOBFS_ref -q $PBS_JOBFS_fastq -x ont | samtools sort -@ $threads -O BAM -o ${name}.ngmlr.sorted.bam
echo "Done Mapping with ngmlr"
date
samtools index ${name}.ngmlr.sorted.bam
#copy to short
cp -r $out_ngmlr ${OUTPUT}/.

#

#call methylation and save to tsv file with log likelihood ratio
out_methyl=${PBS_JOBFS}/"methyl"
mkdir $out_methyl
cd $out_methyl
time /short/sd34/ap5514/myapps/nanopolish/0.9.0/bin/nanopolish/nanopolish call-methylation -t $threads -r $PBS_JOBFS_fastq -b $out_ngmlr/${name}.ngmlr.sorted.bam -g $PBS_JOBFS_ref > ${name}_methylation_calls_ngmlr.tsv

#copy to short
cp -r $out_methyl ${OUTPUT}/.

#

#helper script to make tsv file showing how often each reference position was methylated
out_mfreq=${PBS_JOBFS}/"mfreq"
mkdir $out_mfreq
cd $out_mfreq
time /short/sd34/ap5514/myapps/nanopolish/0.9.0/bin/nanopolish/scripts/calculate_methylation_frequency.py -i ${out_methyl}/${name}_methylation_calls_ngmlr.tsv -s > ${name}_methylation_frequency_ngmlr.tsv

#copy to short
cp -r $out_mfreq ${OUTPUT}/.

#

#delete folders in jobfs 
rm -rf ${PBS_JOBFS}/* 
