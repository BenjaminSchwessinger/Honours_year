#!/bin/bash
#PBS -P xf3
#PBS -q express 
#PBS -l walltime=24:00:00,mem=126GB,ncpus=16
#PBS -l jobfs=250GB

set -vx

#load some modules
module load samtools/v1.7
module load albacore/2.2.2
module load minimap2/2.9

#define the input and output directories
seqrun=20171213_0532_131217_Dave_rapid_blind
INPUT=/short/xf3/bxs800/MinIon/data/Dave
OUTPUT=/short/xf3/bxs800/MinIon/analysis/Dave/${seqrun}
threads=16
genome_file=amil_sf_1.0.fasta


mkdir ${OUTPUT}

#move data over to JOBFS
cd $PBS_JOBFS
mkdir in
mkdir out_1d
mkdir fastq_only
mkdir nanopack
mkdir minimap2

cp $INPUT/${genome_file} $PBS_JOBFS/${genome_file}

cd in

cp -r $INPUT/*Dave_rapid_blind* .

#now unzip all the files on the node
for x in *.tar.gz
do
tar -xopzf ${x} &
done
wait

#now remove all the tarziped files
rm *.tar.gz


cd $PBS_JOBFS
time read_fast5_basecaller.py -i in -t 36 -c r94_450bps_linear.cfg -s out_1d -r -o fastq,fast5 -n 0 -q 0 --disable_pings

#now secure all the fastq data
cat ./out_1d/workspace/fail/*.fastq > $PBS_JOBFS/fastq_only/${seqrun}_fail.fastq
cat ./out_1d/workspace/pass/*.fastq > $PBS_JOBFS/fastq_only/${seqrun}_pass.fastq
cat ./out_1d/sequencing_summary.txt > $PBS_JOBFS/fastq_only/${seqrun}_sequencing_summary.txt

cp -r $PBS_JOBFS/fastq_only $OUTPUT/.

#copy everything back as tar zip and remove all the input folders
time tar czf ${seqrun}_albacore222.tar.gz out_1d
mv ${seqrun}_albacore222.tar.gz $OUTPUT/.
rm -r out_1d
rm -r $PBS_JOBFS/in
cd $PBS_JOBFS/fastq_only

#load nanopack now only
module unload albacore/2.2.2
module load nanopack/14032018

#now loop over the fastq to do stats on the run and map with minimap2
for x in *.fastq
do
out=${x}.${genome}.bam
out_paf=${x}.${genome}.paf
time NanoPlot --fastq_rich ${x} --outdir $PBS_JOBFS/nanopack -p ${x} --threads $threads --loglength
time minimap2 -t $threads -ax map-ont ${PBS_JOBFS}/${genome_file} ${x}  | samtools sort -@ $threads -O BAM -o $PBS_JOBFS/minimap2/${out}
time minimap2 -t $threads -x map-ont ${PBS_JOBFS}/${genome_file} ${x} > $PBS_JOBFS/minimap2/${out_paf}
done

cd $PBS_JOBFS/minimap2

#now loop over the mapped file and do some stats

for x in *.bam
do
samtools index ${x}
time NanoPlot --bam ${x} --outdir $PBS_JOBFS/nanopack --threads $threads --loglength --prefix ${x}
samtools view -h ${x} |     awk -F'\t' '{if((/^@/) || (length($10)>1000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_1k.txt
samtools view -h ${x} |     awk -F'\t' '{if((/^@/) || (length($10)>2000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_2k.txt
samtools view -h ${x} |     awk -F'\t' '{if((/^@/) || (length($10)>10000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_10k.txt
samtools view -h ${x} |     awk -F'\t' '{if((/^@/) || (length($10)>20000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_20k.txt
samtools view -h ${x} |     awk -F'\t' '{if((/^@/) || (length($10)>100000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > ${x}.stats_100k.txt
samtools view -h ${x} |     awk -F'\t' '{if((/^@/) || (length($10)>200000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > ${x}.stats_200k.txt
done


cp -r $PBS_JOBFS/nanopack $OUTPUT/.
cp -r $PBS_JOBFS/minimap2 $OUTPUT/.






