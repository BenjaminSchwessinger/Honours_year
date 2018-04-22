#!/bin/bash
#PBS -P sd34
#PBS -q normal 
#PBS -l walltime=30:00:00,mem=126GB,ncpus=16
#PBS -l jobfs=420GB

set -vx

#define the input and output directories
seqrun=Pst79_1
INPUT=/short/sd34/ap5514/raw_data/${seqrun}_zipped
OUTPUT=/short/sd34/ap5514/basecalling/basecalled_albacore2110/${seqrun}/

#move data over to JOBFS

cd $PBS_JOBFS
mkdir in
mkdir out_1d
mkdir out_1d2
cd in
var=($(ls ${INPUT} --sort=size))
len=${#var[@]}              
for ((i=0; i<$len; i=i+3))                                     
do                                             
cp ${INPUT}/${var[i]} .
done

#now unzip all the files on the node
for x in *.tar.gz
do
tar -xopzf ${x} &
done
wait

#now remove all the tarziped files
rm *.tar.gz

#now do the basecalling 1D basecalling
module load albacore/2.1.10

cd $PBS_JOBFS
time read_fast5_basecaller.py -i in -t 36 -c r95_450bps_linear.cfg -s out_1d -r -o fastq,fast5 -n 0 -q 0 --disable_pings

time tar -czf $OUTPUT/${seqrun}_1d_albacore2110.tar.gz out_1d

rm -r out_1d
rm -r in
