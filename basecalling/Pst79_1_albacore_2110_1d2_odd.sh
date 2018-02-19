#!/bin/bash
#PBS -P sd34
#PBS -q normal 
#PBS -l walltime=48:00:00,mem=126GB,ncpus=16
#PBS -l jobfs=400GB

set -vx

#define the input and output directories
INPUT=/short/sd34/ap5514/raw_data/Pst79_1_zipped
OUTPUT=/short/sd34/ap5514/basecalling/basecalled_albacore2110/Pst79_1_1d2

#move data over to JOBFS

cd $PBS_JOBFS
mkdir in
mkdir out_1d
mkdir out_1d2
cd in
var=($(ls ${INPUT} --sort=size))
len=${#var[@]}
for ((i=1; i<$len; i=i+2))
do
cp ${INPUT}/${var[i]} .
done

#now unzip all the files on the node
for x in *.tar.gz
do
tar -xopzf ${x} &
done
wait

#now remove all the tarzipped files
rm *.tar.gz

#now do the basecalling 1D basecalling
module load albacore/2.1.10
cd $PBS_JOBFS


#now do the 1d2 basecalling
time full_1dsq_basecaller.py -i in -t 20 -k SQK-LSK308 -f FLO-MIN107 -s out_1d2 -r -o fastq -n 0 -q 0 --disable_pings

time tar czf Pst79_1_zipped_1d2_albacore2110.tar.gz out_1d2

mv Pst79_1_zipped_1d2_albacore2110.tar.gz $OUTPUT

rm -r out_1d2   
rm -r in
