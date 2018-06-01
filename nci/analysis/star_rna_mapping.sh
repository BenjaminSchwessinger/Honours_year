#!/bin/bash
#PBS -P sd34 
#PBS -q normal 
#PBS -l walltime=48:00:00,mem=88GB
#PBS -l ncpus=16
#PBS -l jobfs=420GB

#Script to run RNA mapping on all trimmed data files using STAR.

#set all your variables here
short=/short/sd34/ap5514/
genome=Pst_104E_v13_ph_ctg

#variables
DATA=${short}/trimmed_data/RNA/inhouse/
GENOME=${short}/analysis/star_genome/Pst_104E_v13/
OUTPATH=${short}/analysis/rna_mapping_output/${genome}
PBS_G=genome
GFFPATH=${short}/analysis/gene_model_gff_files/

#load modules
module load samtools
module load star
module load stringtie

#here the job starts
set -vx
mkdir ${OUTPATH}
cd $PBS_JOBFS
mkdir $PBS_G
cp $GFFPATH/${genome}_combined_sorted_anno.gff3 .
cp $GENOME/* $PBS_JOBFS/$PBS_G/
cp $DATA/[GS,HE,IT,UG]*.paired*.gz ./
wait
ls

#prepare star input
paired_1=(*_R1_*.paired.*.gz)
paired_2=(*_R2_*.paired.*.gz)
for ((i=0;i<=${#paired_1[@]}-1;i=i+1))
do
len=${#paired_1[i]}
outfolder=${paired_1[i]:0:len-13}
mkdir ${outfolder}

#run STAR aligner
time STAR --runMode alignReads --runThreadN 16 --genomeDir ${PBS_G} --readFilesCommand gunzip -c --readFilesIn ${paired_1[i]} ${paired_2[i]} --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 10000 --alignMatesGapMax 1000000  --outFileNamePrefix ${outfolder}/${outfolder} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --quantMode GeneCounts 

#run stringtie on STAR output
cd ${outfolder}
for bam in *.bam
do 
outname=$(echo $bam | cut -f1 -d .)
time stringtie -p 16 -e -G ../${genome}_combined_sorted_anno.gff3 -l ${outname} -o ${outname}.stringtie.gtf ${bam} 
done
cd ../
done
rm *.gz
mkdir stringtie
mv *001/*.stringtie.gtf stringtie
cd stringtie
ls * > stringtie_mergelist.txt
time stringtie --merge -p 16 -G ../${genome}_combined_sorted_anno.gff3 -o ${genome}.stringtie_merged.gtf stringtie_mergelist.txt
time gffcompare -r ../${genome}_combined_sorted_anno.gff3 -G -o gffcomparemerged ${genome}.stringtie_merged.gtf
mv ../*001/*.bam .
for bam in *.bam
do
label=$(echo $bam | cut -f1 -d .)
time stringtie -e -B -p 16 -G ${genome}.stringtie_merged.gtf -o ballgown/${label}/${label}_${genome}.gtf ${bam}
done
cd ${PBS_JOBFS}
rm -r $PBS_G
cp -r * ${OUTPATH}
echo "STAR done" 1>&2

