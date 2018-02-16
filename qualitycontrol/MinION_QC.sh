#!/bin/bash
#PBS -P xf3
#PBS -q normal 
#PBS -l walltime=48:00:00
#PBS -l mem=120GB
#PBS -l ncpus=16
#PBS -l jobfs=200GB

set -vx

#define some variables at the start
#give this job a name
name=Pst79_run1-4_1d
###
INPUT=$short/
OUTPUT=$short/
ASSEMBLY_BASE_FOLDER=$short/Pst_104_v13_assembly
gff_file=Pst_104E_v13_ph_ctg.anno.gff3
genome_file=Pst_104E_v13_ph_ctg.fa
threads=16
mem_size='120G'
LANFEAR_SCRIPTS=/home/106/ap5514/myapps/minion_qc


#make the output folder
mkdir -p ${OUTPUT}

##this is a script to do basecalling and basic QC of your reads

#now move everything to the node so we can get started
cd $PBS_JOBFS
mkdir TAR_FILES


cp FIXHERE /*tar.gz TAR_FILES/.


mkdir GENOME
cp ${ASSEMBLY_BASE_FOLDER}/${gff_file} GENOME/.
cp ${ASSEMBLY_BASE_FOLDER}/${genome_file} GENOME/.

#go ahead with unziping and basecalling
cd TAR_FILES
for x in *.tar.gz
do
len=${#x}
folder= ${x::len-7}
mv ${x} ${folder}/.
cd ${folder} 
tar -xopf ${x}&
cd $PBS_JOBFS/TAR_FILES
done

wait

rm *.tar.gz


#now caputre all the fastq for pass and fail seperately plus the summary file
cd $PBS_JOBFS

mkdir albacore_fastq 
cd albacore_fastq

cat TAR_FILES/*/out_1d/workspace/fail/*.fastq > ${name}_fail.fastq
cat TAR_FILES/*/out_1d/workspace/pass/*.fastq > ${name}_pass.fastq
cat TAR_FILES/*/out_1d/sequencing_summary.txt  > ${name}_sequencing_summary.txt

#remove TAR_FILES and zip up stuff
cd $PBS_JOBS
rm -r TAR_FILES
tar -cvzf ${name}_albacore_output.tar.gz albacore_output
mv ${name}_albacore_output.tar.gz ${OUTPUT}/.

#quickly check we are in the right spot
#run minion qc script
# minion QC
echo "Running MinION QC R script"
outmqc=$PBS_JOBFS/"minionQC"
seqsum=$PBS_JOBFS/albacore_fastq/${name}_sequencing_summary.txt


#just in case the original minin_QC.R script gets changed we port these changes
rsync -P ${LANFEAR_SCRIPTS}/minion_QC.R minion_QC.R

module load R/3.4.0
mkdir $outmqc

Rscript ./minion_QC.R $seqsum $outmqc

#now map with ngmlr
# map with ngmlr
outngmlr=$PBS_JOBFS/"ngmlr/"
mkdir $outngmlr
cd $outngmlr
echo "Mapping with ngmlr"
date

time /home/106/ap5514/myapps/ngmlr/bin/ngmlr-0.2.6/ngmlr -t ${threads} -r ${PBS_JOBFS}/GENOME/${genome_file} -q ${PBS_JOBFS}/albacore_fastq/${name}_pass.fastq -o ${name}_pass.ngmlr.out.sam -x ont
time /home/106/ap5514/myapps/ngmlr/bin/ngmlr-0.2.6/ngmlr -t ${threads} -r ${PBS_JOBFS}/GENOME/${genome_file} -q ${PBS_JOBFS}/albacore_fastq/${name}_fail.fastq -o ${name}_fail.ngmlr.out.sam -x ont
echo "Done Mapping with ngmlr"
date


#mv all the ngmlr stats analysis
#now make a bam file out of it

FIX THIS
module load samtools/1.4


module load java/jdk1.8.0_60
module load nanoplot/1.0.0
outnano=${PBS_JOBFS}/"nanoplot/"
NanoPlot --fastq_rich ${PBS_JOBFS}/albacore_fastq/${name}_pass.fastq --outdir $outnano --threads $threads --loglength
NanoPlot --fastq_rich ${PBS_JOBFS}/albacore_fastq/${name}_fail.fastq --outdir $outnano --threads $threads --loglength

for x in *.sam
do
samtools view -bS -@ $threads ${x} > ${x}.bam
samtools sort -@ $threads ${x}.bam -o ${x}.out.bam
samtools index ${x}.out.bam
/home/106/ap5514/myapps/qualimap_v2.2.1/qualimap bamqc -bam ${x}.out.bam -outdir ${PBS_JOBFS}/"qualimap_all/" -nt $threads -c --java-mem-size=$mem_size
/home/106/ap5514/myapps/qualimap_v2.2.1/qualimap bamqc -bam ${x}.out.bam -outdir ${PBS_JOBFS}/"qualimap_gff/" -gff ${PBS_JOBFS}/GENOME/${gff_file} -nt $threads -c --java-mem-size=$mem_size
rm ${x}
NanoPlot --bam ${x}.out.bam --outdir $outnano --threads $threads --loglength --prefix bam

# stats on reads > various length (thanks to @gringer here: https://bioinformatics.stackexchange.com/questions/678/get-the-mapping-statistics-of-a-single-read-$
outbam=${x}.out.bam
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>1000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_1k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>2000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_2k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>10000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_10k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>20000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_20k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>100000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_100k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>200000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_200k.txt
done


#now ngmlr is done go on with minimap


#Now map the reads with minimap2 to compare mapping with ngmlr.
#maping with minmap2
outminimap2=${PBS_JOBFS}/"minimap2/"
mkdir $outminimap2
cd $outminimap2
echo "Mapping with minimap2"
date
time 
~/myapps/minimap2/minimap2/minimap2 -ax ava-ont ${PBS_JOBFS}/albacore_fastq/${name}_pass.fastq ${PBS_JOBFS}/GENOME/${genome_file} > ${name}_pass.minimap2.out.sam
~/myapps/minimap2/minimap2/minimap2 -ax ava-ont ${PBS_JOBFS}/albacore_fastq/${name}_fail.fastq ${PBS_JOBFS}/GENOME/${genome_file} > ${name}_fail.minimap2.out.sam
echo "Done mapping with minimap2"
date



for x in *.sam
do
samtools view -bS -@ $threads ${x} > ${x}.bam
samtools sort -@ $threads ${x}.bam -o ${x}.out.bam
samtools index ${x}.out.bam
/home/106/ap5514/myapps/qualimap_v2.2.1/qualimap bamqc -bam ${x}.out.bam -outdir ${PBS_JOBFS}/"qualimap_all/" -nt $threads -c --java-mem-size=$mem_size
/home/106/ap5514/myapps/qualimap_v2.2.1/qualimap bamqc -bam ${x}.out.bam -outdir ${PBS_JOBFS}/"qualimap_gff/" -gff ${PBS_JOBFS}/GENOME/${gff_file} -nt $threads -c --java-mem-size=$mem_size
rm ${x}
NanoPlot --bam ${x}.out.bam --outdir $outnano --threads $threads --loglength --prefix bam

# stats on reads > various length (thanks to @gringer here: https://bioinformatics.stackexchange.com/questions/678/get-the-mapping-statistics-of-a-single-read-$
outbam=${x}.out.bam
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>1000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_1k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>2000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_2k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>10000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_10k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>20000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_20k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>100000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_100k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>200000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_200k.txt
done
#now move everything back to normal

cd ${PBS_JOBFS}

cp $outmqc ${OUTPUT}/.
cp $outminimap2 ${OUTPUT}/.
cp $outngmlr ${OUTPUT}/.


