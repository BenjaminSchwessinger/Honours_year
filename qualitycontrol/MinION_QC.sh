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
name=Pst79_run1_1d2_0245
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
tar -xopf ${x}&
done

wait

rm *.tar.gz

cd $PBS_JOBFS
mkdir albacore_fastq
cp albacore_output/sequencing_summary.txt albacore_fastq/.
cat albacore_output/workspace/*.fastq > albacore_fastq/${name}.all.fastq

#remove TAR_FILES and zip up stuff
rm -r TAR_FILES
tar -cvzf albacore_output.tar.gz albacore_output
rm -r albacore_output

#now move everything from this step down already
mv albacore_output.tar.gz ${OUTPUT}/.

cp -r albacore_fastq ${OUTPUT}/.

#quickly check we are in the right spot
#run minion qc script
# minion QC
echo "Running MinION QC R script"
outmqc=$PBS_JOBFS/"minionQC"
seqsum=$PBS_JOBFS/albacore_fastq/"sequencing_summary.txt"


#just in case the original minin_QC.R script gets changed we port these changes
rsync -P ${LANFEAR_SCRIPTS}/minion_qc/minion_QC.R minion_QC.R

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
time ~/myapps/ngmlr/0.2.6/ngmlr/bin/ngmlr-0.2.6/ngmlr -t ${threads} -r ${PBS_JOBFS}/GENOME/${genome_file} -q ${PBS_JOBFS}/albacore_fastq/${name}.all.fastq -o ${name}.ngmlr.out.sam -x ont
echo "Done Mapping with ngmlr"
date


#mv all the ngmlr stats analysis
#now make a bam file out of it
module load samtools/1.4
samtools view -bS -@ $threads ${name}.ngmlr.out.sam > ${name}.out.bam
samtools sort -@ $threads ${name}.out.bam -o ${name}.out.bam
samtools index ${name}.out.bam

rm ${name}.ngmlr.out.sam

module load java/jdk1.8.0_60
#now do the qualimap
/home/106/rn5305/downloads/qualimap_v2.2.1/qualimap bamqc -bam ${name}.out.bam -outdir ${PBS_JOBFS}/"qualimap_all/" -nt $threads -c --java-mem-size=$mem_size
/home/106/rn5305/downloads/qualimap_v2.2.1/qualimap bamqc -bam ${name}.out.bam -outdir ${PBS_JOBFS}/"qualimap_gff/" -gff ${PBS_JOBFS}/GENOME/${gff_file} -nt $threads -c --java-mem-size=$mem_size


# run nanoplot
#module to load for nanoplot
module unload albacore/1.2.6
module load nanoplot/0.16.2

outnano=${PBS_JOBFS}/"nanoplot/"
NanoPlot --fastq_rich ${PBS_JOBFS}/${name}.all.fastq --outdir $outnano --threads $threads --loglength
NanoPlot --bam ${name}.out.bam --outdir $outnano --threads $threads --loglength --prefix bam

# stats on reads > various length (thanks to @gringer here: https://bioinformatics.stackexchange.com/questions/678/get-the-mapping-statistics-of-a-single-read-$
outbam=$outngmlr/${name}.out.bam
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>1000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_1k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>2000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_2k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>10000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_10k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>20000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_20k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>100000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_100k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>200000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_200k.txt

#now ngmlr is done go on with minimap


#Now map the reads with minimap2 to compare mapping with ngmlr.
#maping with minmap2
outminimap2=${PBS_JOBFS}/"minimap2/"
mkdir $outminimap2
cd $outminimap2
echo "Mapping with minimap2"
date
time 
~/myapps/minimap2/minimap2/minimap2 -ax map10k ${PBS_JOBFS}/albacore_fastq/${name}.all.fastq ${PBS_JOBFS}/GENOME/${genome_file} > ${name}.minimap2.out.sam
echo "Done mapping with minimap2"
date


###Now process the minmap2 sam file
#now make a bam file out of it
samtools view -bS -@ $threads ${name}.minimap2.out.sam > ${name}.minimap2.out.bam
samtools sort -@ $threads ${name}.minimap2.out.bam -o ${name}.minimap2.out.bam
samtools index ${name}.minimap2.out.bam

rm ${name}.minimap2.out.sam

#now do the qualimap
/home/106/rn5305/downloads/qualimap_v2.2.1/qualimap bamqc -bam ${name}.minimap2.out.bam -outdir ${PBS_JOBFS}/"minimpa2.qualimap_all/" -nt $threads -c --java-mem-size=$mem_size

/home/106/rn5305/downloads/qualimap_v2.2.1/qualimap bamqc -bam ${name}.minimap2.out.bam -outdir ${PBS_JOBFS}/"minimap2.qualimap_gff/" -gff ${PBS_JOBFS}/GENOME/${gff_file} -nt $threads -c --java-mem-size=$mem_size

#run nanoplot
#module to load for nanoplot
module load nanoplot

#NanoPlot --bam ${name}.minimap2.out.bam --outdir $outnano --threads $threads --loglength --prefix bam

# stats on reads > various length (thanks to @gringer here: https://bioinformatics.stackexchange.com/questions/678/get-the-mapping-statistics-of-a-single-read-$
outbam=$outminimap2/${name}.minimap2.out.bam
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>1000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_1k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>2000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_2k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>10000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_10k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>20000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_20k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>100000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_100k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>200000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_200k.txt



#now move everything back to normal
cd ${PBS_JOBFS}
mv * ${OUTPUT}/.


