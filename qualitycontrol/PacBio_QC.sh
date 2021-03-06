#!/bin/bash
#PBS -P sd34 
#PBS -q express 
#PBS -l walltime=24:00:00
#PBS -l mem=120GB
#PBS -l ncpus=16
#PBS -l jobfs=200GB

set -vx

#define some variables at the start
#give this job a name
name=PacBio
short=/short/sd34/ap5514
###
INPUT=$short/raw_data/pacbio_fastq/PacBio.fastq.gz # added concatenated file, as nci job couldn't catenate individual fastq.gz
OUTPUT=$short/basecalling/quality_control/PacBio
ASSEMBLY_BASE_FOLDER=$short/Pst_104_v13_assembly
gff_file=Pst_104E_v13_ph_ctg.anno.gff3
genome_file=Pst_104E_v13_ph_ctg.fa
threads=16
mem_size='120G'


#make the output folder
mkdir -p ${OUTPUT}

##this is a script to do basecalling and basic QC of your reads

#now move everything to the node so we can get started
cd $PBS_JOBFS
mkdir FQ_GZ_FILES


cp ${INPUT} FQ_GZ_FILES/. # added basecall files to the node
# cat FQ_GZ_FILES/* >> ${name}.fastq.gz # catenate all PacBio reads

mkdir GENOME
cp ${ASSEMBLY_BASE_FOLDER}/${gff_file} GENOME/.
cp ${ASSEMBLY_BASE_FOLDER}/${genome_file} GENOME/.



module load samtools/1.7
module load java/jdk1.8.0_60
module load nanopack/1.0


#now look at the input data with nanoplot
outnano_raw=${PBS_JOBFS}/"nanopack/raw_data"
time NanoPlot --fastq ${PBS_JOBFS}/FQ_GZ_FILES/${name}.fastq.gz --outdir ${outnano_raw} --threads $threads --loglength

#folder name for NanoPlot analysis on mapped data
outnano_mapped=${PBS_JOBFS}/"nanopack/mapped"

#now map with ngmlr
# map with ngmlr
outngmlr=$PBS_JOBFS/"ngmlr/"
mkdir $outngmlr
cd $outngmlr
echo "Mapping with ngmlr"
date

time /home/106/ap5514/myapps/ngmlr/bin/ngmlr-0.2.6/ngmlr -t ${threads} -r ${PBS_JOBFS}/GENOME/${genome_file} -q ${PBS_JOBFS}/FQ_GZ_FILES/${name}.fastq.gz -o ${name}.ngmlr.out.sam
echo "Done Mapping with ngmlr"
date


#now make a bam file out of it


for x in *.sam
do
samtools view -bS -@ $threads ${x} > ${x}.bam
samtools sort -@ $threads ${x}.bam -o ${x}.out.bam
samtools index ${x}.out.bam
time /home/106/ap5514/myapps/qualimap_v2.2.1/qualimap bamqc -bam ${x}.out.bam -outdir ${PBS_JOBFS}/"qualimap_all/" -nt $threads -c --java-mem-size=$mem_size
time /home/106/ap5514/myapps/qualimap_v2.2.1/qualimap bamqc -bam ${x}.out.bam -outdir ${PBS_JOBFS}/"qualimap_gff/" -gff ${PBS_JOBFS}/GENOME/${gff_file} -nt $threads -c --java-mem-size=$mem_size
rm ${x}

NanoPlot --bam ${x}.out.bam --outdir $outnano_mapped --threads $threads --loglength --prefix ${x}

# stats on reads > various length (thanks to @gringer here: https://bioinformatics.stackexchange.com/questions/678/get-the-mapping-statistics-of-a-single-read-$
outbam=${x}.out.bam
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>1000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_1k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>2000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_2k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>10000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_10k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>20000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_20k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>100000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > ${x}.stats_100k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>200000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > ${x}.stats_200k.txt
done



#now ngmlr is done go on with minimap





#Now map the reads with minimap2 to compare mapping with ngmlr.
#mapping with minmap2
outminimap2=${PBS_JOBFS}/"minimap2/"
mkdir $outminimap2
cd $outminimap2
echo "Mapping with minimap2"
date

time #changed below to have reference before fastq 

~/myapps/minimap2/v2.7/minimap2/minimap2 -t $threads -ax map-pb ${PBS_JOBFS}/GENOME/${genome_file} ${PBS_JOBFS}/FQ_GZ_FILES/${name}.fastq.gz | samtools sort -@ $threads -O BAM -o ${name}.minimap2.out.bam

~/myapps/minimap2/v2.7/minimap2/minimap2 -x map-pb ${PBS_JOBFS}/GENOME/${genome_file} ${PBS_JOBFS}/FQ_GZ_FILES/${name}.fastq.gz > ${name}.minimap2.out.paf
echo "Done mapping with minimap2"
date



for x in *.bam
do
samtools index ${x}
time /home/106/ap5514/myapps/qualimap_v2.2.1/qualimap bamqc -bam ${x} -outdir ${PBS_JOBFS}/"qualimap_all/" -nt $threads -c --java-mem-size=$mem_size
time /home/106/ap5514/myapps/qualimap_v2.2.1/qualimap bamqc -bam ${x} -outdir ${PBS_JOBFS}/"qualimap_gff/" -gff ${PBS_JOBFS}/GENOME/${gff_file} -nt $threads -c --java-mem-size=$mem_size


NanoPlot --bam ${x} --outdir $outnano_mapped --threads $threads --loglength --prefix ${x}

# stats on reads > various length (thanks to @gringer here: https://bioinformatics.stackexchange.com/questions/678/get-the-mapping-statistics-of-a-single-read-$
outbam=${x}
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>1000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_1k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>2000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_2k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>10000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_10k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>20000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > ${x}.stats_20k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>100000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > ${x}.stats_100k.txt
samtools view -h $outbam |     awk -F'\t' '{if((/^@/) || (length($10)>200000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > ${x}.stats_200k.txt
done
#now move everything back to normal

cd ${PBS_JOBFS}

cp -r $outminimap2 ${OUTPUT}/.
cp -r $outngmlr ${OUTPUT}/.
cp -r $outnano_mapped ${OUTPUT}/.

