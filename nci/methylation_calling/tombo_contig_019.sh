#!/bin/bash
#PBS -P sd34
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=128GB
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
input=$short/methylation_calling/input/test_pcontig_019/input

albacore_output_fastq=$input/pcontig_019_aln.fastq
fast5_dir=$input/pcontig_019_aln_fast5.tar.gz
reference_fasta=$input/ref_pcontig_019.fasta

PBS_JOBFS_fastq=$PBS_JOBFS/fastq/${name}_aln.fastq
PBS_JOBFS_fast5=$PBS_JOBFS/fast5/${name}_aln_fast5.tar.gz
PBS_JOBFS_ref=$PBS_JOBFS/ref/ref_${name}.fasta

OUTPUT=$short/$input/${name}_mc/tombo

threads=16
mem_size='128G'

#OUTPUT
#make the output folder
mkdir -p ${OUTPUT}


#now move everything to the node so we can get started
cd $PBS_JOBFS
mkdir fastq fast5 ref alt_model de_novo

cp $albacore_output_fastq $PBS_JOBFS_fastq
cp $fast5_dir $PBS_JOBFS_fast5
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
PBS_JOBFS_fast5=$PBS_JOBFS/fast5/${name}_aln_fast5

#

#load modules
cd $PBS_JOBFS
module load tombo

#

#Pre-process reads

tombo annotate_raw_with_fastqs --fast5-basedir $PBS_JOBFS_fast5 --fastq-filenames $PBS_JOBFS_fastq

#

#Re-squiggle Algorithm (Raw Data to Genome Alignment)

tombo resquiggle $PBS_JOBFS_fast5 $PBS_JOBFS_ref --processes $threads 

### Remember to copy index file over to NCI.

#

#Comparing to an alternative 5mC and 6mA model (recommended method)
cd $PBS_JOBFS/alt_model
mkdir 5mC 6mA
tombo test_significance --fast5-basedirs $PBS_JOBFS_fast5 --alternate-bases 5mC 6mA --statistics-file-basename $name # --per-read-statistics-filename plot_stats 

#Outputs the following files to PBS_JOBS: sample.5mC.tombo.stats, sample.6mA.tombo.stats.
cp ${name}.5mC.tombo.stats 5mC
cp ${name}.6mA.tombo.stats 6mA

#Extract sequences Surrounding Modified Positions
# output the genomic sequence surrounding locations with the largest fraction of modified reads
cd $PBS_JOBFS/alt_model/5mC
tombo write_most_significant_fasta --statistics-filename sample.5mC.tombo.stats --genome-fasta $PBS_JOBFS_ref 
cd $PBS_JOBFS/alt_model/6mA 
tombo write_most_significant_fasta --statistics-filename sample.6mA.tombo.stats --genome-fasta $PBS_JOBFS_ref

#Plot it --> can't plot yet because chromosome information is needed for "--genome-locations" flag, but Pst-104E only has contigs, not chromosomes. Plotting code is commented out of test_significance line until fixed.
#tombo plot_per_read --per-read-statistics-filename plot_stats --genome-locations pcontig_019 --genome-fasta $PBS_JOBFS

#Copy back to NCI
cp -r $PBS_JOBFS/alt_model $OUTPUT

#

#Identifying de novo base modifications
cd $PBS_JOBFS/de_novo
tombo test_significance --fast5-basedirs $PBS_JOBFS_fast5 --statistics-file-basename de_novo #--per-read-statistics-basename denovo_plot_stats

#Extract sequences Surrounding Modified Positions
# output the genomic sequence surrounding locations with the largest fraction of modified reads
tombo write_most_significant_fasta --statistics-filename de_novo.tombo.stats --genome-fasta $PBS_JOBFS_ref

#Plot it --> can't plot yet because chromosome information is needed for "--genome-locations" flag, but Pst-104E only has contigs, not chromosomes. Plotting code is commented out of test_significance line until fixed.
#tombo plot_per_read --per-read-statistics-filename denovo_plot_stats --genome-locations pcontig_019 --genome-fasta $PBS_JOBFS 

#Copy back to NCI
cp -r $PBS_JOBFS/de_novo $OUTPUT
