#!/bin/bash
#PBS -P sd34
#PBS -q normal
#PBS -l walltime=24:00:00,mem=64GB
#PBS -l ncpus=16
#PBS -l jobfs=300GB

#Script to run featureCounts on all BAM files generated during STAR RNA-mapping on the primary (p), haplotig (h), and combined (ph) genomes.

## Primary contig

#set all your variables here
short=/short/sd34/ap5514
genome=Pst_104E_v13_p_ctg

bams=$(ls ${short}/analysis/rna_mapping_output/${genome}/stringtie/*.bam)
ref=${short}/Pst_104_v13_assembly/${genome}.fa
anno=${short}/analysis/gene_model_gff_files/${genome}_combined_sorted_anno.gff3
outpath=${short}/analysis/counts

#load modules
module load subread 

#here the job starts
set -vex
mkdir $outpath
cd $PBS_JOBFS

# for genes
featureCounts -T 16 -t gene -g Name -a $anno -G $ref -o ${genome}_counts_gene.txt $bams
# for exons
featureCounts -T 16 -t exon -g Parent -a $anno -G $ref -o ${genome}_counts_exon.txt $bams
wait

mv *.txt $outpath

##

## Haplotig

genome=Pst_104E_v13_h_ctg

bams=$(ls ${short}/analysis/rna_mapping_output/${genome}/stringtie/*.bam)
ref=${short}/Pst_104_v13_assembly/${genome}.fa
anno=${short}/analysis/gene_model_gff_files/${genome}_combined_sorted_anno.gff3
outpath=${short}/analysis/counts

# for genes
featureCounts -T 16 -t gene -g Name -a $anno -G $ref -o ${genome}_counts_gene.txt $bams
# for exons
featureCounts -T 16 -t exon -g Parent -a $anno -G $ref -o ${genome}_counts_exon.txt $bams
wait

mv *.txt $outpath

##

## Combined

genome=Pst_104E_v13_ph_ctg

bams=$(ls ${short}/analysis/rna_mapping_output/${genome}/stringtie/*.bam)
ref=${short}/Pst_104_v13_assembly/${genome}.fa
anno=${short}/analysis/gene_model_gff_files/${genome}_combined_sorted_anno.gff3
outpath=${short}/analysis/counts


# for genes
featureCounts -T 16 -t gene -g Name -a $anno -G $ref -o ${genome}_counts_gene.txt $bams
# for exons
featureCounts -T 16 -t exon -g Parent -a $anno -G $ref -o ${genome}_counts_exon.txt $bams
wait

mv *.txt $outpath
