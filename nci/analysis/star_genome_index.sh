#This is a script to generate a genome index for the ph, p and h assemblies as an interactive job.
#Interactive job
qsub -I -P sd34 -q express -l ncpus=4 -l mem=30GB -l wd

#Load module
module load star

# may need to change permissions on genome fasta file in case of error: "EXITING because of INPUT ERROR: could not open genomeFastaFile:"
cd /short/sd34/ap5514/Pst_104_v13_assembly/ # go to genome directory
mkdir STAR_genome # make a genome directory for use in STAR
cd STAR_genome

# copy the genome fasta to the STAR directory
cp Pst_104E_v13_ph_ctg.fa STAR_genome/
cp Pst_104E_v13_p_ctg.fa STAR_genome/
cp Pst_104E_v13_h_ctg.fa STAR_genome/

# change all permissions for genome file
chmod a+rwx Pst_104E_v13_ph_ctg.fa
chmod a+rwx Pst_104E_v13_p_ctg.fa
chmod a+rwx Pst_104E_v13_h_ctg.fa


#Run STAR genome indexing
#Might need to put all flags on the same lines without "\" because STAR doesn't seem to recognise them and considers them another input -____- 

# ph genome assembly
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /short/sd34/ap5514/analysis/star_genome/Pst_104E_v13_ph1 \
--genomeFastaFiles /short/sd34/ap5514/Pst_104_v13_assembly/STAR_genome/Pst_104E_v13_ph_ctg.fa \
--sjdbGTFfile /short/sd34/ap5514/analysis/gene_model_gff_files/Pst_104E_v13_ph_ctg_combined_sorted_anno.gff3 \
--sjdbGTFtagExonParentTranscript Parent

# p genome assembly
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /short/sd34/ap5514/analysis/star_genome/Pst_104E_v13_p \
--genomeFastaFiles /short/sd34/ap5514/Pst_104_v13_assembly/STAR_genome/Pst_104E_v13_p_ctg.fa \
--sjdbGTFfile /short/sd34/ap5514/analysis/gene_model_gff_files/Pst_104E_v13_p_ctg_combined_sorted_anno.gff3 \
--sjdbGTFtagExonParentTranscript Parent

# h genome assembly
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /short/sd34/ap5514/analysis/star_genome/Pst_104E_v13_h \
--genomeFastaFiles /short/sd34/ap5514/Pst_104_v13_assembly/STAR_genome/Pst_104E_v13_h_ctg.fa \
--sjdbGTFfile /short/sd34/ap5514/analysis/gene_model_gff_files/Pst_104E_v13_h_ctg_combined_sorted_anno.gff3 \
--sjdbGTFtagExonParentTranscript Parent
