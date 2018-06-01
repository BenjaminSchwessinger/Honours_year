#Interactive job
qsub -I -P sd34 -q express -l ncpus=4 -l mem=30GB -l wd

#Load module
module load star

#Run STAR genome indexing
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /short/sd34/ap5514/analysis/star_genome/Pst_104E_v13 \
--genomeFastaFiles /short/sd34/ap5514/Pst_104_v13_assembly/STAR_genome/Pst_104E_v13_ph_ctg.fa \ # copy the genome file to a new folder as STAR changes it
--sjdbGTFfile /short/sd34/ap5514/analysis/gene_model_gff_files/Pst_104E_v13_ph_combined_sorted_anno.gff3 \
--sjdbGTFtagExonParentTranscript Parent

# may need to change permissions on genome fasta file in case of error: "EXITING because of INPUT ERROR: could not open genomeFastaFile:"
cd /short/sd34/ap5514/Pst_104_v13_assembly/ # go to genome directory
mkdir STAR_genome # make a genome directory for use in STAR
cp Pst_104E_v13_ph_ctg.fa STAR_genome/ # copr the genome fasta to the STAR directory
chmod a+rwx Pst_104E_v13_ph_ctg.fa # change all permissions for genome file

