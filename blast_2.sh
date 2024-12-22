#!/bin/bash --login

#SBATCH --account tschwand_miya
#SBATCH --mail-type ALL
#SBATCH --mail-user kinga.kaszap@unil.ch

#SBATCH --job-name blast_for_dupgenefinder
#SBATCH --output blast_for_dupgenefinder.out

#SBATCH --partition cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=2:00:00
#SBATCH --export NONE


module load gcc/11.4.0
module load blast-plus/2.14.1


makeblastdb -in Tps_gene_protein_sequences.fasta -input_type fasta -dbtype prot -title Tps_db -out Tps_db

blastp -query Tps_gene_protein_sequences.fasta -db Tps_db -evalue 1e-10 -max_target_seqs 5 -out Tps_Tps_protein.blast -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -num_threads 12

blastp -query Bacillus_protein_sequences.fasta -db Tps_db -evalue 1e-10 -max_target_seqs 5 -out Tps_Brsri_protein.blast -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -num_threads 12

