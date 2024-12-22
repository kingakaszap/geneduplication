#!/bin/bash --login

#SBATCH --account tschwand_alt_splicing
#SBATCH --mail-type ALL
#SBATCH --mail-user iulia.darolti@unil.ch

#SBATCH --job-name 01.blastp_Tps_Brsri
#SBATCH --output 01.blastp_Tps_Brsri.out

#SBATCH --partition cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=0:30:00
#SBATCH --export NONE

module load gcc/10.4.0
module load blast-plus/2.12.0

blastp -query Bacillus_protein_sequences.fasta -db Tps_db -evalue 1e-10 -max_target_seqs 5 -out Tps_Brsri_protein.blast -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -num_threads 12