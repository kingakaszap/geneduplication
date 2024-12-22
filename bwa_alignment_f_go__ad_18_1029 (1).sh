#!/bin/bash --login

#SBATCH --account tschwand_miya
#SBATCH --mail-type ALL
#SBATCH --mail-user kinga.kaszap@unil.ch

#SBATCH --job-name bwa_tps_F_Go_Ad_18-1029_lib59
#SBATCH --output bwa_tps_F_Go_Ad_18-1029_lib59.out

#SBATCH --partition cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --export NONE

module spider bwa
module load bwa

input_dir="/work/FAC/FBM/DEE/tschwand/miya/gene_duplication_files"
output_dir="/work/FAC/FBM/DEE/tschwand/miya/kkaszap/outputs/bwa"

bwa index ${input_dir}/Tps_gene_nucleotide_sequences.fasta
bwa mem ${input_dir}/Tps_gene_nucleotide_sequences.fasta \
${input_dir}/Tps_F_Go_Ad_18-1029_lib59_R1_qtrimmed.fq.gz \
${input_dir}/Tps_F_Go_Ad_18-1029_lib59_R2_qtrimmed.fq.gz \
 -t 12 > ${output_dir}/Tps_F_Go_Ad_18-1029_lib59.sam


