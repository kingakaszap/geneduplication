#!/bin/bash --login

#SBATCH --account tschwand_miya
#SBATCH --mail-type ALL
#SBATCH --mail-user kinga.kaszap@unil.ch

#SBATCH --job-name bwa_tps_multiple
#SBATCH --output bwa_tps_multiple.out

#SBATCH --partition cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=5:30:00
#SBATCH --export NONE

module load gcc/11.4.0
module load bwa/0.7.17

input_dir="/work/FAC/FBM/DEE/tschwand/miya"
output_dir="/work/FAC/FBM/DEE/tschwand/miya"

bwa mem ${input_dir}/Tps_gene_nucleotide_sequences.fasta \
${input_dir}/Tps_F_Go_Ad_18-1030_lib60_R1_qtrimmed.fq.gz \
${input_dir}/Tps_F_Go_Ad_18-1030_lib60_R2_qtrimmed.fq.gz \
 -t 12 > ${output_dir}/Tps_F_Go_Ad_18-1030_lib60.sam

bwa mem ${input_dir}/Tps_gene_nucleotide_sequences.fasta \
${input_dir}/Tps_F_Go_Ad_18-1031_lib61_R1_qtrimmed.fq.gz \
${input_dir}/Tps_F_Go_Ad_18-1031_lib61_R2_qtrimmed.fq.gz \
 -t 12 > ${output_dir}/Tps_F_Go_Ad_18-1031_lib61.sam

bwa mem ${input_dir}/Tps_gene_nucleotide_sequences.fasta \
${input_dir}/Tps_F_Go_Ad_18-1041_lib38_R1_qtrimmed.fq.gz \
${input_dir}/Tps_F_Go_Ad_18-1041_lib38_R2_qtrimmed.fq.gz \
 -t 12 > ${output_dir}/Tps_F_Go_Ad_18-1041_lib38.sam

bwa mem ${input_dir}/Tps_gene_nucleotide_sequences.fasta \
${input_dir}/Tps_M_Te_Ad_18-1033_lib31_R1_qtrimmed.fq.gz \
${input_dir}/Tps_M_Te_Ad_18-1033_lib31_R2_qtrimmed.fq.gz \
 -t 12 > ${output_dir}/Tps_M_Te_Ad_18-1033_lib31.sam

bwa mem ${input_dir}/Tps_gene_nucleotide_sequences.fasta \
${input_dir}/Tps_M_Te_Ad_18-1035_lib43_R1_qtrimmed.fq.gz \
${input_dir}/Tps_M_Te_Ad_18-1035_lib43_R2_qtrimmed.fq.gz \
 -t 12 > ${output_dir}/Tps_M_Te_Ad_18-1035_lib43.sam

bwa mem ${input_dir}/Tps_gene_nucleotide_sequences.fasta \
${input_dir}/Tps_M_Te_Ad_18-1036_lib33_R1_qtrimmed.fq.gz \
${input_dir}/Tps_M_Te_Ad_18-1036_lib33_R2_qtrimmed.fq.gz \
 -t 12 > ${output_dir}/Tps_M_Te_Ad_18-1036_lib33.sam

bwa mem ${input_dir}/Tps_gene_nucleotide_sequences.fasta \
${input_dir}/Tps_M_Te_Ad_18-1042_lib44_R1_qtrimmed.fq.gz \
${input_dir}/Tps_M_Te_Ad_18-1042_lib44_R2_qtrimmed.fq.gz \
 -t 12 > ${output_dir}/Tps_M_Te_Ad_18-1042_lib44.sam