#!/bin/bash --login

#SBATCH --account tschwand_miya
#SBATCH --mail-type ALL
#SBATCH --mail-user kinga.kaszap@unil.ch

#SBATCH --job-name sam_filter
#SBATCH --output sam_filter.out

#SBATCH --partition cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=2:30:00
#SBATCH --export NONE

module load gcc/11.4.0
module load samtools/1.17

input_dir="/work/FAC/FBM/DEE/tschwand/miya"
output_dir="/work/FAC/FBM/DEE/tschwand/miya"

samtools view -q 10 ${input_dir}/Tps_F_Go_Ad_18-1029_lib59.sam \
 > ${output_dir}/Tps_F_Go_Ad_18-1029_lib59.uniq.sam

samtools view -q 10 ${input_dir}/Tps_F_Go_Ad_18-1030_lib60.sam \
 > ${output_dir}/Tps_F_Go_Ad_18-1030_lib60.uniq.sam

samtools view -q 10 ${input_dir}/Tps_F_Go_Ad_18-1031_lib61.sam \
 > ${output_dir}/Tps_F_Go_Ad_18-1031_lib61.uniq.sam

samtools view -q 10 ${input_dir}/Tps_F_Go_Ad_18-1041_lib38.sam \
 > ${output_dir}/Tps_F_Go_Ad_18-1041_lib38.uniq.sam

samtools view -q 10 ${input_dir}/Tps_M_Te_Ad_18-1033_lib31.sam \
 > ${output_dir}/Tps_M_Te_Ad_18-1033_lib31.uniq.sam

samtools view -q 10 ${input_dir}/Tps_M_Te_Ad_18-1035_lib43.sam \
 > ${output_dir}/Tps_M_Te_Ad_18-1035_lib43.uniq.sam

samtools view -q 10 ${input_dir}/Tps_M_Te_Ad_18-1036_lib33.sam \
 > ${output_dir}/Tps_M_Te_Ad_18-1036_lib33.uniq.sam

samtools view -q 10 ${input_dir}/Tps_M_Te_Ad_18-1042_lib44.sam \
 > ${output_dir}/Tps_M_Te_Ad_18-1042_lib44.uniq.sam


