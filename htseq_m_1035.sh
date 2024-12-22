#!/bin/bash --login

#SBATCH --account tschwand_miya
#SBATCH --mail-type ALL
#SBATCH --mail-user kinga.kaszap@unil.ch

#SBATCH --job-name htseq_m_1035
#SBATCH --output htseq_m_1035.out

#SBATCH --partition cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --export NONE

input_dir="../"
output_dir_default="outputs/htseq_default"
output_dir_nonunique="outputs/htseq_nonunique"

echo started
date

module purge
dcsrsoft use arolle
module load gcc/10.4.0
module load htseq/0.11.2

htseq-count -f bam \
 -r name -s reverse \
 ${input_dir}/Tps_M_Te_Ad_18-1035_lib43_genome_namesort.bam ${input_dir}/gene_duplication_files/for_htseq/Tps.v1_spname_htseq.gtf > ${output_dir_default}/Tps_M_Te_Ad_18-1035_lib43_htseq_stranded_reverse_default.txt

htseq-count -f bam \
 -r name -s reverse \
 -m intersection-nonempty \
 --nonunique all \
 ${input_dir}/Tps_M_Te_Ad_18-1035_lib43_genome_namesort.bam ${input_dir}/gene_duplication_files/for_htseq/Tps.v1_spname_htseq.gtf > ${output_dir_nonunique}/Tps_M_Te_Ad_18-1035_lib43_htseq_stranded_reverse_nonunique.txt