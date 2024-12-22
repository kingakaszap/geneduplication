#!/bin/bash -l

#SBATCH --account tschwand_alt_splicing
#SBATCH --mail-type ALL
#SBATCH --mail-user iulia.darolti@unil.ch

#SBATCH --job-name test_htseq_stranded_reverse
#SBATCH --output test_htseq_stranded_reverse.out

#SBATCH --partition cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=1:00:00
#SBATCH --export NONE

echo started
date

module purge
dcsrsoft use arolle
module load gcc/10.4.0
module load htseq/0.11.2

htseq-count --idattr transcript_id -f bam -r name -s no Tps_F_Go_Ad_18-1029_lib59.namesort.bam Tps.v1_spname_htseq.gtf > Tps_F_Go_Ad_18-1029_lib59_htseq_stranded_reverse.txt

echo finished
date