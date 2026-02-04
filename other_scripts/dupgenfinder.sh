#!/bin/bash --login

#SBATCH --account tschwand_miya
#SBATCH --mail-type ALL
#SBATCH --mail-user kinga.kaszap@unil.ch

#SBATCH --job-name dupgenfinder
#SBATCH --output dupgenfinder_corrected.out

#SBATCH --partition cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=2:00:00
#SBATCH --export NONE


# add filepath to software before calling the function after perl 
perl software/DupGen_finder/DupGen_finder.pl -i inputs/dupgenfinder -t Tps -c Brsri -o outputs/dupgenfinder_corrected/

