#!/bin/bash
#SBATCH -A b1081
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -p b1081
#SBATCH --mem=0G
#SBATCH --job-name="surface_batch"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arianafei2024@u.northwestern.edu
#SBATCH -o "%x.o%j"

cd /projects/b1081/member_directories/aporter/Ariana_Captures/GrattonLab-General-Repo/SurfacePipeline
module load matlab
matlab -nosplash -nodesktop -singleCompThread -r surface_batch

