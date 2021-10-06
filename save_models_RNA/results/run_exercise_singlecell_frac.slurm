#!/bin/bash
#SBATCH -n 1 # Number of cores requested 
#SBATCH -J exercise # job name
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 7200 # Runtime in minutes 
#SBATCH -p batch # Partition to submit to 
#SBATCH --account=mpsnyder
#SBATCH --mem=150000 # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH -o exercise_%A.out # Standard out goes to this file 
#SBATCH -e exercise_%A.err # Standard err goes to this filehostname 

module load r/3.5
Rscript exercise_predict_fraction_singlecell.R