#!/bin/bash
#SBATCH -n 1 # Number of cores requested 
#SBATCH -J combined_microglia # job name
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 7200 # Runtime in minutes 
#SBATCH -p batch # Partition to submit to 
#SBATCH --account=mpsnyder
#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH -o combined_mural_%A.out # Standard out goes to this file 
#SBATCH -e combined_mural_%A.err # Standard err goes to this filehostname 

module load r/3.5
Rscript CLUSTER_PredictCombined_Mural.R