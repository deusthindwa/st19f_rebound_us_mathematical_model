#!/bin/bash
#SBATCH -J fit_ipd
#SBATCH --mem=5G
#SBATCH --cpus-per-task=4
#SBATCH --partition=week
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deus.thindwa@yale.edu
#SBATCH -o fit_ipd.out

Rscript 04_fit_model.R