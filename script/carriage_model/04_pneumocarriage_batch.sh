#!/bin/bash
#SBATCH -J fit_carriage
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --partition=week
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deus.thindwa@yale.edu
#SBATCH -o fit_carriage.out

Rscript 01_pneumocarriage_run.R