#!/bin/bash
#SBATCH -J fit1
#SBATCH --mem=10G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deus.thindwa@yale.edu
#SBATCH -o fit1.out

Rscript script/IPD/07_lhsParrun.R