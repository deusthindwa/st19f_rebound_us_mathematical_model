#!/bin/bash
#SBATCH -J fit_ipd_scenarios
#SBATCH --array=1-3              # 3 jobs, one per scenario
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2        # 2 CPUs → 2 parallel chains
#SBATCH --partition=long
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deus.thindwa@yale.edu
#SBATCH -o fit_ipd_scenarios_%a.out       # separate log per job

scripts=(04a_fit_scenario1.R 04b_fit_scenario2.R 04c_fit_scenario3.R)
Rscript ${scripts[$SLURM_ARRAY_TASK_ID - 1]}