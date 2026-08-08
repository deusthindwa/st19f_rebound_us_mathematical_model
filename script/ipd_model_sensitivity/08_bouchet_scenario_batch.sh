#!/bin/bash
#SBATCH -J fit_ipd_scenarios
#SBATCH --array=1-3              # 3 jobs, one per scenario
#SBATCH --mem=15G
#SBATCH --cpus-per-task=4        # 4 CPUs → 4 parallel chains
#SBATCH --partition=week
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deus.thindwa@yale.edu
#SBATCH -o fit_ipd_scenarios_%a.out       # separate log per job

module load R/4.4.1-foss-2022b
module load GCC/12.2.0

#remove the ODE integrator stack size limit (prevents chain crashes on)
ulimit -s unlimited

scripts=(04a_fit_scenario1.R 04b_fit_scenario2.R 04c_fit_scenario3.R)
Rscript ${scripts[$SLURM_ARRAY_TASK_ID - 1]}