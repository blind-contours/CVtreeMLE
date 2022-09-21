#!/bin/bash
# Job name:
#SBATCH --job-name=2D_CVtreeMLE_sim
#
# Partition:
#SBATCH --partition=savio3
#

#SBATCH --qos=biostat_savio3_normal
#SBATCH --account=co_biostat

# Number of nodes for use case:
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --exclusive
#
# Wall clock limit:
#SBATCH --time 48:00:00
#
## Command(s) to run (example):
module load r/4.0.3

### for foreach+doSNOW ###
R CMD BATCH --no-save 03_run_2d_simulation.R 2D_sim.Rout
