#!/bin/bash
# Job name:
#SBATCH --job-name=3Dr4CVtreeMLE
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
#SBATCH --time 72:00:00
#
## Command(s) to run (example):
module load r/4.0.3

### for foreach+doSNOW ###
R CMD BATCH --no-save ../03_run_3d_simulation_r4.R 3D_run4.Rout
