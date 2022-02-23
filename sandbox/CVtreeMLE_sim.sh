#!/bin/bash
# Job name:
#SBATCH --job-name=CVtreeMLE_sim
#
# Partition:
#SBATCH --partition=savio2
#

#SBATCH --qos=biostat_savio2_normal
#SBATCH --account=co_biostat

# Number of nodes for use case:
#SBATCH --nodes=1
#
# Wall clock limit:
#SBATCH --time 48:00:00
#
## Command(s) to run (example):
module load r/3.6.3

### for foreach+doSNOW ###
R CMD BATCH --no-save 03_run_simulation.R CV_treeMLE_sim.Rout
