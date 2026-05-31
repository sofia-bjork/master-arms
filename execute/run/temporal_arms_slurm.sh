#!/bin/bash -l

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A naiss2025-22-761 

# The name of the script is myjob
#SBATCH -J temparms1-all

# The partition
#SBATCH -p main

# 24 hours wall clock time will be given to this job
#SBATCH -t 04:00:00

# Number of nodes
#SBATCH --nodes=1

#SBATCH --mail-user=bjorkso@chalmers.se
#SBATCH --mail-type=ALL

## Set the names for the error and output files. 
## It can be smart to set a path to these to your project directory, which you can do by adding that path right after the '=' sign
#SBATCH --error=/cfs/klemming/home/b/bjorkso/Private/temporal-arms-main/execute/sbatch_jobs/job.%J.err
#SBATCH --output=/cfs/klemming/home/b/bjorkso/Private/temporal-arms-main/execute/sbatch_jobs/job.%J.out

module purge 

module load PDC/24.11
module load apptainer/1.4.0-cpeGNU-24.11

WORKDIR=/cfs/klemming/home/b/bjorkso/Private/temporal-arms-main;
TMPDIR=/cfs/klemming/scratch/b/bjorkso;
SANDBOX="$TMPDIR/temparms-sandbox"

apptainer exec --bind "$WORKDIR":/work --bind "$TMPDIR/conda_envs":/envs \
    --bind /cfs/klemming/projects/supr/naiss2025-23-46/ARMS_Genoscope_data_full/projet_DBB/ARMS-Macro/:/cfs/klemming/projects/supr/naiss2025-23-46/ARMS_Genoscope_data_full/projet_DBB/ARMS-Macro/ \
    "$SANDBOX" bash /work/execute/run/temporal_arms_run2.sh



