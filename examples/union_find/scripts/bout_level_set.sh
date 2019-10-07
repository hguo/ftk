#!/bin/bash

#SBATCH --job-name=bout_ls_1024
#SBATCH --account=PEDAL
#SBATCH --partition=bdw
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=32
#SBATCH --output=bout_ls_1024.out
#SBATCH --error=bout_ls_1024.error
#SBATCH --time=2:00:00

# Setup My Environment
export I_MPI_FABRICS=shm:tmi

# Run My Program
srun -n 1024 /home/jiayixu/Programs/ftk/build/examples/union_find/distributed_level_set_tracking_2d -i "/home/jiayixu/Data/bout_all_normalize_new.raw" -w 425 -h 880 -t 701 --write-level-set "/home/jiayixu/Data/bout_1024.level_set" --threshold 2.5