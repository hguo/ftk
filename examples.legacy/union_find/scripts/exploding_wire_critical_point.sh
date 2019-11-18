#!/bin/bash

#SBATCH --job-name=ew_cpt_1024
#SBATCH --account=PEDAL
#SBATCH --partition=bdw
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=32
#SBATCH --output=ew_cpt_1024.out
#SBATCH --error=ew_cpt_1024.error
#SBATCH --time=2:00:00

# Setup My Environment
export I_MPI_FABRICS=shm:tmi

# Run My Program
srun -n 1024 /home/jiayixu/Programs/ftk/build/examples/union_find/distributed_critical_point_tracking_2d -i "/home/jiayixu/Data/left_fg_ftk.raw" -w 384 -h 384 -t 4745 --write-traj "/home/jiayixu/Data/exploding_wires_1024.traj"