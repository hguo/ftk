#!/bin/bash

#SBATCH --job-name=ew_cpt_1
#SBATCH --account=PEDAL
#SBATCH --partition=bdw
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=ew_cpt_1.out
#SBATCH --error=ew_cpt_1.error
#SBATCH --time=2:00:00

# Setup My Environment
export I_MPI_FABRICS=shm:tmi

# Run My Program
srun -n 1 /home/jiayixu/Programs/ftk/build/examples/union_find/distributed_critical_point_tracking_2d --read-traj "/home/jiayixu/Data/exploding_wires_1024.traj" --write-traj-vtk "/home/jiayixu/Data/exploding_wires_1024.traj.vtp"