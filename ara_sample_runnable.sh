#!/bin/bash
#SBATCH --job-name=tsunami-lab
#SBATCH --partition=s_hadoop
#SBATCH --nodes=1
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=36
#SBATCH --ntasks=1
#SBATCH --time=15:00
export OMP_NUM_THREADS=36
# could be tested in the future
# export OMP_PROC_BIND=true
# export OMP_PLACES=threads
export LD_LIBRARY_PATH=./ara/software/lib:$LD_LIBRARY_PATH
./build/tsunami_lab ./config/tsunami2d-tohoku-250.ara.yaml
