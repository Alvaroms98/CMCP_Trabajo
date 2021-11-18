#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00
#SBATCH --partition=mcpd

OMP_NUM_THREADS=1 ./poisson_openmp.exe 200 200