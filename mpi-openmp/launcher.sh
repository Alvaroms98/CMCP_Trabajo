#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=40:00
#SBATCH --partition=mcpd

OMP_NUM_THREADS=1 mpiexec ./poisson_mpi_openmp.exe 200 200
