#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00
#SBATCH --partition=mcpd

mpiexec ./poisson_mpi.exe 400 400