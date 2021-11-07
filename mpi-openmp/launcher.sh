#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --time=20:00
#SBATCH --partition=mcpd

mpiexec ./poisson_mpi_openmp.exe 400 400
