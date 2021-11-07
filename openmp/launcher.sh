#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00
#SBATCH --partition=mcpd

./poisson_openmp.exe 400 400