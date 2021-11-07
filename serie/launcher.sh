#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00
#SBATCH --partition=mcpd

./poisson_origin.exe 400 400
