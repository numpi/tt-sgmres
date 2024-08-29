#!/bin/bash
#SBATCH --job-name=ex_5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=127000mb
#SBATCH --partition=gpu
#SBATCH --time=72:00:00
#SBATCH --output=%x.log

module load matlab

command="cd ~/TT-Toolbox; setup; \
  cd ~/tt-sgmres; addpath(pwd); \
  cd examples; addpath(pwd); \
  ex_5; quit"

matlab -r "$command"
