#!/bin/bash -l

#allocate 1 node for 24 hour on LB2 partition
#SBATCH --nodes=1 --ntasks=2 --time=23:30:00
#SBATCH --mem-per-cpu=1750
###stderror file
#SBATCH --error=error.err
#SBATCH -J Adv/Uni 

module load git
module load gcc/9.2.0
module load openmpi/4.0.3
module load boost
module load eigen/3.3.7
module load cmake
module load python/3.7.4
module load gnuplot
module load texlive

interFlow > log
#srun -n 24 interFlow -parallel > log
#reconstructPar

