#!/bin/bash -l

#allocate 1 node for 24 hour on LB2 partition
#SBATCH -A special00005 
#SBATCH -n 144 
#SBATCH --mem-per-cpu=3000 
#SBATCH -t 23:30:00     
#SBATCH --error=error.err
#SBATCH -J SphSpr 

module load git
module load gcc/9.2.0
module load openmpi/4.0.3
module load boost
module load eigen/3.3.7
module load cmake
module load python/3.7.4
module load gnuplot
module load texlive

srun -n 24 interFlow -parallel -fileHandler collated > log

#reconstructPar -fileHandler collated

