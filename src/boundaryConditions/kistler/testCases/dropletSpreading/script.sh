#!/bin/bash -l
  
#allocate 1 node for 24 hour on LB2 partition
#SBATCH -A special00005 
#SBATCH -n 96 
#SBATCH --mem-per-cpu=3700 
#SBATCH -t 120:30:00     
#SBATCH --error=error.err
#SBATCH -J CV-C1 

module load git
module load gcc/9.2.0
module load openmpi/4.0.3
module load boost
module load eigen/3.3.7
module load cmake
module load python/3.7.4
module load gnuplot
module load texlive

srun -n 32 interFlow -parallel -fileHandler collated > log
#srun -n 144 foamToVTK -parallel  > log.vtk

#reconstructPar -fileHandler collated

