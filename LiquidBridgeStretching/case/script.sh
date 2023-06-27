#!/bin/bash -l
  
#allocate 1 node for 24 hour on LB2 partition
#SBATCH -A special00004 
#SBATCH -n 96 
#SBATCH --mem-per-cpu=3000 
#SBATCH -t 06:00:00     
#SBATCH --error=error.err
#SBATCH -J str_liq

module load git
module load gcc/9.2.0
module load openmpi/4.0.3
module load boost
module load eigen/3.3.7
module load cmake
module load python/3.7.4
module load gnuplot
module load texlive
source /work/projects/special00005/B01/OpenFOAM-v2212/etc/bashrc

srun -n 32 interFlow -parallel  > log
#srun -n 144 foamToVTK -parallel  > log.vtk

#reconstructPar -fileHandler collated

