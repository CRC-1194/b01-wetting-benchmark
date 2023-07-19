#!/bin/bash -l

#allocate 1 node for 24 hour on LB2 partition
#SBATCH --nodes=1 --ntasks=2 --time=20:30:00
#SBATCH --mem-per-cpu=1750
###stderror file
#SBATCH --error=error.err
#SBATCH -A special00005
#SBATCH -J Adv/Ref 
#SBATCH --mail-user=hassan.asghar@tu-darmstadt.de
#SBATCH --mail-type=ALL

module load git
module load gcc/9.2.0
module load openmpi/4.0.3
module load boost
module load eigen/3.3.7
module load cmake
module load python/3.7.4
module load gnuplot
module load texlive
#source /work/projects/project01456/openfoam-v2012/etc/bashrc
#source /work/projects/project01456/k71/openfoam/etc/bashrc 
source /work/projects/special00005/B01/OpenFOAM-v2112/etc/bashrc
#source $HOME/.bashrc
#cp 0.orig/alpha.water 0/
#setFields

#If !decomposed
#rm -rf processor*
#decomposePar

interFlow > log
#srun -n 24 interFlow -parallel > log
#reconstructPar

#interIsoFoam
#for case in advectionTestGradedMesh_000*; do cd $case; interIsoFoam > log.interIsoFoam && cd .. ; done
