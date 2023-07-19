#!/bin/bash

pyFoamRunParameterVariation.py --list-variation case case.parameter > variation_file
Help()
{
   # Display Help
   echo 
   echo "command to run the case study"
   echo "./reproduce_AdvectionTest_RefinedMesh.sh "
   echo
}

while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done


#source the OpenFOAM-v2112 environment

#Parametrize  using pyFoam
./create-study.py -s test -c case -p case.parameter

#intialization of cases
for case in test_*; do cd $case; ./Allclean; cd ..; done
for case in test_*; do cd $case; ./Allrun; cd ..; done

#local run
for case in test_*; do cd $case; interFlow >log.interFlow; cd ..; done

#remote run
#for case in test_*; do cd $case; sbatch script.sh; cd ..; done



