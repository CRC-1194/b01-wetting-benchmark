#!/bin/bash

pyFoamRunParameterVariation.py --list-variation case case.parameter > variation_file

Help()
{
   # Display Help
   echo 
   echo "command to run the case study"
   echo "*./reproduce_dynamicsConvergence.sh "
   echo
}

while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done


testFileName=test


#Parametrize  using pyFoam
./create-study$2.py -s $testFileName -c case -p case.parameter

#intialization of cases
for case in $testFileName*; do cd $case; ./Allclean; cd ..; done
#for case in $testFileName*; do cd $case; ./Allrun; cd ..; done

#local run
#for case in $testFileName*; do cd $case; interFlow >log.interFlow; cd ..; done

#remote run
#for case in $testFileName*; do cd $case; sbatch script.sh; cd ..; done

