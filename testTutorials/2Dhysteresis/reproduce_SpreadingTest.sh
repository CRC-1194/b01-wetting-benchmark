#!/bin/bash

pyFoamRunParameterVariation.py --list-variation case case.parameter > variation_file
fileName="test"

Help()
{
   # Display Help
   echo 
   echo "Commad to run to the case study"
   echo "./reproduce_SpreadingTest.sh"
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
./create-study.py -s $fileName -c case -p case.parameter

#intialization of cases
for case in $fileName*; do cd $case; ./Allclean; cd ..; done
for case in $fileName*; do cd $case; rm -rf 0*/*.template; cd ..; done
for case in $fileName*; do cd $case; rm -rf system/*.template; cd ..; done
#local run //Single Processor
#for case in $fileName*; do cd $case; interFlow >log.interFlow; cd ..; done

#remote run
#for case in $fileName*; do cd $case; sbatch script.sh; cd ..; done



