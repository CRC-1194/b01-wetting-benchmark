#!/bin/bash


testFileName='test'

cp  case.parameter.backup case.parameter
pyFoamRunParameterVariation.py --list-variations case case.parameter > variation_file
#./.sh [fluid]

        #Parametrize  using pyFoam
./create-study.py -s $testFileName -c case -p case.parameter

        #intialization of cases
for case in $testFileName*; do cd $case; ./Allclean; cd ..; done
for case in $testFileName*; do cd $case; ./Allrun; cd ..; done

        #local run
#for case in $testFileName*; do cd $case; interFlow >log.interFlow; cd ..; done

        #remote run
for case in $testFileName*; do cd $case; sbatch script.sh; cd ..; done

