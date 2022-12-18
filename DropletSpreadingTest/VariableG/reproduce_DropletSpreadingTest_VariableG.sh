#!/bin/bash

pyFoamRunParameterVariation.py --list-variation case case.parameter > variation_file

Help()
{
   # Display Help
   echo 
   echo "command to run the case study"
   echo "./reproduce_DropletSpreadingTest_VariableG.sh [fluid] [contact_angle]"
   echo
   echo "E.g., for water_gylcerol, run ./reproduce_DropletSpreadingTest_VariableG.sh water_glycerol 50"
   echo "options:"
   echo "fluid   water_glycerol | water"
}

while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done


cp  case.parameter.backup case.parameter

#./sh [fluid] [contact_angle]

testFileName=water_glycerol
if [ $1 = "water_glycerol" ]
then
	sed -i 's/aNu/2.507e-05/g' case.parameter
	sed -i 's/aRho/1194.9/g' case.parameter
	sed -i 's/bNu/1.516e-05/g' case.parameter
	sed -i 's/bRho/1.2040/g' case.parameter
	sed -i 's/surTen/0.0635/g' case.parameter
	sed -i 's/eTime/0.05/g' case.parameter
	testFileName="wgtest"$2

	#source the OpenFOAM-v2112 environment

	#Parametrize  using pyFoam
	./create-study$2.py -s $testFileName -c case -p case.parameter

	#intialization of cases
	for case in $testFileName*; do cd $case; ./Allclean; cd ..; done
	for case in $testFileName*; do cd $case; ./Allrun; cd ..; done

	#local run
	#for case in $testFileName*; do cd $case; interFlow >log.interFlow; cd ..; done

	#remote run
	for case in $testFileName*; do cd $case; sbatch script.sh; cd ..; done

elif [ $1 = "water" ]
then
	sed -i 's/aNu/1.0e-06/g' case.parameter
        sed -i 's/aRho/1000/g' case.parameter
        sed -i 's/bNu/1.48e-05/g' case.parameter
        sed -i 's/bRho/1.0/g' case.parameter
        sed -i 's/surTen/0.072/g' case.parameter 
        sed -i 's/eTime/0.25/g' case.parameter
	testFileName="wtest"$2

	#source the OpenFOAM-v2112 environment

        #Parametrize  using pyFoam
        ./create-study$2.py -s $testFileName -c case -p case.parameter

        #intialization of cases
        for case in $testFileName*; do cd $case; ./Allclean; cd ..; done
        for case in $testFileName*; do cd $case; ./Allrun; cd ..; done

        #local run
        #for case in $testFileName*; do cd $case; interFlow >log.interFlow; cd ..; done

        #remote run
        for case in $testFileName*; do cd $case; sbatch script.sh; cd ..; done

fi
#echo $1
#echo $testFileName


