#!/bin/bash

pyFoamRunParameterVariation.py --list-variation case case.parameter > variation_file
nProc='40'

Help()
{
   # Display Help
   echo 
   echo "Command to run the case studies:"
   echo "./reproduce_2DCapillaryRiseTest.sh [U_boundarycondition] [Omega]"	
   echo "For example:"
   echo "./reproduce_2DCapillaryRiseTest.sh slip omega1"
   echo "options:"
   echo " U_boundarycondition:  slip|noslip"
   echo " Omega:    omega1| omega0.1 |omega0.5"
   echo 
}

while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done


cp  case.parameter.backup case.parameter
#./reproduce_2DCapillaryRiseTest.sh [U_boundarycondition (slip/noslip)] [Omega (omega1/ omega.1 / omega.5)]

test=testOmega1
if [ $1 = "noslip" ]
then
	cp  case/0.org/U.noslip.template case/0.org/U.template
elif [ $1 = "slip" ]
then
	cp  case/0.org/U.slip.template case/0.org/U.template
fi

if [ $2 = "omega1" ]
then
	sed -i 's/gy/4.17/g' case.parameter
        sed -i 's/surfTen/0.04/g' case.parameter
	sed -i 's/alpha1rho/83.1/g' case.parameter
        sed -i 's/maxca/0.029/g' case.parameter
	sed -i 's/etime/0.7/g' case.parameter
	sed -i 's/Slip/0.0001/g' case.parameter
	test="testOmega1"
        if [ $1 = "noslip" ]
	then
		test="testU0"
	fi
	#source the OpenFOAM-v2112 environment
	#Parametrize  using pyFoam
	./create-study.py -s $test -c case -p case.parameter

	#intialization of cases
	for case in $test*; do cd $case; ./Allclean; cd ..; done
	for case in $test*; do cd $case; ./Allrun; cd ..; done

	#local run
	#for case in $test*; do cd $case; mpirun -n 4 interFlow -parallel > log.interFlow; cd ..; done

	#remote run
	for case in $test*; do cd $case; sbatch script.sh; cd ..; done

elif [ $2 = "omega0.1" ]
then
        sed -i 's/gy/1.04/g' case.parameter
        sed -i 's/surfTen/0.2/g' case.parameter
        sed -i 's/alpha1rho/1663.1/g' case.parameter
        sed -i 's/maxca/0.0033/g' case.parameter
        sed -i 's/etime/10/g' case.parameter
        sed -i 's/Slip/0.001/g' case.parameter
	test="testOmega0.1"

        #source the OpenFOAM-v2112 environment
        #Parametrize  using pyFoam
        ./create-study.py -s $test -c case -p case.parameter

        #intialization of cases
        for case in $test*; do cd $case; ./Allclean; cd ..; done
        for case in $test*; do cd $case; ./Allrun; cd ..; done

        #local run
        #for case in $test*; do cd $case; mpirun -n 4 interFlow -parallel > log.interFlow; cd ..; done

        #remote run
        for case in $test*; do cd $case; sbatch script.sh; cd ..; done

elif [ $2 = "omega0.5" ]
then
        sed -i 's/gy/6.51/g' case.parameter
        sed -i 's/surfTen/0.1/g' case.parameter
        sed -i 's/alpha1rho/133.0/g' case.parameter
        sed -i 's/maxca/0.015/g' case.parameter
	sed -i 's/etime/1.2/g' case.parameter
  	sed -i 's/Slip/0.001/g' case.parameter
        test="testOmega0.5"

        #source the OpenFOAM-v2112 environment
        #Parametrize  using pyFoam
        ./create-study.py -s $test -c case -p case.parameter

        #intialization of cases
        for case in $test*; do cd $case; ./Allclean; cd ..; done
        for case in $test*; do cd $case; ./Allrun; cd ..; done

        #local run
        #for case in $test*; do cd $case; mpirun -n 4 interFlow -parallel > log.interFlow; cd ..; done

        #remote run
        for case in $test*; do cd $case; sbatch script.sh; cd ..; done
fi
