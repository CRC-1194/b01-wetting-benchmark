import numpy as np
import math
import ast
import os 
import csv

sin =np.sin
cos= np.cos
pi = np.pi

class height_vector:
    x_coord = 0.0
    y_coord = 0.0
    z_coord = 0.0  

class Funcs:
    #Formula from the Patel et al. paper
    @classmethod
    def get_V(self, a, b,R0):
        return (pi/3*(R0*sin(a)/sin(b))**(3) * (1+cos(b))**(2) *(2-cos(b))) \
               +(pi/3*R0**(3) * (1+cos(a))**(2) *(2-cos(a))) \
               -(4/3*pi*R0**3) 

    # Solving the geometrical relations from the Patel et.al. paper to get contact radius r_f and height h_f
    #using bisection method to get a,b that minmize |V^{exact} - V(a,b, R_0) |
    @classmethod
    def getAngles(self, thetaDeg, V, R0):
        thetaDeg = 180.0-thetaDeg
        a, b = np.radians(thetaDeg/2), np.radians(thetaDeg/2)

        isNewLess = False
        isNewLessPrev = False
        if(self.get_V(a, b,R0) < V):
            isNewLess = True
            isNewLessPrev = True

        count =0
        boolSub = True
        deltaV = V - self.get_V(a,b,R0)
        deltaVPrev = deltaV

        while (isNewLess == isNewLessPrev ):
            #print(math.fabs(get_V(a, b)))
            if(count ==0):
                count+=0.1
                a, b = np.radians(((thetaDeg-count)/2) ), np.radians(((thetaDeg+count)/2))
                deltaV = V - self.get_V(a,b,R0)
                if (deltaV<deltaVPrev):
                    boolSub = True
                else:
                    boolSub = False
                continue

            deltaV = V - self.get_V(a,b,R0)
            #check for Delta_V 
            if (boolSub): #getting closer to the solution -> continue decreasing a and increasing b 
                a, b = np.radians(((thetaDeg-count)/2)), np.radians(((thetaDeg+count)/2))
                deltaVPrev =deltaV
            else:
                a, b = np.radians(((thetaDeg+count)/2)), np.radians(((thetaDeg-count)/2))
                deltaVPrev =deltaV
            
            if(self.get_V(a, b,R0) < V):
                isNewLess = True
            elif(self.get_V(a, b,R0) > V):
                isNewLess = False
            count+=0.1
        
        hf= ((R0*sin(a)/sin(b)*(1+cos(b))) - (R0*(1-cos(a))))/R0
        rf = sin(a)
        return hf, rf, (rf/sin(b))


    # Calculate the height of the equilibrium shape from the isoAlpha.vtk file and write to postProcessing/height.csv file
    @classmethod
    def writeHeightFile(self,calcHeights,vtk_folders, heightFileNames):
        if(calcHeights):
            for idx, vtk_folder in enumerate(vtk_folders):
                #getting the time step folder names
                sub_folders = [name for name in os.listdir(vtk_folder) if os.path.isdir(os.path.join(vtk_folder, name))]
                sub_folders.sort()
                numberOfCoordinates = 0.0
                ifCoordArray = 0 # just to take POINTS array
                heightObject = height_vector()
                counter =0

                #file-content
                # header = ['Time(s)', 'Height(m)']
                with open(heightFileNames[idx], 'w') as f:
                    writer = csv.writer(f)

                for sub_folder in sub_folders:
                    heightObject.z_coord = 0
                    with open(vtk_folder + sub_folder +"/isoAlpha.vtk") as reader:
                        for line in reader:
                            if "POINTS" in line: #start of coordinates
                                #print(line, end='')
                                breakUp = line.split()
                                numberOfCoordinates = int(breakUp[1])
                                ifCoordArray = 1;

                            if "\n" not in line[0] and ifCoordArray ==1 and "POINTS" not in line:
                                coordBreakUp = line.split()

                                if (counter ==0):
                                    heightObject.z_coord = float (coordBreakUp[2])
                                    counter =1
                                else:
                                #check on Z coordinate
                                    if (len(coordBreakUp)>3):
                                        heightObject.z_coord = max (heightObject.z_coord, float(coordBreakUp[2]))
                                        heightObject.z_coord = max (heightObject.z_coord, float(coordBreakUp[5]))
                                    else:
                                        heightObject.z_coord = max (heightObject.z_coord, float(coordBreakUp[2]))

                            if "\n" in line[0] and ifCoordArray ==1:
                                ifCoordArray =2

                    data_to_write_to_csv = [sub_folder, str(heightObject.z_coord)]

                    with open(heightFileNames[idx], 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow(data_to_write_to_csv)
