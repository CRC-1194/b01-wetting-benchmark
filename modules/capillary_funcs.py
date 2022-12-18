import numpy as np
import math
import ast
import vtk
import os
import csv
class height_vector:
    x_coord = 0.0
    y_coord = 0.0
    z_coord = 0.0

class Funcs:
#calculates the height of the interface
# by reading .vtk surface files
# and write it in .csv file
    @classmethod
    def writeHeightFile(self, calcHeights,vtk_folders, heightFileNames):
        if(calcHeights):
            for idx, vtk_folder in enumerate(vtk_folders):
                #getting the time step folder names
                sub_folders = [name for name in os.listdir(vtk_folder) if os.path.isdir(os.path.join(vtk_folder, name))]
                sub_folders.sort()

                #file-content
                # header = ['Time(s)', 'Height(m)']
                with open(heightFileNames[idx], 'w') as f:
                    writer = csv.writer(f)

                for sub_folder in sub_folders:
                    heightObject = height_vector()
                    heightObject.z_coord = 0
                    numberOfCoordinates = 0.0
                    ifCoordArray = 0 # just to take POINTS array
                    counter =0
                    with open(vtk_folder + sub_folder +"/isoAlpha.vtk") as reader:
                        for line in reader:
                            if "POINTS" in line: #start of coordinates
                                breakUp = line.split()
                                numberOfCoordinates = int(breakUp[1])
                                ifCoordArray = 1;
 
                            elif "\n" not in line[0] and ifCoordArray ==1 and "POINTS" not in line:
                                coordBreakUp = line.split()
                                coordPerLine = int(len(coordBreakUp) / 3)
                                if (counter ==0):
                                    heightObject.z_coord = float (coordBreakUp[1])
                                    counter =1
                                if(counter==1):
                                #check on Z coordinate
                                    for i in range(0, coordPerLine):
                                        coordinateIndex = ((int(len(coordBreakUp)-1)%3)-1) + (i*3)
                                        heightObject.z_coord = min (heightObject.z_coord, float(coordBreakUp[coordinateIndex]))

                            if "\n" in line[0] and ifCoordArray ==1:
                                ifCoordArray =2
    
                    data_to_write_to_csv = [sub_folder, str(heightObject.z_coord*1000.0)]
                
                    with open(heightFileNames[idx], 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow(data_to_write_to_csv)    

    @classmethod
    #Function gives dictionary mapping variation number to parameter vector
    def fileMapping(self,Omega):
        cwd = os.getcwd()
        casefolders = [folder for folder in os.listdir(cwd) if Omega in folder] 
        casefolders.sort()
        mapNumber = [int(num.split('_')[1]) for num in  casefolders ]

        var_map = {}
        var_lines = open("variation_file").readlines()
        nX =[]

        for line in var_lines:
            # Skip lines without mapping
            if '{' not in line:
                continue
            var_num = int(line.split()[1])
            if(var_num in mapNumber):
                dict_start = line.find('{')
        # Mappings in variation file can directly be interpreted by Python
                var_map[var_num] = ast.literal_eval(line[dict_start:-1])

        for key, value in var_map.items():
            ms = str(value).split(', ')[1][:-1]
            meshSize = ms.split(': ')[1]
            nX.append( int(meshSize)) 
        return (nX)

#Form a specific file structure for the case. It makes parsing easy.
    @classmethod
    def fileStructure(self,dataFile, pattern):
        #File structure    
        dataFolder = "/postProcessing/"
        dataFile = dataFile
        cwd = os.getcwd()
        casefolders = [cwd + "/" + folder for folder in os.listdir(cwd) if pattern in folder]    
        datafolders = [df+dataFolder for df in casefolders]
        datafolders.sort()
        fileNames = [fN+dataFile for fN in datafolders]
        return(fileNames)
