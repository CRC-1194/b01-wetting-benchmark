import numpy as np
import math
import ast
import vtk
import os
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN
import csv

pi = np.pi
sin = np.sin
cos = np.cos


class height_vector:
    x_coord = 0.0
    y_coord = 0.0
    z_coord = 0.0


class Funcs:
    #Formulas taken from Mathis et al. preprint
    @classmethod
    def getAnalyticalWettedArea(self, theta, dropletRadius):
        dropletVolume = 4 / (2*3) * pi * dropletRadius**(3)
        g_theta = sin(theta) * ((pi* ((1-cos(theta))**2) * (2+cos(theta))) / 3)**(-1/3)
        wettedRadius = dropletVolume**(1/3) * g_theta
        wettedArea = pi * wettedRadius**(2)
        r = wettedRadius / sin(theta)
        sphericalCapHeight = r * (1-cos(theta))
        return wettedArea*1000000

    @classmethod
    def getAnalyticalWettedRadius(self,theta,dropletRadius):
        dropletVolume = 4 / (2*3) * pi * dropletRadius**(3)
        g_theta = sin(theta) * ((pi* ((1-cos(theta))**2) * (2+cos(theta))) / 3)**(-1/3)
        wettedRadius = dropletVolume**(1/3) * g_theta
        wettedArea = pi * wettedRadius**(2)
        r = wettedRadius / sin(theta)
        sphericalCapHeight = r * (1-cos(theta))
        return wettedRadius, sphericalCapHeight

    @classmethod
    def getAnalyticalWettedRadiusByVolume(self,theta,dropletVolume):
        g_theta = sin(theta) * ((pi* ((1-cos(theta))**2) * (2+cos(theta))) / 3)**(-1/3)
        wettedRadius = dropletVolume**(1/3) * g_theta
        wettedArea = pi * wettedRadius**(2)
        r = wettedRadius / sin(theta)
        sphericalCapHeight = r * (1-cos(theta))
        return wettedRadius, sphericalCapHeight

    @classmethod
    def getR(self,vofFile):
        dropRadius = 0.0
        var_lines = open (vofFile).readlines()
        for line in var_lines:
            if 'radius' in line:
                dropRadius = float(line.split()[1][:-1])
        return dropRadius
    
#Form a specific file structure for the case. It makes parsing easy.
    @classmethod
    def fileStructure(self,dataFolder, dataFile, pattern, var_list):
        fileNames = []
        #File structure    
        dataFile = dataFile
        cwd = os.getcwd()
        for idx, id in enumerate(var_list):
            casefolders = [cwd + "/" + folder for folder in os.listdir(cwd) if pattern in folder and var_list[idx] in folder]    
            datafolders = [df+dataFolder for df in casefolders] 
            datafolders.sort()
            files = [fN+dataFile for fN in datafolders]
            fileNames.append(files[0])
        return(fileNames)

    @classmethod
    def fileMapping(self,label, sm):
        cwd = os.getcwd()
        casefolders = [folder for folder in os.listdir(cwd) if label in folder] 
        casefolders.sort()
        mapNumber = [int(num.split('_')[1]) for num in  casefolders ]
        var_map = {}
        var_lines = open("variation_file").readlines()
        for line in var_lines:
                # Skip lines without mapping
            if '{' not in line:
                continue
            if sm in line:
                var_num = int(line.split()[1])
                if(var_num in mapNumber):
                    dict_start = line.find('{')
                    # Mappings in variation file can directly be interpreted by Python
                    var_map[var_num] = ast.literal_eval(line[dict_start:-1])
        return (var_map)

#reads the surface .vtk file
#calculates the interface height
#writes the height in .csv file
    @classmethod
    def writeHeightFile(self,calcHeights, wettedAreaFileNames, vtk_folders, heightFileNames):
        if(calcHeights):
            for idx, vtk_folder in enumerate(vtk_folders):
                #getting the time step folder names
                sub_folders = [name for name in os.listdir(vtk_folder) if os.path.isdir(os.path.join(vtk_folder, name))]
                sub_folders = sorted(sub_folders, key=float)

                #file-content
                # header = ['Time(s)', 'Height(m)']
                with open(heightFileNames[idx], 'w') as f:
                    writer = csv.writer(f)

                for sub_folder in sub_folders:
                    heightObject = height_vector()
                    heightObject.z_coord = 0
                    #ifCoordArray = 0
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
                                #print(coordBreakUp)
                                coordPerLine = int(len(coordBreakUp) / 3)
                                
                                if (counter ==0):
                                    heightObject.z_coord = float (coordBreakUp[2])
                                    counter =1
                                if(counter==1):
                                #check on Z coordinate
                                    for i in range(0, coordPerLine):
                                        coordinateIndex = ((int(len(coordBreakUp)-1)%3)) + (i*3)
                                        heightObject.z_coord = max (heightObject.z_coord, float(coordBreakUp[coordinateIndex]))

                            if "\n" in line[0] and ifCoordArray ==1:
                                ifCoordArray =2
    
                    data_to_write_to_csv = [sub_folder, str(heightObject.z_coord)]
                    
                    with open(heightFileNames[idx], 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow(data_to_write_to_csv)
