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

#Form a specific file structure for the case. It makes parsing easy.
    @classmethod
    def fileStructure(self,caseFolder,dataFile, pattern):
        #File structure    
        dataFolder = "/postProcessing/"
        dataFile = dataFile
        cwd = os.getcwd()+'/'+caseFolder
        casefolders = [cwd + "/" + folder for folder in os.listdir(cwd) if pattern in folder]    
        datafolders = [df+dataFolder for df in casefolders]
        datafolders.sort()
        fileNames = [fN+dataFile for fN in datafolders]
        return(fileNames)
    
    @classmethod
    #Function gives dictionary mapping variation number to parameter vector
    def fileMapping(self,label):
        cwd = os.getcwd()
        casefolders = [folder for folder in os.listdir(cwd) if label in folder] 
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
            ms = str(value).split(', ')[2][:-1]
            meshSize = ms.split(': ')[-1]
            nX.append( int(0.001/float(meshSize))) 
        return (nX)


