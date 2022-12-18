#!/usr/bin/python

from optparse import OptionParser
import sys
from subprocess import call
import os


def set_up_dimensions (dim_1, dim_2): 
    """Renames 0/field.dim_1 to 0/field and removes 0/field.dim_2"""
    files_dim_1 = [os.path.join("0",fname) for fname in os.listdir("0") \
                if os.path.isfile(os.path.join("0", fname)) and \
                fname.endswith(dim_1)] 
    for file_dim_1 in files_dim_1:
        base_name = file_dim_1.rstrip("." + dim_1)
        call(["mv", file_dim_1, base_name])
        call(["rm", "-rf", base_name + "." + dim_2])

usage = """A wrapper for pyFoamRunParameterVariation.py that generates the
directory structure for a parameter study. 

Meshes are not generated and preprocessing is not done. 
Used to prepare for execution on a cluster.

Usage: ./create-study.py -c templateCase -p paramFile -s studyName"""

study_name=""

parser = OptionParser(usage=usage)

parser.add_option("-c", "--case", dest="casedir",
                  help="Template case directory.", 
                  metavar="CASEDIR")

parser.add_option("-p", "--parameter-file", dest="paramfile", 
                  help="PyFoam parameter file used by pyFoamRunParameterVariation.py.", 
                  metavar="PARAMFILE")

parser.add_option("-s", "--study-name", dest="studyname", 
                  help="Name of the parameter study.", 
                  metavar="STUDYNAME")
if __name__ == "__main__":

    (options, args) = parser.parse_args()

    if ((options.casedir == None) or  
        (options.paramfile == None) or 
        (options.studyname == None)): 
        print ("Error: case or parameter option not used. Use --help option for more information.") 
        sys.exit(1)

    (options, args) = parser.parse_args()
    #numberOfVariations = [25, 31, 37, 43, 49, 50, 56, 62, 68, 74, 75, 81, 87, 93, 99]
    numberOfVariations = [0, 1, 2, 3, 4]
    study_name = options.studyname
    
    for nVar in numberOfVariations:
        # Generate parameter study simulation cases 
        call(["pyFoamRunParameterVariation.py", "--no-execute-solver", "--no-server-process", 
              "--no-mesh-create", "--no-case-setup", "--cloned-case-prefix=%s" % options.studyname, 
              "--every-variant-one-case-execution", "--single-variation=%i" %nVar, 
              "--create-database", options.casedir, options.paramfile])


os.system("echo "+ study_name + "* | xargs -n 1")          
