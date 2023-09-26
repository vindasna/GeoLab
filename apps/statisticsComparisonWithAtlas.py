#!${PYTHON_BINARY}
import sys, os, shutil
import numpy as np

import argparse
import textwrap
from argparse import RawTextHelpFormatter

#---- For printing colors in terminal ----#
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#-----------------------------------------#


DOC = """
---------------------------

Plot the statistical analysis of a comparisonWithAtlas.tsv

Command example:
----------------
python3 statisticsComparisonWithAtlas.py \
-i input.tsv
-v 1 \

"""

def is_file(filepath):
    """ Check file's existence - argparse 'file' argument.
    """
    if not os.path.isfile(filepath):
        raise argparse.ArgumentError("File does not exist: %s" % filepath)
    return filepath

def is_dir(filepath):
    """ Check file's existence - argparse 'dir' argument.
    """
    if not os.path.isdir(filepath):
        raise argparse.ArgumentError("Directory does not exist: %s" % filepath)
    return filepath


def get_cmd_line_args():
    """
    Create a command line argument parser and return a dict mapping
    <argument name> -> <argument value>.
    """
    parser = argparse.ArgumentParser(
        prog="python3 statisticsComparisonWithAtlas.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input",
        type=is_file, required=True, metavar="<path>",
        help=( "Input file (.tsv)" ) )

    # Optional arguments
    parser.add_argument(
        "-v", "--verbose",
        type=int, choices=[0, 1, 2], default=0,
        help="Increase the verbosity level: 0 silent, 1 verbose.")


    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)
    verbose = kwargs.pop("verbose")

    return kwargs, verbose


def readComparisonWithAtlas( path ) :
    if not os.path.isfile( path ) :
        print( f"ERROR in readComparisonWithAtlas : file {path} does not exist")
        exit( 1 )

    outDict = {}
    with open( path, 'r' ) as f :
        for line in f :
            words = line.split( "\t" )
            if words[ 0 ] == "bundleName" or words[ 0 ] == "Bundle_Name":
                outDict[ "bundleName" ] = []
                outDict[ "Coverage" ] = []
                outDict[ "Adjacency" ] = []
                outDict[ "Overlap" ] = []
                outDict[ "NbFibers" ] = []
            else :
                outDict[ "bundleName" ].append( words[ 0 ] )
                outDict[ "Coverage" ].append( float( words[ 1 ] ) )
                outDict[ "Adjacency" ].append( float( words[ 2 ] ) )
                outDict[ "Overlap" ].append( float( words[ 3 ] ) )
                outDict[ "NbFibers" ].append( float( words[ 5 ] ) )

    return( outDict )

def main() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    input_tsv_path = inputs[ "input" ]

    dictComparison = readComparisonWithAtlas( input_tsv_path )

    coverageMean = np.mean( dictComparison[ "Coverage" ] )
    coverageMedian = np.median( dictComparison[ "Coverage" ] )
    coverageStd = np.std( dictComparison[ "Coverage" ] )
    adjacencyMean = np.mean( dictComparison[ "Adjacency" ] )
    adjacencyMedian = np.median( dictComparison[ "Adjacency" ] )
    adjacencyStd = np.std( dictComparison[ "Adjacency" ] )
    overlapMean = np.mean( dictComparison[ "Overlap" ] )
    overlapMedian = np.median( dictComparison[ "Overlap" ] )
    overlapStd = np.std( dictComparison[ "Overlap" ] )
    nbFibersMean = np.mean( dictComparison[ "NbFibers" ] )
    nbFibersMedian = np.median( dictComparison[ "NbFibers" ] )
    nbFibersStd = np.std( dictComparison[ "NbFibers" ] )

    # PBE-1 and PBE-10
    counter = 0
    pbe1 = 0
    pbe10 = 0
    for tmp in dictComparison[ "NbFibers" ] :
        if tmp >= 1 :
            pbe1 += 1
        if tmp >= 10 :
            pbe10 += 1
        counter += 1 
    pbe1 /= counter
    pbe10 /= counter


    print( f"Scores per bundles ( mean (median) +- std ) :\n"
           f"Coverage : {coverageMean}({coverageMedian}) +- {coverageStd} \n"
           f"Adjacency : {adjacencyMean}({adjacencyMedian}) +- {adjacencyStd}\n"
           f"Overlap : {overlapMean}({overlapMedian}) +- {overlapStd}\n" 
           f"NbFibers : {nbFibersMean}({nbFibersMedian}) +- {nbFibersStd}\n" 
           f"PBE-1 : {pbe1}\n" 
           f"PBE-10 : {pbe10}" )

if __name__ == "__main__" :
    main()
