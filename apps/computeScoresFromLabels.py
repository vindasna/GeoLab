#!/usr/bin/python3

import sys, os, subprocess, shutil

import numpy as np
import scipy.stats
import scipy.special
import matplotlib.pyplot as plt

import argparse
import textwrap
from argparse import RawTextHelpFormatter

import multiprocessing



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

Compute scores for a segmentation

Command example:
----------------
python3 computeScoresFromLabels.py \
-pl predicted.labels \
-pd predicted.dict \
-tl predicted.labels \
-td predicted.dict \
-v 1 \

"""

def is_file_or_directory(filepath):
    """ Check file's existence - argparse 'type' argument.
    """
    if not os.path.isfile(filepath):
        if not os.path.isdir(filepath):
            raise argparse.ArgumentError( "File or directory does not exist: %s"
                                                                     % filepath)
    return filepath

def is_file(filepath):
    """ Check file's existence - argparse 'type' argument.
    """
    if not os.path.isfile(filepath):
        raise argparse.ArgumentError( "File does not exist: %s" % filepath )

    return filepath


def get_cmd_line_args():
    """
    Create a command line argument parser and return a dict mapping
    <argument name> -> <argument value>.
    """
    parser = argparse.ArgumentParser(
        prog="python3 analyseAtlasBundle.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group( "required arguments" )
    required.add_argument(
        "-pl", "--predicted-labels",
        type=is_file, required=True, metavar="<path>",
        help="Path to the predicted labels .txt" )

    required.add_argument(
        "-pd", "--predicted-dictionary",
        type=is_file, required=True, metavar="<path>",
        help="Path to the predicted labels dictionary .dict" )

    required.add_argument(
        "-tl", "--true-labels",
        type=is_file, required=True, metavar="<path>",
        help="Path to the true labels .txt" )

    required.add_argument(
        "-td", "--true-dictionary",
        type=is_file, required=True, metavar="<path>",
        help="Path to the true labels dictionary .dict" )

    required.add_argument(
        "-o", "--output",
        type=str, metavar="<path>",
        help="Output file .tsv " )

    # Optional arguments
    parser.add_argument(
        "-force", "--overwrite",
        action='store_true', default=False,
        help="Overwrite output file even if it already exists" )
    parser.add_argument(
        "-parallel", "--multiprocessing",
        action='store_true', default=False,
        help="Use multiprocessing for computations" )
    parser.add_argument(
        "-v", "--verbose",
        type=int, choices=[0, 1, 2], default=1,
        help="Increase the verbosity level: 0 silent, 1 verbose." )


    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)
    verbose = kwargs.pop("verbose")

    return kwargs, verbose


def readLabels( pathLabels, pathDict ) :
    if not ( pathLabels.endswith( ".txt" ) and pathDict.endswith( ".dict" ) ) :
        print( f"ERROR : in readLabels, pathLabels must be .txt and pathDict "
               f"must be .dict" )
        sys.exit( 1 )



    labelsDict = {}
    with open( pathDict, "r" ) as f :
        for line in f :
            words = line.split( ":" )
            bundleName = words[ 0 ].replace( " ", "" )
            labelBundle = int( words[ 1 ] )
            labelsDict[ labelBundle ] = bundleName

    labelsOut = {}
    with open( pathLabels, "r" ) as f :
        for line in f :
            words = line.split( ":" )
            fiberIndexBundle = int( words[ 0 ].replace( " ", "" ) )
            labelBundle = int( words[ 1 ] )
            if labelBundle > len( labelsDict.keys() ) :
                print( f"ERROR : in readLabels, len( labelsDict.keys() ) = "
                       f"{len( labelsDict.keys() )} smaller than labelBundle = "
                       f"{labelBundle})"  )
                sys.exit( 1 )
            if ( fiberIndexBundle not in labelsOut.keys() ) :
                labelsOut[ fiberIndexBundle ] = [ labelsDict[ labelBundle ] ]
            else :
                labelsOut[ fiberIndexBundle ].append(
                                                     labelsDict[ labelBundle ] )

    return labelsOut




def main():
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()
    predictedLabelsPath = inputs[ "predicted_labels" ]
    predictedDictPath = inputs[ "predicted_dictionary" ]
    trueLabelsPath = inputs[ "true_labels" ]
    trueDictPath = inputs[ "true_dictionary" ]
    outputPath = inputs[ "output" ]
    force = inputs[ "overwrite" ]
    useMultiprocessing = inputs[ "multiprocessing" ]

    if ( not force and os.path.isfile( outputPath ) ):
        print( f"ERROR : output file {outputPath} already exists and -force "
               f"was not used" )
        sys.exit( 1 )


    predictedLabels = readLabels( predictedLabelsPath, predictedDictPath )
    trueLabels = readLabels( trueLabelsPath, trueDictPath )

    
