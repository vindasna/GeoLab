#!${PYTHON_BINARY}

import sys, os, subprocess

import time

import numpy as np

import nibabel as nib
from nibabel.streamlines import Field
from nibabel.orientations import aff2axcodes

import argparse
import textwrap
from argparse import RawTextHelpFormatter

import multiprocessing

import json


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

Plot the statistical analysis of a bundle

Command example:
----------------
python3 createHCPTractograms.py \
-i input.ana \
-v 1 \

"""

def is_file( filepath ):
    """ Check file's existence - argparse 'file' argument.
    """
    if not os.path.isfile( filepath ):
        raise argparse.ArgumentError( "File does not exist: %s" % filepath )
    return filepath

def is_dir( filepath ):
    """ Check file's existence - argparse 'directory' argument.
    """
    if not os.path.isdir( filepath ):
        raise argparse.ArgumentError( "File does not exist: %s" % filepath )
    return filepath


def get_cmd_line_args():
    """
    Create a command line argument parser and return a dict mapping
    <argument name> -> <argument value>.
    """
    parser = argparse.ArgumentParser(
        prog="python3 createHCPTractograms.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input",
        type=is_dir, required=True, metavar="<path>",
        help="Directory with with recognized bundles")
    required.add_argument(
        "-a", "--atlas",
        type=is_dir, required=True, metavar="<path>",
        help="Path to the directory with the atlas" )

    # Optional arguments
    parser.add_argument(
        "-o", "--output",
        type=str, metavar="<path>",
        help="Output filename (.txt) ")
    parser.add_argument(
        "-minNbFibers", "--minNbFibers",
        type=int, metavar="<int>", default=1,
        help=( "Minimum number of fibers to count a bundle as recognized "
                                                               "(default : 1)"))
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
        type=int, choices=[0, 1, 2], default=0,
        help="Increase the verbosity level: 0 silent, 1 verbose.")


    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)
    verbose = kwargs.pop("verbose")

    return kwargs, verbose


def run_sh_process( cmd, shell = False ) :
    command_name = cmd[ 0 ]
    if shell :
        cmd = " ".join( cmd )
    process = subprocess.Popen( cmd, stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE, shell = shell )
    output_cmd, error_cmd = process.communicate()
    returncode_cmd = process.returncode
    if returncode_cmd != 0 :
        print( f"\n{bcolors.FAIL}ERROR with command {command_name}"
               f"{bcolors.ENDC}" )
        output_cmd = output_cmd.decode( "utf-8" )
        error_cmd = error_cmd.decode( "utf-8" )
        print( f"{bcolors.OKBLUE}Output:{bcolors.ENDC} \n{output_cmd}"
               f"{bcolors.OKBLUE}Error:{bcolors.ENDC} \n{error_cmd}" )
        sys.exit( 1 )
        return( 0 )
    return( output_cmd )


def read_bundle( bundle_filename ) :
    # Checking extension bundle_filename
    if ( ( not bundle_filename.endswith( ".bundlesdata" ) ) and
                              ( not bundle_filename.endswith( ".bundles" ) ) ) :
        print( "ERROR : Wrong format for the input of readBundlesFile, got "
              f"{bundle_filename} which is not .bundles/bundlesdata " )
        sys.exit()
    else:
        if bundle_filename.endswith( ".bundlesdata" ):
            bundle_filename = bundle_filename.replace( ".bundlesdata",
                                                                    ".bundles" )

    attributes = [ 'binary',
                   'bundles',
                   'byte_order',
                   'curves_count',
                   'radio',
                   'length',
                   'nSubjects',
                   'data_file_name',
                   'format',
                   'io_mode',
                   'item_count',
                   'label_type',
                   'labels',
                   'object_type',
                   'averageRadius',
                   'minRadius',
                   'maxRadius',
                   'averageAngle',
                   'minAngle',
                   'maxAngle',
                   'averageDirectionAngle',
                   'minDirectionAngle',
                   'maxDirectionAngle',
                   'averageShapeAngle',
                   'minShapeAngle',
                   'maxShapeAngle',
                   'averageLength',
                   'minLength',
                   'maxLength',
                   'averageDisimilarity',
                   'minDisimilarity',
                   'maxDisimilarity',
                   'density',
                   'centerBundleX',
                   'centerBundleY',
                   'centerBundleZ',
                   'resolutionX',
                   'resolutionY',
                   'resolutionZ',
                   'sizeX',
                   'sizeY',
                   'sizeZ',
                   'space_dimension' ]

    minf_dict = dict()
    for i in range( len( attributes ) ) :
        minf_dict[ attributes[ i ] ] = -1


    with open( bundle_filename, "r" ) as f:
        for line in f:
            words = line.split( ":" )
            attribute_name = words[ 0 ].replace( "'", "").replace( " ", "" )
            if attribute_name in minf_dict :
                minf_dict[ attribute_name ] = words[ 1 ].replace( ",\n", "" )

    return minf_dict


def main() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()
    #################################################
    input_dir = inputs[ "input" ]
    #################################################
    atlas_dir = inputs[ "atlas" ]
    #################################################
    minNbFibers = inputs[ "minNbFibers" ]
    #################################################
    force = inputs[ "overwrite" ]
    #################################################
    parallel = inputs[ "multiprocessing" ]
    #################################################
    output_file = inputs[ "output" ]
    if output_file:
        if os.path.isfile( output_file ) and not force :
            print( "ERROR : Output file already exists, use the -force flag "
                   "if you want to overwrite the file" )
            sys.exit( 1 )
    else :
        output_file = os.path.join( input_dir, "countFoundBundles.txt" )
        if os.path.isfile( output_file ) and not force :
            print( "ERROR : Output file already exists, use the -force flag "
                   "if you want to overwrite the file" )
            sys.exit( 1 )

    ## ---------------------------------------------------------------------- ##
    fibers_per_recognized_bundle = {}
    for bundle in os.listdir( atlas_dir ) :
        if bundle.endswith( ".bundles" ) :
            bundle_name = bundle.replace( ".bundles", "" )
            fibers_per_recognized_bundle[ bundle_name ] = 0

    nbRecognizedBundles = 0
    nbFibersPerBundles = []
    for bundle_name in fibers_per_recognized_bundle :
        bundle_filename = os.path.join( input_dir, f"{bundle_name}.bundles" )
        if not (os.path.isfile( bundle_filename ) ) :
            continue
        bundleInfo = read_bundle( bundle_filename )
        nbFibersInBundle = int( bundleInfo[ "curves_count" ] )
        if ( nbFibersInBundle > minNbFibers ) :
            nbRecognizedBundles += 1
            fibers_per_recognized_bundle[ bundle_name ] = nbFibersInBundle
            nbFibersPerBundles.append( nbFibersInBundle )


    nbFibersPerBundlesMean = np.mean( nbFibersPerBundles )
    nbFibersPerBundlesMedian = np.median( nbFibersPerBundles )
    nbFibersPerBundlesStd = np.std( nbFibersPerBundles )

    percentageRecognized = ( nbRecognizedBundles /
                                    len( fibers_per_recognized_bundle.keys() ) )

    print( f"Scores per bundles ( mean (median) +- std ) :\n"
           f"Number of fibers per bundle : {nbFibersPerBundlesMean}"
           f"({nbFibersPerBundlesMedian}) +- {nbFibersPerBundlesStd} \n"
           f"Percentage recognized : {percentageRecognized}" )

    file = open( output_file, 'w' )
    for bundle_name in fibers_per_recognized_bundle :
        file.write( f"{bundle_name} : "
                             f"{fibers_per_recognized_bundle[bundle_name]} \n" )


if __name__ == "__main__" :
    main()
