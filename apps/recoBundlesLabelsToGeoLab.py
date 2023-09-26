#!${PYTHON_BINARY}
import sys, os, subprocess, shutil

import time

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
measures_options = [ "FA", "MD", "OD", "ICVF", "ISOVF" ]
slope_dict = {}
verbose = 0
# To have reproductible results
_seed = 100
np.random.seed( _seed )
standartOut = sys.stdout
standartErr = sys.stderr


DOC = """
---------------------------

Plot the statistical analysis of a bundle

Command example:
----------------
python3 recoBundlesLabelsToGeoLab.py \
-rd RecoBundlesDir \
-ln GeoLabLabelsNames.dict \
-o outLabels.txt \
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
        prog="python3 recoBundlesLabelsToGeoLab.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-rd", "--recobundles-dir",
        type=is_dir, required=True, metavar="<path>",
        help=( "Output directory of RecoBundles with the .npy files" ) )

    required.add_argument(
        "-ln", "--labels-names",
        type=is_file, required=True, metavar="<path>",
        help=( "Dictionary of labels names (.dict in GeoLab)" ) )

    required.add_argument(
        "-o", "--out",
        type=str, required=True, metavar="<path>",
        help=( "Output labels in .txt format" ) )

    # Optional arguments
    parser.add_argument(
        "-tn", "--template-name",
        type=str, default = "moved_{}_recognized",
        help=("Template for .trk bundles names in the form of \"prefix_{}_sufix\""
              " where {} is replaced by the bundle name (default : "
              "moved_{}_recognized). ") )

    parser.add_argument(
        "-v", "--verbose",
        type=int, choices=[0, 1, 2], default=0,
        help="Increase the verbosity level: 0 silent, 1 verbose." )


    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)
    verbose = kwargs.pop("verbose")

    return kwargs, verbose



def readDict( path ) :
    _dict = {}
    with open( path, "r" ) as f:
        for line in f:
            words = line.split( ":" )
            _dict[ int( words[ 1 ] ) ] = words[ 0 ].replace( " ", "" )

    return( _dict )

def getLabelFromBundleName( bundlesDict, bundleName ) :
    for label in bundlesDict.keys() :
        if bundlesDict[ label ] == bundleName :
            return( label )
    print( f"Error in getLabelFromBundleName() : \'{bundleName}\' not found in "
           f"bundlesDict" )
    print( f"Bundles names in bundlesDict : " )
    for label in bundlesDict.keys() :
        print( f"\'{bundlesDict[ label ]}\'" )
    sys.exit( 1 )


def saveLabels( labels, path ) :
    # Function where labels is a dictionary and NOT a list
    nbFibers = max( list( labels.keys() ) ) + 1
    with open( path, "w" ) as f :
        for fiber in range( nbFibers ) :
            if fiber in labels.keys() :
                for _label in labels[ fiber ] :
                    f.write( f"{fiber} : {_label}\n" )
            else :
                f.write( f"{fiber} : -1\n" )

def saveDict( inDict, path ) :
    with open( path, 'w' ) as f :
        for _key in inDict :
            f.write( f"{inDict[ _key ]} : {_key}\n" )

def main() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    recoBundlesDir = inputs[ "recobundles_dir" ]

    labelsNamesPath = inputs[ "labels_names" ]

    outputLabelsFile= inputs[ "out" ]
    if not outputLabelsFile.endswith( ".txt" ) :
        print( "ERROR : output file must be .txt" )
        sys.exit( 1 )

    template_bundle_name = inputs[ "template_name" ]

    tmpFiles = os.listdir( recoBundlesDir )
    recoLabelsPaths = []
    for _file in tmpFiles :
        if _file.endswith( ".npy" ) :
            recoLabelsPaths.append( os.path.join( recoBundlesDir, _file ) )
    tmpFiles = None # Free memory


    labelNames = readDict( labelsNamesPath )

    nbBundles = len( recoLabelsPaths )

    outRecoLabels = {} # Use dict instead of list for simplicity
    counter = 1
    tmpPrefix = template_bundle_name.split( "_{}_" )[ 0 ]
    tmpSufix = template_bundle_name.split( "_{}_" )[ 1 ]
    max_value_recobundles_dict = int( np.max( list( labelNames.keys() ) ) )
    for labelsPath in recoLabelsPaths :
        print( f"Processing : [{counter}/{nbBundles}]", end = "\r" )
        tmp = np.load( labelsPath )
        _bundleName = os.path.basename( labelsPath ).replace( "__labels.npy",
                                                                            "" )
        
        _bundleName = _bundleName.replace( f"{tmpPrefix}_", "" )
        _bundleName = _bundleName.replace( f"_{tmpSufix}", "" )

        if _bundleName not in labelNames.keys() :
            max_value_recobundles_dict += 1
            labelNames[ max_value_recobundles_dict ] = _bundleName

        _labelValue = getLabelFromBundleName( labelNames, _bundleName )
        for _fiber in tmp :
            if _fiber not in outRecoLabels.keys() :
                outRecoLabels[ _fiber ] = [ _labelValue ]
            else :
                outRecoLabels[ _fiber ].append( _labelValue )
        counter += 1

    print( "\nDone" )

    saveLabels( outRecoLabels, outputLabelsFile )

    outputDictFile = outputLabelsFile.replace( ".txt", ".dict" )

    saveDict( labelNames, outputDictFile )

    return



if __name__ == "__main__" :
    t1 = time.time()
    main()
    duration = time.time() - t1
    print( f"Duration : {duration} s" )
