import os, sys


import numpy as np

import argparse
import textwrap
from argparse import RawTextHelpFormatter

DOC = """
---------------------------

Plot the statistical analysis of a bundle

Command example:
----------------
python3 showColorCodedBundlesByAge.py \
-tsv subjectsForTractogram.tsv \
-v 1 \
"""
####################################################################################

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
        prog="python3 computeMeanAndStdFeaturesAtlas.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-a", "--atlas-dir",
        type=is_dir, required=True, metavar="<path>",
        help=( "Path to the atlaas directory WITH .minf files" ) )


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

def readMinf( path ) :
    out_dict = {}
    with open( path, "r" ) as f :
        for line in f :
            words = line.split( " : " )
            words = [ tmp.replace( " ", "" ).replace( "'", "" ).replace( ",", "" ) for tmp in words ]
            if words[ 0 ] == "averageRadius" :
                out_dict[ words[ 0 ] ] = float( words[ 1 ] )
            if words[ 0 ] == "averageAngle" :
                out_dict[ words[ 0 ] ] = float( words[ 1 ] )
            if words[ 0 ] == "averageDirectionAngle" :
                out_dict[ words[ 0 ] ] = float( words[ 1 ] )
            if words[ 0 ] == "averageShapeAngle" :
                out_dict[ words[ 0 ] ] = float( words[ 1 ] )
            if words[ 0 ] == "averageLength" :
                out_dict[ words[ 0 ] ] = float( words[ 1 ] )
            if words[ 0 ] == "averageDistanceBetweenMedialPoints" :
                out_dict[ words[ 0 ] ] = float( words[ 1 ] )
    
    return( out_dict )



"""
attributes = {
    'averageRadius' : 8.50211482952917,
    'minRadius' : 2.6980835231151055,
    'maxRadius' : 16.550452281003675,
    'averageAngle' : 44.94805129018535,
    'minAngle' : 6.482599559069678,
    'maxAngle' : 83.58257100176971,
    'averageDirectionAngle' : 48.31607149845826,
    'minDirectionAngle' : 7.911592912053183,
    'maxDirectionAngle' : 117.13312562553148,
    'averageShapeAngle' : 144.18866689570052,
    'minShapeAngle' : 116.30650579361702,
    'maxShapeAngle' : 166.9266331677386,
    'averageLength' : 102.21255788039441,
    'minLength' : 75.46484357567138,
    'maxLength' : 129.38835202152492,
    'averageDisimilarity' : 9.162441055967609,
    'minDisimilarity' : 0,
    'maxDisimilarity' : 15.816149502903867,
    'density' : 0.00737555,
    'centerBundleX' : -36.8488,
    'centerBundleY' : -23.7169,
    'centerBundleZ' : 33.1087,
    'averageDistanceBetweenMedialPoints' : 11.564442778632502,
    'minDistanceBetweenMedialPoints' : 3.3664218254902343,
    'maxDistanceBetweenMedialPoints' : 23.68103662361494,
    'disimilarityWithAtlas' : -1,
    'coverageWithAtlas' : -1,
    'overlapWithAtlas' : -1,
    'adjacencyWithAtlas' : -1
,
}
"""
####################################################################
def main() :
    inputs, verbose = get_cmd_line_args()

    # atlas_dir_path = "/volatile/Atlases/Pamela_tck_transformed_15p/"
    atlas_dir_path = inputs[ "atlas_dir" ]

    #-------------------------------------------------------------#
    minf_list = []
    for tmpFile in os.listdir( atlas_dir_path ) :
        if tmpFile.endswith( ".minf" ) :
            minf_list.append( tmpFile )

    attributes_names = [ "averageRadius",
                        "averageAngle",
                        "averageDirectionAngle",
                        "averageShapeAngle",
                        "averageLength",
                        "averageDistanceBetweenMedialPoints" ]

    attributesValuesDict = {}
    means_dict = {}
    std_dict = {}
    for tmpAttribute in attributes_names :
        attributesValuesDict[ tmpAttribute ] = []
        means_dict[ tmpAttribute ] = -1
        std_dict[ tmpAttribute ] = -1

    nbFiles = len( minf_list )
    counter = 0
    for tmpFile in minf_list :
        tmpPath = os.path.join( atlas_dir_path, tmpFile )
        tmpDict = readMinf( tmpPath )
        if len( list( tmpDict.keys() ) ) != 6 :
            print( f"ERROR for file {tmpFile}: The dict does not contain 6 attributes" )
            sys.exit( 1 )
        for tmpAttribute in attributes_names :
            attributesValuesDict[ tmpAttribute ].append( tmpDict[ tmpAttribute ] )

    for tmpAttribute in attributes_names :
        means_dict[ tmpAttribute ] = np.mean( attributesValuesDict[ tmpAttribute ] )
        std_dict[ tmpAttribute ] = np.std( attributesValuesDict[ tmpAttribute ] )

    print( f"Means : \n{means_dict}" )
    print( f"Std : \n{std_dict}" )

if __name__ == "__main__" :
    main()
