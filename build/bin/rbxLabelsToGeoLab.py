#!/usr/bin/env python3

import os, sys, shutil

import json

import subprocess


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

Convert RecoBundlesX labels to GeoLab type labels

Command example:
----------------
python3 scorePredictionSGT.py \
-i rbxLabels.json \
-ol out.labels \
-od out.dict \

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

def is_dir(dirpath):
    """ Check file's existence - argparse 'type' argument.
    """
    if not os.path.isdir(dirpath):
        raise argparse.ArgumentError( "Directory does not exist: %s" % dirpath )

    return dirpath


def get_cmd_line_args():
    """
    Create a command line argument parser and return a dict mapping
    <argument name> -> <argument value>.
    """
    parser = argparse.ArgumentParser(
        prog="python3 rbsLabelsToGeoLab.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group( "required arguments" )
    required.add_argument(
        "-i", "--input",
        type=is_file, required=True, metavar="<path>",
        help="Path to the predicted labels (.json)" )
    
    required.add_argument(
        "-t", "--tractogram",
        type=is_file, required=True, metavar="<path>",
        help="Input tractogram given to RecoBundleX during labeling (.tck)" )

    required.add_argument(
        "-ol", "--out-labels",
        type=str, required=True, metavar="<path>",
        help="Path where to save labels .txt" )

    required.add_argument(
        "-od", "--out-dictionary",
        type=str, required=True, metavar="<path>",
        help=("Path where to dictionary .dict ") )



    # Optional arguments
    parser.add_argument(
        "-force", "--overwrite",
        action='store_true', default=False,
        help="Overwrite output file even if it already exists" )
    parser.add_argument(
        "-v", "--verbose",
        type=int, choices=[0, 1, 2], default=1,
        help="Increase the verbosity level: 0 silent, 1 verbose." )


    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)
    verbose = kwargs.pop("verbose")

    return kwargs, verbose


def saveLabels( in_list, path ) :
    nbStreamlines = len( in_list )
    with open( path, 'w' ) as f :
        for tmpFiberIndex in range( nbStreamlines ) :
            tmpCounter = 0
            for tmpLabel in in_list[ tmpFiberIndex ] :
                f.write( f"{tmpFiberIndex} : {tmpLabel}" )
                if tmpFiberIndex != ( nbStreamlines - 1 ) or tmpCounter != ( len( in_list[ tmpFiberIndex ] ) - 1 ) :
                    f.write( "\n" )


def saveDict( inDict, path ) :
    with open( path, 'w' ) as f :
        for _key in inDict :
            f.write( f"{inDict[ _key ]} : {_key}\n" )


def readRecoBundlesX( path, nbStreamlines ) :
    out_labels_dict = {}
    out_labels = {}
    for tmpFiber in range( nbStreamlines ) :
        out_labels[ tmpFiber ] = [ -1 ]

    with open( path, 'r' ) as f :
        tmpDict = json.load( f )
    
    counter = 0
    for tmpFile in tmpDict.keys() :
        tmpBundleName = os.path.basename( tmpFile )
        tmpBundleName = tmpBundleName.replace( ".trk", "" )
        tmpLabel = counter
        if tmpLabel not in out_labels_dict.keys() :
            out_labels_dict[ counter ] = tmpBundleName
            counter += 1
        
        for tmpIndex in tmpDict[ tmpFile ][ "indices" ] :
            if len( out_labels[ tmpIndex ] ) == 1 and out_labels[ tmpIndex ][ 0 ] == -1 :
                out_labels[ tmpIndex ] = [ tmpLabel ]
            else :
                out_labels[ tmpIndex ].append( tmpLabel )
    

   
    
    
    return( out_labels_dict, out_labels )

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
    return( output_cmd.decode( "utf-8" ) )


def getNbFibersBundles( input_file ) :
    tmpCommand = [ f"tckinfo {input_file} "
                   f"| grep \" count\"" ]
    outCmd = run_sh_process( tmpCommand, shell = True )
    nbFibers = int( outCmd.split( " " )[ -1 ] )
    return( nbFibers )



def main() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    input_rbx_path = inputs[ "input" ]

    input_tractogram_path = inputs[ "tractogram" ]
    if not input_tractogram_path.endswith( ".tck" ) :
        print( f"ERROR : input tractogram must be in .tck format" )
        sys.exit( 1 )

    out_labels_path = inputs[ "out_labels" ]

    out_dict_path = inputs[ "out_dictionary" ]

    force = inputs[ "overwrite" ]


    ###########################################################################################

    if os.path.isfile( out_labels_path ) and os.path.isfile( out_dict_path ) and not force :
        print( "Output files already exists and the -force flag was not use, skipping"
               "computations" )
        return

    nbFibers = getNbFibersBundles( input_tractogram_path )

    out_labels_dict, out_labels = readRecoBundlesX( input_rbx_path, nbFibers )

    saveDict( out_labels_dict, out_dict_path )

    saveLabels( out_labels, out_labels_path )











if __name__ == "__main__" :
    main()
