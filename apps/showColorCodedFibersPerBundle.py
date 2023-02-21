#!/usr/bin/python3
import os, sys, shutil
import numpy as np
from dipy.viz import window, actor
from dipy.data import fetch_bundles_2_subjects, read_bundles_2_subjects
from dipy.tracking.streamline import transform_streamlines
from dipy.io.stateful_tractogram import Space, StatefulTractogram
from dipy.io.streamline import load_tractogram, save_tractogram

import pickle

import matplotlib.pyplot as plt


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

Plot the statistical analysis of a bundle

Command example:
----------------
python3 showColorCodedFibersPerBundle.py \
-a atlasDir \
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
        prog="python3 showColorCodedFibersPerBundle.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-a", "--atlas-dir",
        type=is_dir, required=True, metavar="<path>",
        help=( "Path to the atlaas directory WITH .TRK files" ) )

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

def main() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    bundles_dir = inputs[ "atlas_dir" ]

    ############################################################################
    bundles_names = []
    for tmpBundleName in os.listdir( bundles_dir ) :
        if tmpBundleName.endswith( ".trk" ) :
            bundles_names.append( tmpBundleName.replace( ".trk", "" ) )


    totalNbStreamlines = 0
    for _bundle_name in bundles_names :
        bundle_path = os.path.join( bundles_dir, f"{_bundle_name}.trk" )
        bundle_data = load_tractogram( bundle_path , "same" )
        totalNbStreamlines += len( bundle_data )


    _streamlines = []
    _measure_values = []

    for _bundle_name in bundles_names :
        bundle_path = os.path.join( bundles_dir, f"{_bundle_name}.trk" )
        bundle_data = load_tractogram( bundle_path , "same" )
        nbStreamlines = len( bundle_data )
        tmpMeasureValue = nbStreamlines / totalNbStreamlines
        for fiber in bundle_data.streamlines :
            _streamlines.append( fiber )
        _measure_values += nbStreamlines * [ tmpMeasureValue ]

    measure_min = min( _measure_values )
    measure_max = max( _measure_values )

    scene = window.Scene()

    _streamlines = np.array( _streamlines )
    _measure_values = np.array( _measure_values )


    lut_cmap = actor.colormap_lookup_table( scale_range = ( measure_min,
                                                                 measure_max ) )
    stream_actor2 = actor.line( _streamlines, _measure_values, linewidth=0.1,
                                                      lookup_colormap=lut_cmap )

    bar = actor.scalar_bar( lut_cmap )

    scene.add(stream_actor2)
    scene.add(bar)
    window.show(scene, size=(600, 600), reset_camera=False)

if __name__ == "__main__" :
    main()
