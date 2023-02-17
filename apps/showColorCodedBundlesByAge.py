#!/usr/bin/python3
import os, sys, shutil
import numpy as np
from dipy.viz import window, actor
from dipy.data import fetch_bundles_2_subjects, read_bundles_2_subjects
from dipy.tracking.streamline import transform_streamlines
from dipy.io.stateful_tractogram import Space, StatefulTractogram
from dipy.io.streamline import load_tractogram, save_tractogram

import pickle


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
python3 showColorCodedBundlesByAge.py \
-tsv subjectsForTractogram.tsv \
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
        prog="python3 showColorCodedBundlesByAge.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-a", "--atlas-dir",
        type=is_dir, required=True, metavar="<path>",
        help=( "Path to the atlaas directory WITH .TRK files" ) )

    required.add_argument(
        "-m", "--measure-dict",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the .pickle file containing the dictionary with "
               "the values of the measure slopes per bundle " ) )

    # required.add_argument(
    #     "-ma", "--measure-age",
    #     type=is_file, required=True, metavar="<path>",
    #     help=( "Path to the .txt containing the measure with age "
    #                                 "(usually called ${Measure}AndAge.txt) " ) )

    required.add_argument(
        "-age", "--age",
        type=int, required=True, metavar="<path>",
        help=( "Age to show color coded" ) )

    # Optional arguments
    parser.add_argument(
        "-wn", "--window-title",
        type=str, default="showColorCodedBundlesByAge",
        help="Title for the window (default : showColorCodedBundlesByAge )")

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


def readFileWithMeansMeasure( path ) :
    data_dict = {}
    i = 0
    with open( path, "r" ) as f :
        for line in f :
            if line.endswith( ":\n" ) :
                bundleName = line.replace( ":\n", "" )
                if bundleName not in data_dict.keys() :
                    data_dict[ bundleName ] = [ [], [], [], 0 ]
            else:
                if line != "" and line != "\n\n" and line != "\n" :
                    infoSubject = line.split( "\t" )
                    try :
                        data_dict[ bundleName ][ 0 ].append( int(
                                                            infoSubject[ 0 ] ) )
                        data_dict[ bundleName ][ 1 ].append( float(
                                                            infoSubject[ 1 ] ) )
                        data_dict[ bundleName ][ 2 ].append(
                                          infoSubject[ 2 ].replace( "\n", "" ) )
                        data_dict[ bundleName ][ 3 ] += 1
                    except :
                        print( f"|{line}|" )
                        sys.exit()
            i += 1

    return( data_dict )



def main() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    bundles_dir = inputs[ "atlas_dir" ]

    dict_values_path = inputs[ "measure_dict" ]

    # measure_with_age_path = inputs[ "measure_age" ]

    age = inputs[ "age" ]

    window_title = inputs[ "window_title" ]


    ############################################################################
    # measure_with_age_dict = readFileWithMeansMeasure( measure_with_age_path )

    # measure_population_mean_per_bundle = {}
    # for bundleName in measure_with_age_dict.keys() :
    #     measure_population_mean_per_bundle[ bundleName ] = np.mean(
    #                                   measure_with_age_dict[ bundleName ][ 1 ] )

    bundles_names = os.listdir( bundles_dir )
    with open(dict_values_path, 'rb') as handle:
        dict_values = pickle.load( handle )

    # measure_values_1 = []
    # measure_values_2 = []
    # breakpoints = []
    # for bundle_name in dict_values :
    #     try :
    #         if dict_values[ bundle_name ] == -1 :
    #             continue
    #         else :
    #             measure_values_1.append( float( dict_values[ bundle_name ][
    #                                                                   "a1" ] ) )
    #             measure_values_2.append( float( dict_values[ bundle_name ][
    #                                                                   "a2" ] ) )
    #             breakpoints.append( float( dict_values[ bundle_name ][
    #                                                           "breakpoint" ] ) )
    #     except :
    #         measure_values_1.append( float( dict_values[ bundle_name ][
    #                                                                   "a1" ] ) )
    #         measure_values_2.append( float( dict_values[ bundle_name ][
    #                                                                   "a2" ] ) )
    #         breakpoints.append( float( dict_values[ bundle_name ][
    #                                                           "breakpoint" ] ) )

    # age_min = min( breakpoints )
    # age_max = max( breakpoints )
    # if age < age_min :
    #     print( f"ERROR : age is smaller than minimum age {age_min}" )
    #     sys.exit( 1 )
    # if age > age_max :
    #     print( f"ERROR : age is greater than maximum age {age_max}" )
    #     sys.exit( 1 )

    # measure_min = min( [ min( measure_values_1 ), min( measure_values_2 ) ] )
    # measure_max = max( [ max( measure_values_1 ), max( measure_values_2 ) ] )

    _streamlines = []
    _measure_values = []
    for _bundle in bundles_names :
        bundle_path = os.path.join( bundles_dir, _bundle )
        _bundle_name = _bundle.replace( ".trk", "" )
        if _bundle_name in dict_values.keys() :
            bundle_data = load_tractogram( bundle_path , "same" )
            nbStreamlines = len( bundle_data )
            try :
                if dict_values[ _bundle_name ] == -1 :
                    continue
                else :
                    for fiber in bundle_data.streamlines :
                        _streamlines.append( fiber )
                    if dict_values[ _bundle_name ][ "breakpoint" ] < age :
                        tmpMeasureValue = dict_values[ _bundle_name ][ "a2" ]
                        # tmpMeasureValue = ( 100 *
                        #             ( dict_values[ _bundle_name ][ "a2" ] /
                        #             measure_population_mean_per_bundle[ bundleName ] ) )

                    else :
                        tmpMeasureValue = dict_values[ _bundle_name ][ "a1" ]
                        # tmpMeasureValue = ( 100 *
                        #             ( dict_values[ _bundle_name ][ "a1" ] /
                        #             measure_population_mean_per_bundle[ bundleName ] ) )

                    _measure_values += nbStreamlines * [ tmpMeasureValue ]
            except :
                for fiber in bundle_data.streamlines :
                    _streamlines.append( fiber )
                if dict_values[ _bundle_name ][ "breakpoint" ] < age :
                    tmpMeasureValue = dict_values[ _bundle_name ][ "a2" ]
                    # tmpMeasureValue = ( 100 *
                    #             ( dict_values[ _bundle_name ][ "a2" ] /
                    #             measure_population_mean_per_bundle[ bundleName ] ) )

                else :
                    tmpMeasureValue = dict_values[ _bundle_name ][ "a1" ]
                    # tmpMeasureValue = ( 100 *
                    #             ( dict_values[ _bundle_name ][ "a1" ] /
                    #             measure_population_mean_per_bundle[ bundleName ] ) )

                _measure_values += nbStreamlines * [ tmpMeasureValue ]

    for _bundle_name in dict_values :
        try :
            if dict_values[ _bundle_name ] == -1 :
                continue
            else :
                _a1 = dict_values[ _bundle_name ][ "a1" ]
                _a2 = dict_values[ _bundle_name ][ "a2" ]
                print( f"{_bundle_name} : {_a1}\t|\t{_a2}" )
        except :
            _a1 = dict_values[ _bundle_name ][ "a1" ]
            _a2 = dict_values[ _bundle_name ][ "a2" ]
            print( f"{_bundle_name} : {_a1}\t|\t{_a2}" )
            # tmpMeasure = measure_population_mean_per_bundle[ _bundle_name ]
            # print( f"{_bundle_name} : {tmpMeasure}" )

    measure_min = min( _measure_values )
    measure_max = max( _measure_values )

    scene = window.Scene()


    _streamlines = np.array( _streamlines )
    _measure_values = np.array( _measure_values )
    zeroInNormalizedSlope = ( - measure_min ) / ( measure_max - measure_min )

    # tmpCounter = 0
    # for tmpValue in np.unique( _measure_values ) :
    #     print( f"{tmpCounter} : {tmpValue}" )
    #     tmpCounter += 1


    print( "\nValue in colormap\tTrue slope value" )
    print( f"0                \t{measure_min}" )
    print( f"1                \t{measure_max}" )
    print( f"{zeroInNormalizedSlope}                \t0" )

    # hue = (0.5, 0.5)  # blue only
    # saturation = (0.0, 1.0)  # black to white
    lut_cmap = actor.colormap_lookup_table( scale_range = ( measure_min,
                                                                 measure_max ) )

    stream_actor2 = actor.line( _streamlines, _measure_values, linewidth=0.5,
                                                      lookup_colormap=lut_cmap )
    # stream_actor2 = actor.line(bundle.streamlines, fa, linewidth=0.1)

    bar = actor.scalar_bar( lut_cmap )

    scene.add(stream_actor2)
    scene.add(bar)
    window.show(scene, title = window_title, size=(600, 600), reset_camera=False)

if __name__ == "__main__" :
    main()
