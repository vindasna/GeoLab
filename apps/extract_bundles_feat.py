#!/usr/bin/python3
import sys, os, subprocess, shutil

import time
import numpy as np
import h5py

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
python3 countSuccessfulSubjects.py \
-in inputDir
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
        prog="python3 countSuccessfulSubjects.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input",
        type=is_file, required=True, metavar="<path>",
        help=( "Input tractogram in .bundles format" ) )

    required.add_argument(
        "-o", "--output",
        type=str, required=True, metavar="<path>",
        help=( "Output file path" ) )

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


def readBundlesFile( bundle_filename, verbose ):

    if verbose > 1 :
        print (f'Reading ... ... {bundle_filename}')

    # Checking extension bundle_filename
    if bundle_filename.endswith( ".bundles" ) :
        bundle_data_filename = bundle_filename.replace( ".bundles",
                                                                 "bundlesdata" )
    else if ( bundle_filename.endswith( ".bundlesdata" ) ) :
        bundle_data_filename = bundle_filename
        bundle_filename = bundle_filename.replace( ".bundlesdata", "bundles" )
    else :
        print( "ERROR : the only tractogram format supported is .bundles/"
                                                                ".bundlesdata" )
        sys.exit( 1 )

    ns = dict()
    exec(open(bundle_filename).read(), ns)
    curves_count = (ns[ 'attributes' ][ 'curves_count' ])
    labels = (ns[ 'attributes' ][ 'bundles' ])
    try:
        resolution = [ ns[ 'attributes' ][ 'resolutionX' ],
                        ns[ 'attributes' ][ 'resolutionY' ],
                        ns[ 'attributes' ][ 'resolutionZ' ] ]
    except:
        print("Not resolution info in .bundles file...")
        resolution = [-1, -1, -1]
    try:
        size = [ ns[ 'attributes' ][ 'sizeX' ],
                    ns[ 'attributes' ][ 'sizeY' ],
                    ns[ 'attributes' ][ 'sizeZ' ] ]
    except:
        print("Not size info in .bundles file...")
        size = [-1, -1, -1]


    bundle = []
    nPoints = []
    if curves_count > 0:
       f = open(bundle_data_filename,'rb')
       for curve in range(curves_count):
           p = np.fromfile(f, dtype=np.int32, count=1)
           nPoints.append(p[0])
           bundle.append(np.fromfile(f, dtype=np.float32, count=p[0]*3))
       f.close()

       if int(max(nPoints)) == int(min(nPoints)):
          n_points = int(max(nPoints))
       else:
             print ("Number of points per curve not equal...")
             n_points = nPoints
    else:
          print( f"ERROR : problem reading {bundle_filename} or its "
                 f"corresponding .bundlesdata file" )
          sys.exit( 1 )

    return bundle, curves_count, n_points, labels, resolution, size


def convertBundleVectorToMatrix( bundle, curves_count, nbPoints, bundle_name ) :
    bundle_matrix = np.array( bundle ).reshape( ( curves_count, nbPoints, 3 ) )
    return( bundle_matrix )

def main() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    input_tractogram_path = inputs[ "input" ]

    if ( not input_tractogram_path.endswith( ".bundles" ) and
                        not input_tractogram_path.endswith( ".bundlesdata" ) ) :
        print( f"ERROR : the only input tractogram supported is .bundles" )
        sys.exit( 1 )


    output_file = inputs[ "output" ]
    if ( not input.endswith( ".h5" ) ) :
        print( f"ERROR : the only output tractogram supported is .h5" )
        sys.exit( 1 )


    output_feat_marix = f"{output_file}"

    # Transforming .bundles -> .h5
    if verbose :
        print( f"Transforming .bundles -> .h5" )

    bundle, curves_count, nbPoints, _, _, _ = readBundlesFile(
                                                          input_tractogram_path,
                                                          verbose )
    bundle_name = os.path.basename( input_tractogram_path )
    bundle_matrix = convertBundleVectorToMatrix( bundle, curves_count,
                                                         nbPoints, bundle_name )
    bundle_matrix = np.array( bundle_matrix, dtype = "float32" )
    bundle_matrix = np.reshape( bundle_matrix,
                                              ( curves_count, nbPoints, 3, 1 ) )



    dt = h5py.special_dtype( vlen=str )
    hf_feat_matrix = h5py.File( output_feat_marix, 'w' )
    hf_feat_matrix.create_dataset( 'feat', data = bundle_matrix )
    hf_feat_matrix.close()

    if verbose :
        print( 'Feature matrix shape:', bundle_matrix.shape )



if __name__ == "__main__" :
    t1 = time.time()
    main()
    elapsed_time = time.time() - t1
    print( f"Elapsed time : {elapsed_time}" )
