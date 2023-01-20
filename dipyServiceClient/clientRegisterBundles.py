#!/usr/bin/python3
import sys, os, subprocess, shutil

import time


from argparse import ArgumentError, ArgumentParser, RawTextHelpFormatter
from textwrap import dedent

import socket

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

Compute centroids of a bundle

Command example:
----------------
python3 clientRegisterBundles.py \
-m moving.bundles \
-s static.bundles \
-ra reference.nii \
-o output.bundles \
-b moving.bundles \
-n 15 \
-xfm rigid \
-lf outLogFile.txt \
-p 5000 \
-v 1 \

"""

def is_file(filepath):
    """ Check file's existence - argparse 'file' argument.
    """
    if not os.path.isfile(filepath):
        raise ArgumentError("File does not exist: %s" % filepath)
    return filepath

def is_dir(filepath):
    """ Check file's existence - argparse 'dir' argument.
    """
    if not os.path.isdir(filepath):
        raise ArgumentError("Directory does not exist: %s" % filepath)
    return filepath


def get_cmd_line_args():
    """
    Create a command line argument parser and return a dict mapping
    <argument name> -> <argument value>.
    """
    parser = ArgumentParser(
        prog="python3 clientRegisterBundles.py",
        description=dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-m", "--moving",
        type=is_file, required=True, metavar="<path>",
        help=( "Moving streamlines ( .bundles/.tck/.trk )" ) )

    required.add_argument(
        "-s", "--static",
        type=is_file, required=True, metavar="<path>",
        help=( "Static streamlines ( .bundles/.tck/.trk )" ) )

    required.add_argument(
        "-ra", "--reference-anatomy",
        type=is_file, required=True, metavar="<path>",
        help=( "Reference anatomy ( .nii )" ) )

    required.add_argument(
        "-o", "--output",
        type=str, required=True, metavar="<path>",
        help=( "Output path ( .bundles/.tck/.trk )" ) )

    # Optional arguments
    parser.add_argument(
        "-b", "--bundle",
        type=is_file, metavar="<path>",
        help=( "Path to the bundle you want to apply the transform, if "
               "argument not given, apply transform to moving streamlines "
               "file" ) )

    parser.add_argument(
        "-n", "--nbPoints",
        type=int, metavar="<int>", default=21,
        help=( "Number of points to resampling for streamline registration "
               "( default = 21 )" ) )

    parser.add_argument(
        "-xfm", "--xfm-type",
        type=str, metavar="<int>", default="affine", choices=["rigid","affine"],
        help=( "Type of transformation for registering (options : \"affine\", "
               "\"rigid\", default = \"affine\")" ) )

    parser.add_argument(
        "-cv", "--convert-bundles",
        type=is_file, metavar="<path>",
        help=( "Path to the ConvertBundleFormat binary" ) )

    parser.add_argument(
        "-lf", "--log-file",
        type=str, metavar="<path>",
        help="Path where to save logFile" )

    parser.add_argument(
        "-p", "--port",
        type=int, metavar="<int>", default=5000,
        help="Port of the server to connect (default = 5000 )" )

    parser.add_argument(
        "-v", "--verbose",
        type=int, choices=[0, 1, 2], default=0,
        help="Increase the verbosity level: 0 silent, 1 verbose.")


    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)
    verbose = kwargs.pop("verbose")

    return kwargs, verbose

def printInLog( inputString, logFilePath ) :
    if ( logFilePath ) :
        logFile =  open( logFilePath, "a" )
    else :
        logFile = sys.stdout

    print( inputString, file = logFile )

    if ( logFilePath ) :
        logFile.close()


def client_program() :
    t_init = time.time()
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    # Log file
    logFilePath = inputs[ "log_file" ]
    if logFilePath :
        sys.stdout = open(os.devnull, 'w' )

    printInLog( "clientRegisterBundles.py", logFilePath )

    # Reference anatomy
    reference_anatomy = inputs[ "reference_anatomy" ]

    # Atlas neighborhood
    static = inputs[ "static" ]
    if ( not static.endswith( ".bundles" ) and not static.endswith( ".trk" ) and
                                               not static.endswith( ".tck" ) ) :
        printInLog( "ERROR : Static streamlines must be in .bundles/.trk/.tck "
                                                         "format", logFilePath )
        exit( 1 )

    # Tractogram neighborhood
    moving = inputs[ "moving" ]
    if ( not moving.endswith( ".bundles" ) and not moving.endswith( ".trk" ) and
                                               not moving.endswith( ".tck" ) ) :
        printInLog( "ERROR : Moving streamlines must be in .bundles/.trk/.tck ",
                                                         "format", logFilePath )
        exit( 1 )

    # Bundle to which apply transform
    bundle_filename = inputs[ "bundle" ]
    isBundle = False
    if ( bundle_filename ) :
        isBundle = True
        if not os.path.isfile( bundle_filename ) :
            printInLog( f"ERROR : -b argument file {bundle_filename} does not "
                                                         "exists", logFilePath )
            exit( 1 ) ;

    # Output
    output_moved = inputs[ "output" ]
    if ( not output_moved.endswith( ".bundles" ) and
                                         not output_moved.endswith( ".trk" ) and
                                         not output_moved.endswith( ".tck" ) ) :
        printInLog( "ERROR : Output streamlines must be in .bundles/.trk/.tck ",
                                                         "format", logFilePath )
        exit( 1 )

    # Transformation type
    xfmType = inputs[ "xfm_type"]

    # Number of points for alignement
    nbPoints = inputs[ "nbPoints"]

    # ConvertBundleFormat binary
    convert_bundle_format = inputs[ "convert_bundles" ]
    if not convert_bundle_format :
        convert_bundle_format = "ConvertBundleFormat"

    printInLog( f"ConvertBundleFormat binary path : {convert_bundle_format}",
                                                                   logFilePath )

    # Port number
    port = inputs[ "port" ]
    if port < 5000 :
        printInLog( f"The port must be greater or equal to 5000, got port {port}",
                                                                     logFilePath )
        sys.exit()

    # Input parameters
    if verbose :
        printInLog( f"Trasnform type : {xfmType}\nNumber of points : {nbPoints}",
                                                                     logFilePath )

    argsToSend = [ static, moving, bundle_filename, reference_anatomy,
                   output_moved, xfmType, nbPoints, convert_bundle_format ]
    argsToSend = [ "registerBundles" ] + argsToSend
    argsToSend = [ str( tmp ) for tmp in argsToSend ]

    printInLog( " ".join( argsToSend ), logFilePath )

    ################################## Client ##################################
    printInLog( f"Connecting to server in port {port}", logFilePath )
    host = socket.gethostname()  # as both code is running on same pc

    client_socket = socket.socket()  # instantiate
    client_socket.connect( ( host, port ) )  # connect to the server

    maxBytesReceive = 1024
    closeClient = False
    time_beginning = time.time()
    client_socket.send( "\t".join( argsToSend ).encode() )  # send message
    timeOut = 50 # In s
    duration = 0
    printInLog( "Connection successful", logFilePath )
    # tmpDebug = 0
    t1 = time.time()
    while not closeClient :
        # if tmpDebug < 10 :
        #     printInLog( tmpDebug, logFilePath )
        data = client_socket.recv( maxBytesReceive ).decode()
        # printInLog( f"Data received : {data}", logFilePath )
        if ( duration > timeOut or data == "-1" ) :
            closeClient = True
            client_socket.close()  # close the connection
            printInLog( "registerBundles : closed connection because "
                        "of timeout or error in dipyServer.py", logFilePath )

        if ( data == "0" ) :
            closeClient = True
            client_socket.close()  # close the connection

        duration = time.time() - t1
        # tmpDebug += 1
        # printInLog( tmpDebug, logFilePath )

    if ( data == "0" ) :
        printInLog( "Connection closed by client", logFilePath )

    executionTime = time.time() - t_init
    # printInLog( f"\nDuration : {executionTime} s\n", logFilePath )
    printInLog( f"\nDuration : {duration} s\n", logFilePath )


if __name__ == '__main__':
    client_program()
