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
python3 computeCentroids.py \
-i input.bundles \
-o output.bundles \
-r reference.nii \
-thr 5.0 \
-nbPoints 15 \
-mnc 200 \
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
        prog="python3 clientComputeCentroids.py",
        description=dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group( "required arguments" )
    required.add_argument(
        "-i", "--input",
        type=is_file, required=True, metavar="<path>",
        help="Path to the input tractogram" )

    required.add_argument(
        "-o", "--output",
        type=str, required=True, metavar="<path>",
        help="Saving path for the centroids" )

    required.add_argument(
        "-r", "--reference",
        type=is_file, required=True, metavar="<path>",
        help="Reference image from which extract resolution and volume size" )

    # Optional arguments
    parser.add_argument(
        "-thr", "--threshold",
        type=float, metavar="<float>",
        help="QuickBundles treshold ( default = 5.0 )" )
    parser.add_argument(
        "-nbPoints", "--number-points",
        type=int, metavar="<int>",
        help="QuickBundles resample number of points ( default = 21 )" )
    parser.add_argument(
        "-mnc", "--max-nb-clusters",
        type=int, metavar="<int>",
        help="QuickBundles max number of clusters ( default = 200 )" )
    parser.add_argument(
        "-cv", "--convert-bundles",
        type=is_file, metavar="<path>",
        help="Path to the ConvertBundleFormat binary" )
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

def client_program() :
    t_init = time.time()
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    input_filename = inputs[ "input" ]

    output_filename = inputs[ "output" ]

    reference_anatomy = inputs[ "reference" ]

    qb_threshold = inputs[ "threshold" ]
    if not qb_threshold :
        qb_threshold = 5.0

    nbPoints = inputs[ "number_points" ]
    if not nbPoints :
        nbPoints = 21

    max_nb_clusters = inputs[ "max_nb_clusters" ]
    if not max_nb_clusters :
        max_nb_clusters = 200

    convert_bundle_format = inputs[ "convert_bundles" ]
    if not convert_bundle_format :
        convert_bundle_format = ( "ConvertBundleFormat" )

    logFilePath = inputs[ "log_file" ]
    if logFilePath :
        sys.stdout = open( os.devnull, 'w' )

    # Port number
    port = inputs[ "port" ]
    if port < 5000 :
        print( f"The port must be greater or equal to 5000, got port {port}" )
        sys.exit()

    # Sanity checks
    """
    if not input_filename.endswith( ".bundles" ) :
        print( f"{bcolors.FAIL}ERROR : {bcolors.ENDC} the only supported input "
                                                          "format is .bundles" )
        sys.exit( 1 )

    if not output_filename.endswith( ".bundles" ) :
        print( f"{bcolors.FAIL}ERROR : {bcolors.ENDC} the only supported output"
                                                         " format is .bundles" )
        sys.exit()
    """
    argsToSend = [ input_filename, output_filename, reference_anatomy,
                                        qb_threshold, nbPoints, max_nb_clusters,
                                                         convert_bundle_format ]
    argsToSend = [ "computeCentroids" ] + argsToSend
    argsToSend = [ str( tmp ) for tmp in argsToSend ]

    # print( "xxxxxxxxxxxxxxxxxxxxxxxxxxx" )
    # print( " ".join( argsToSend ) )
    # print( "xxxxxxxxxxxxxxxxxxxxxxxxxxx" )

    ################################## Client ##################################
    host = socket.gethostname()  # as both code is running on same pc
    port = 5000  # socket server port number

    client_socket = socket.socket()  # instantiate
    print( "Connecting to server... ", end = "" )
    client_socket.connect( ( host, port ) )  # connect to the server
    print( "Done" )

    maxBytesReceive = 1024
    closeClient = False
    print( "Sending message to server... ", end = "" )
    client_socket.send( "\t".join( argsToSend ).encode() )  # send message
    print( "Done" )
    timeOut = 50 # In s
    duration = 0
    t1 = time.time()
    while not closeClient :
        data = client_socket.recv( maxBytesReceive ).decode()
        # if ( duration > timeOut or os.path.isfile( output_filename ) ) :
        if ( duration > timeOut or data == "-1" ) :
            closeClient = True
            client_socket.close()  # close the connection
            printInLog( "computeCentroids : closed connection because "
                        "of timeout or error in dipyServer.py", logFilePath )

        if ( data == "0" ) :
            closeClient = True
            client_socket.close()  # close the connection

        duration = time.time() - t1

    if ( data == "0" ) :
        print( "Connection closed by client" )

    executionTime = time.time() - t_init
    if ( logFilePath ) :
        logFile =  open( logFilePath, "a" )
    else :
        logFile = sys.stdout

    # print( " ".join( argsToSend ) + f"\nDuration : {executionTime} s\n",
    #                                                             file = logFile )
    print( " ".join( argsToSend ) + f"\nDuration : {duration} s\n",
                                                                file = logFile )

    if ( logFilePath ) :
        logFile.close()

if __name__ == '__main__':
    client_program()
