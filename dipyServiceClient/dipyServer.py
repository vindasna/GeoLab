#!/usr/bin/python3
import sys, os, subprocess, shutil

import time

from argparse import ArgumentError, ArgumentParser, RawTextHelpFormatter
from textwrap import dedent

import numpy as np
from dipy.io.streamline import load_tractogram, save_trk, save_tractogram
from dipy.tracking.streamline import Streamlines, set_number_of_points
from dipy.segment.clustering import QuickBundles
from dipy.io.stateful_tractogram import Space, StatefulTractogram

from dipy.align.streamlinear import StreamlineLinearRegistration

from dipy.segment.featurespeed import ResampleFeature # DIPY 1.5
from dipy.segment.metric import AveragePointwiseEuclideanMetric

from nibabel.streamlines.array_sequence import ArraySequence
from nibabel.streamlines.trk import TrkFile

import socket

# from threading import Thread

from multiprocessing import Process

from setproctitle import setproctitle

# Set process name
setproctitle( "dipyServer.py" )

# Global variables
nbTotalConnections = 0

#------------------------------------------------------------------------------#
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

Launch dipy service

Command example:
----------------
python3 dipyServer.py \
-log-file logFilePath.txt \
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
        prog="python3 computeTractogramMrtrix.py",
        description=dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group( "required arguments" )

    # Optional arguments
    parser.add_argument(
        "-lf", "--log-file",
        type=str, metavar="<path>",
        help="Path where to save logFile" )

    parser.add_argument(
        "-v", "--verbose",
        type=int, choices=[0, 1, 2], default=0,
        help="Increase the verbosity level: 0 silent, 1 verbose.")


    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)
    verbose = kwargs.pop("verbose")

    return kwargs, verbose

def close_server( server, logFile ) :
    server.shutdown( socket.SHUT_RDWR )
    server.close()
    print( "Server closed", file = logFile )

def run_sh_process( cmd, shell = False ) :
    command_name = cmd[ 0 ]
    if shell :
        cmd = " ".join( cmd )
    process = subprocess.Popen( cmd, stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE, shell = shell,
                                                              encoding="utf-8" )
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


def computeCentroids( input_filename, output_filename, reference_anatomy,
                      qb_threshold, nbPoints, max_nb_clusters,
                                                 convert_bundle_format, conn ) :

    if ( input_filename.endswith( ".bundles" ) and
                                  not output_filename.endswith( ".bundles" ) ) :
        print( f"ERROR : input and output bundles must have the same format" )
        conn.send( "-1".encode() )
        conn.close()
        sys.exit( 1 )
    elif ( input_filename.endswith( ".trk" ) and
                                  not output_filename.endswith( ".trk" ) ) :
        print( f"ERROR : input and output bundles must have the same format" )
        conn.send( "-1".encode() )
        conn.close()
        sys.exit( 1 )
    elif ( input_filename.endswith( ".tck" ) and
                                  not output_filename.endswith( ".tck" ) ) :
        print( f"ERROR : input and output bundles must have the same format" )
        conn.send( "-1".encode() )
        conn.close()
        sys.exit( 1 )

    print( "Computing centroids" )
    # Converting .bundles to .trk
    bundleName = os.path.basename( input_filename ).replace( ".bundles", "")
    tmp_dir = os.path.join( os.path.dirname( output_filename ),
                                                        f"tmpDir_{bundleName}" )
    if not os.path.isdir( tmp_dir ) :
        os.mkdir( tmp_dir )

    if ( input_filename.endswith( ".bundles" ) ) :
        tractogram_trk_path = os.path.join( tmp_dir,  os.path.basename(
                                input_filename ).replace( ".bundles", ".trk" ) )
        tractogram_trk_minf_path = os.path.join( tmp_dir,  os.path.basename(
                               input_filename ).replace( ".bundles", ".minf" ) )
        convert_bundles_command = [ convert_bundle_format,
                                    "-i", input_filename,
                                    "-r", reference_anatomy,
                                    "-o", tractogram_trk_path ]
        run_sh_process( cmd = convert_bundles_command, shell = True )
    elif ( input_filename.endswith( ".trk" ) or input_filename.endswith( ".tck" ) ) :
        tractogram_trk_path = input_filename
        if input_filename.endswith( ".trk" ) :
            tractogram_trk_minf_path = input_filename.replace( ".trk", ".minf" )
        else :
            tractogram_trk_minf_path = input_filename.replace( ".tck", ".minf" )

    if not os.path.isfile( tractogram_trk_path ) :
        print( f"{bcolors.FAIL}ERROR : {bcolors.ENDC} could not convert input "
                                                   "tractogram to .trk format" )
        if os.path.isdir( tmp_dir ) :
            shutil.rmtree( tmp_dir )

        conn.send( "-1".encode() )
        conn.close()
        sys.exit( 1 )


    # Reading streamlines.
    # print( "Reading streamlines... ", end = " ", flush = True )
    if ( tractogram_trk_path.endswith( ".trk" ) ) :
        tractogram = load_tractogram( tractogram_trk_path, 'same',
                                                      bbox_valid_check = False )
    else :
        tractogram = load_tractogram( tractogram_trk_path, reference_anatomy,
                                                      bbox_valid_check = False )
    streamlines = tractogram.streamlines
    # print( "Done", flush = True )

    # Streamlines will be resampled to 21 points on the fly.
    # print( "Computing centroids... ", end = " ", flush = True )
    feature = ResampleFeature( nb_points = nbPoints )
    metric = AveragePointwiseEuclideanMetric( feature = feature )  # a.k.a. MDF
    qb = QuickBundles(threshold = qb_threshold, metric = metric,
                                             max_nb_clusters = max_nb_clusters )
                                                       # max_nb_clusters = 200 )
                                                       # max_nb_clusters = 500 )

    clusters = qb.cluster( streamlines )

    # print( "Nb. clusters:", len( clusters ) )

    centroids = []
    for i in range( len( clusters ) ) :
        centroids.append( clusters[ i ].centroid )

    centroids = np.array( centroids, dtype = np.ndarray )

    streamlines = ArraySequence( iterable = centroids )

    trk_file = TrkFile( streamlines )

    # tractogram_trk = StatefulTractogram( streamlines,
    #                                      reference_anatomy,
    #                                      Space.RASMM )
    if ( tractogram_trk_path.endswith( ".trk" ) ) :
        tractogram_trk = StatefulTractogram( streamlines, tractogram_trk_path,
                                                                   Space.RASMM )
    else :
        tractogram_trk = StatefulTractogram( streamlines, reference_anatomy,
                                                                   Space.RASMM )

    if ( output_filename.endswith( ".bundles" ) ) :
        output_trk = os.path.join( tmp_dir,
                 os.path.basename( output_filename ).replace( ".bundles",
                                                                      ".trk" ) )
    else :
        output_trk = output_filename
    # save_trk( tractogram_trk, output_trk, bbox_valid_check = False )
    save_tractogram( tractogram_trk, output_trk, bbox_valid_check = False )

    if not os.path.isfile( output_trk ) :
        print( f"{bcolors.FAIL}ERROR : {bcolors.ENDC} could not save "
                                                                   "centroids" )
        if os.path.isdir( tmp_dir ) :
            shutil.rmtree( tmp_dir )

        conn.send( "-1".encode() )
        conn.close()
        sys.exit( 1 )

    if ( output_filename.endswith( ".bundles" ) ) :
        output_trk_minf = output_trk.replace( ".trk", ".minf" )
        if ( tractogram_trk_minf_path != output_trk_minf
                              and os.path.isfile( tractogram_trk_minf_path ) ) :
            shutil.copy2( tractogram_trk_minf_path, output_trk_minf )
        else :
            print( f"{bcolors.FAIL}ERROR : {bcolors.ENDC} file "
                   f"{tractogram_trk_minf_path} does not exists" )
            if os.path.isdir( tmp_dir ) :
                shutil.rmtree( tmp_dir )

            conn.send( "-1".encode() )
            conn.close()
            sys.exit( 1 )
        convert_bundles_command = [ convert_bundle_format,
                                    "-i", output_trk,
                                    "-r", reference_anatomy,
                                    "-o", output_filename ]
        run_sh_process( cmd = convert_bundles_command, shell = True )
        if not os.path.isfile( output_filename ) :
            print( f"{bcolors.FAIL}ERROR : {bcolors.ENDC} could not convert "
                "output tractogram to .trk format" )

        if os.path.isdir( tmp_dir ) :
            shutil.rmtree( tmp_dir )

    print( f"Closing connection {conn}" )
    # conn.shutdown()
    conn.send( "0".encode() )
    conn.close()

    return( 0 )


def registerBundles( static, moving, bundle_filename, reference_anatomy,
                     output_moved, xfmType, nbPoints, convert_bundle_format,
                                                           conn, logFilePath ) :
    # For the moment no need to check bundle_filename as it is the same as static
    if ( static.endswith( ".bundles" ) and ( not moving.endswith( ".bundles" )
                                or not output_moved.endswith( ".bundles" ) ) ) :
        print( f"ERROR : static, moving and output_moved bundles must have the "
               f"same format" )
        conn.send( "-1".encode() )
        conn.close()
        sys.exit( 1 )
    elif ( static.endswith( ".trk" ) and ( not moving.endswith( ".trk" ) or
                                       not output_moved.endswith( ".trk" ) ) ) :
        print( f"ERROR : static, moving and output_moved bundles must have the "
               f"same format" )
        conn.send( "-1".encode() )
        conn.close()
        sys.exit( 1 )
    elif ( static.endswith( ".tck" ) and ( not moving.endswith( ".tck" ) or
                                       not output_moved.endswith( ".tck" ) ) ) :
        print( f"ERROR : static, moving and output_moved bundles must have the "
               f"same format" )
        conn.send( "-1".encode() )
        conn.close()
        sys.exit( 1 )

    isBundle = False
    if os.path.isfile( bundle_filename ) :
        isBundle = True


    print( "Registering bundles" )
    if not os.path.isfile( bundle_filename ) :
        print( f"ERROR : -b argument file {bundle_filename} does not exists" )
        conn.send( "-1".encode() )
        conn.close()
        sys.exit( 1 )

    if ( static.endswith( ".bundles" ) ) :
        # Converting .bundles -> .trk
        # print( "Converting static to trk...", end = " " )
        static_trk_filename = static.replace( ".bundles", ".trk" )
        static_trk_minf_filename = static.replace( ".bundles", ".minf" )
        convert_bundle_command = [ convert_bundle_format,
                                   "-i", static,
                                   "-o", static_trk_filename,
                                   "-r", reference_anatomy,
                                   "-v", str( 0 ) ]
        run_sh_process( cmd = convert_bundle_command, shell = True )
        # print( "Done" )

        # print( "Converting moving to trk...", end = " " )
        moving_trk_filename = moving.replace( ".bundles", ".trk" )
        moving_trk_minf_filename = moving.replace( ".bundles", ".minf" )
        convert_bundle_command = [ convert_bundle_format,
                                   "-i", moving,
                                   "-o", moving_trk_filename,
                                   "-r", reference_anatomy,
                                   "-v", str( 0 ) ]
        run_sh_process( cmd = convert_bundle_command, shell = True )
        # print( "Done" )
    elif ( static.endswith( ".trk" ) or static.endswith( ".tck" ) ) :
        static_trk_filename = static
        moving_trk_filename = moving
        if static.endswith( ".trk" ) :
            static_trk_minf_filename = static.replace( ".trk", ".minf" )
            moving_trk_minf_filename = moving.replace( ".trk", ".minf" )
        else :
            static_trk_minf_filename = static.replace( ".tck", ".minf" )
            moving_trk_minf_filename = moving.replace( ".tck", ".minf" )
    else :
        print( f"ERROR : the only supported formats are .bundles/.bundlesdata"
               f", .trk and .tck" )
        conn.send( "-1".encode() )
        conn.close()
        sys.exit( 1 )


    # Loading tractograms
    # static_trk = load_tractogram( static_trk_filename, reference_anatomy,
    #                                                   bbox_valid_check = False )
    if ( static_trk_filename.endswith( ".trk" ) ) :
        static_trk = load_tractogram( static_trk_filename, 'same',
                                                      bbox_valid_check = False )
    else :
        static_trk = load_tractogram( static_trk_filename, reference_anatomy,
                                                      bbox_valid_check = False )

    static_trk_streamlines = set_number_of_points( static_trk.streamlines,
                                                                      nbPoints )
    # moving_trk = load_tractogram( moving_trk_filename, reference_anatomy,
    #                                                   bbox_valid_check = False )
    if ( moving_trk_filename.endswith( ".trk" ) ) :
        moving_trk = load_tractogram( moving_trk_filename, 'same',
                                                      bbox_valid_check = False )
    else :
        moving_trk = load_tractogram( moving_trk_filename, reference_anatomy,
                                                      bbox_valid_check = False )

    moving_trk_streamlines = set_number_of_points( moving_trk.streamlines,
                                                                      nbPoints )

    isBundleFormat = False
    if isBundle :
        if bundle_filename.endswith( ".bundles" ) :
            bundle_filename_trk = bundle_filename.replace( ".bundles", ".trk" )
            bundle_filename_minf_trk = bundle_filename.replace( ".bundles",
                                                                       ".minf" )
            isBundleFormat = True
        elif bundle_filename.endswith( ".bundlesdata" ) :
            bundle_filename_trk = bundle_filename.replace( ".bundlesdata",
                                                                        ".trk" )
            bundle_filename_minf_trk = bundle_filename.replace( ".bundlesdata",
                                                                       ".minf" )
            isBundleFormat = True
        elif bundle_filename.endswith( ".trk" ) :
            bundle_filename_trk = bundle_filename
            bundle_filename_minf_trk = bundle_filename.replace( ".trk",
                                                                       ".minf" )
        elif bundle_filename.endswith( ".tck" ) :
            bundle_filename_trk = bundle_filename
            bundle_filename_minf_trk = bundle_filename.replace( ".tck",
                                                                       ".minf" )
        else :
            print( f"ERROR : the only supported formats are "
                   f".bundles/.bundlesdata, .trk and .tck" )
            conn.send( "-1".encode() )
            conn.close()
            sys.exit( 1 )


        # print( "Converting output to bundles...", end = " " )
        if isBundleFormat :
            convert_bundle_command = [ convert_bundle_format,
                                       "-i", bundle_filename,
                                       "-o", bundle_filename_trk,
                                       "-r", reference_anatomy ]
            run_sh_process( cmd = convert_bundle_command, shell = True )
        # print( "Done" )
        # bundleData = load_tractogram( bundle_filename_trk, reference_anatomy,
        #                                               bbox_valid_check = False )
        if ( bundle_filename_trk.endswith( ".trk" ) ) :
            bundleData = load_tractogram( bundle_filename_trk, 'same',
                                                      bbox_valid_check = False )
        else :
            bundleData = load_tractogram( bundle_filename_trk, reference_anatomy,
                                                      bbox_valid_check = False )

        bundleData_streamlines = set_number_of_points(
                                              bundleData.streamlines, nbPoints )


    # Registering
    # print( "Aligning...", end = " ", flush = True )

    optimizationOptions = {'maxcor': 10, 'ftol': 1e-7,
                                'gtol': 1e-5, 'eps': 1e-8,
                                'maxiter': 200}
    srr = StreamlineLinearRegistration( num_threads = -1, x0 = xfmType,
                                 options = optimizationOptions, verbose = True )
    srm = srr.optimize( static = static_trk_streamlines,
                                               moving = moving_trk_streamlines )

    if isBundle :
        aligned = srm.transform( bundleData_streamlines )
    else :
        aligned = srm.transform( moving_trk_streamlines )


    # print( "Done", flush = True )


    # Saving results
    # print( "Saving..." , end = " ", flush = True )

    if ( output_moved.endswith( ".bundles" ) ) :
        output_moved_trk = output_moved.replace( ".bundles", ".trk" )
        output_moved_minf_trk = output_moved.replace( ".bundles", ".minf" )
        if ( bundle_filename_minf_trk != output_moved_minf_trk
                              and os.path.isfile( bundle_filename_minf_trk ) ) :
            shutil.copy2( bundle_filename_minf_trk, output_moved_minf_trk )
        else :
            print( f"{bcolors.FAIL}ERROR : {bcolors.ENDC} file "
                   f"{bundle_filename_minf_trk} does not exists" )

            conn.send( "-1".encode() )
            conn.close()
            sys.exit( 1 )
    elif ( output_moved.endswith( ".trk" ) or
                                             output_moved.endswith( ".tck" ) ) :
        output_moved_trk = output_moved
        if ( output_moved.endswith( ".trk" ) ) :
            output_moved_minf_trk = output_moved.replace( ".trk", ".minf" )
        else :
            output_moved_minf_trk = output_moved.replace( ".tck", ".minf" )


    # moved_tractogram = StatefulTractogram( aligned, reference_anatomy,
    #                                                                Space.RASMM )
    if ( moving_trk_filename.endswith( ".trk" ) ) :
        moved_tractogram = StatefulTractogram( aligned, moving_trk_filename,
                                                                   Space.RASMM )
    else :
        moved_tractogram = StatefulTractogram( aligned, reference_anatomy,
                                                                   Space.RASMM )
    save_tractogram( moved_tractogram, output_moved_trk, False  )
    # print( "Done", flush = True )

    # print( "Converting output to bundles...", end = " " )
    if ( output_moved.endswith( ".bundles" ) ) :
        convert_bundle_command = [ convert_bundle_format,
                                   "-i", output_moved_trk,
                                   "-o", output_moved,
                                   "-r", reference_anatomy,
                                   "-v", str( 0 ) ]
        run_sh_process( cmd = convert_bundle_command, shell = True )
        # print( "Done" )

        # print( "Cleaning...", end = " " )
        if isBundleFormat :
            if os.path.isfile( static_trk_filename ) :
                os.remove( static_trk_filename )

            if os.path.isfile( moving_trk_filename ) :
                os.remove( moving_trk_filename )

            if os.path.isfile( output_moved_trk ) :
                os.remove( output_moved_trk )

            if isBundle :
                if os.path.isfile( bundle_filename_trk ) :
                    os.remove( bundle_filename_trk )

    print( 7 )

    # print( "Done" )

    print( f"Closing connection {conn}" )
    conn.send( "0".encode() )
    conn.close()

    return( 0 )

#------------------------------ Server functions ------------------------------#
def server_program() :
    global nbTotalConnections

    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()
    logFilePath = inputs[ "log_file" ]
    if logFilePath :
        sys.stdout = open( os.devnull, 'w' )
        if os.path.isfile( logFilePath ) :
            os.remove( logFilePath )


    ########################### Initialising server ############################
    if ( logFilePath ) :
        logFile = open( logFilePath, 'a' )
    else :
        logFile = sys.stdout

    # Get the host name
    print( "Getting host name... ", end = "", file = logFile )
    host = socket.gethostname()
    print( f"Done", file = logFile )


    # Creating server socket
    port = 5000 # Ports smaller than 1024 are used for internet protocols
    print( f"Creating server socket... Trying port {port}", end = "\r",
                                                                file = logFile )
    isServerLaunched = False
    while ( not isServerLaunched ) :
        try :
            server_socket = socket.socket()
            server_socket.bind( ( host, port ) )
            isServerLaunched = True
        except :
            port += 1
            print( f"Creating server socket... Trying port {port}", end = "\r",
                                                                file = logFile )
    print( "Done", file = logFile )
    print( f"Server port : {port}", file = logFile )

    # Configure how many clients the server can listen at the same time
    print( "Configuring server... ", end = "", file = logFile )
    max_nb_listeners = 10000
    server_socket.listen( max_nb_listeners )
    print( f"Max number of listeners : {max_nb_listeners}", file = logFile )

    if ( logFilePath ) :
        logFile.close()


    ########################### Listening to clients ###########################

    closeConnection = False
    processList = []
    while not closeConnection :
        if ( logFilePath ) :
            logFile = open( logFilePath, 'a' )
            sys.stderr = logFile
        else :
            logFile = sys.stdout
        print( f"Number of total connections : {nbTotalConnections}",
                                                                file = logFile )
        conn, adress = server_socket.accept() # Accept new connection
        nbTotalConnections += 1
        print( f"Connection from : {adress}", file = logFile )
        maxBytesReceive = 1024
        # Receive data stream (limited to maxBytesReceive bytes)
        data = conn.recv( maxBytesReceive ).decode()
        if data :
            if data == "close" :
                closeConnection = True
                dataToSend = "0"
                conn.send( dataToSend.encode() ) # send data to client
                close_server( server_socket, logFile )
                break
            argumentsReceived = data.split( "\t" )
            print( " ".join( argumentsReceived ) + "\n", file = logFile )
            command = argumentsReceived[ 0 ]
            if command == "computeCentroids" :
                input_filename = argumentsReceived[ 1 ]
                output_filename = argumentsReceived[ 2 ]
                reference_anatomy = argumentsReceived[ 3 ]
                qb_threshold = float( argumentsReceived[ 4 ] )
                nbPoints = int( argumentsReceived[ 5 ] )
                max_nb_clusters = int( argumentsReceived[ 6 ] )
                convert_bundle_format = argumentsReceived[ 7 ]
                p = Process( target = computeCentroids, args = [
                        input_filename, output_filename, reference_anatomy,
                        qb_threshold, nbPoints, max_nb_clusters,
                        convert_bundle_format, conn ], daemon = True )
                p.start()
                processList.append( p )
            elif command == "registerBundles" :
                static = argumentsReceived[ 1 ]
                moving = argumentsReceived[ 2 ]
                bundle_filename = argumentsReceived[ 3 ]
                reference_anatomy = argumentsReceived[ 4 ]
                output_moved = argumentsReceived[ 5 ]
                xfmType = argumentsReceived[ 6 ]
                nbPoints = int( argumentsReceived[ 7 ] )
                convert_bundle_format = argumentsReceived[ 8 ]
                p = Process( target = registerBundles, args = [
                        static, moving, bundle_filename, reference_anatomy,
                        output_moved, xfmType, nbPoints,
                        convert_bundle_format, conn, logFilePath ], daemon = True )
                p.start()
                processList.append( p )
            elif command == "test" :
                conn.send(
                    f"Number total connections : {nbTotalConnections}".encode() )
                time.sleep( 1 )
                conn.send( "0".encode() )
                print( "Connection to client-server : OK", file = logFile )
                conn.close()
        processIndex = 0
        for processTmp in processList :
            if not processTmp.is_alive() :
                processTmp.terminate()
                del processList[ processIndex ]
            processIndex += 1

        if ( logFilePath ) :
            logFile.close()





if __name__ == '__main__':
    server_program()
