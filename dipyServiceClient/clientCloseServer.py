#!/usr/bin/python3
import sys, os, subprocess, shutil

import time

import socket

from argparse import ArgumentError, ArgumentParser, RawTextHelpFormatter
from textwrap import dedent


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
python3 clientCloseServer.py \
-p 5000 \
-v 1 \

"""

def get_cmd_line_args():
    """
    Create a command line argument parser and return a dict mapping
    <argument name> -> <argument value>.
    """
    parser = ArgumentParser(
        prog="python3 clientCloseServer.py",
        description=dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-p", "--port",
        type=int, required=True, metavar="<int>",
        help=( "Port of the server to close" ) )


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

def client_program() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    # Port number
    port = inputs[ "port" ]
    if port < 5000 :
        print( f"The port must be greater or equal to 5000, got port {port}" )
        sys.exit()

    ################################## Client ##################################
    host = socket.gethostname()  # as both code is running on same pc

    client_socket = socket.socket()  # instantiate
    client_socket.connect( ( host, port ) )  # connect to the server

    maxBytesReceive = 1024
    closeClient = False
    client_socket.send( "close".encode() )  # send message
    while not closeClient :
        data = client_socket.recv( maxBytesReceive ).decode()  # receive response
        if ( data == "0" ) :
            closeClient = True
            client_socket.close()  # close the connection


if __name__ == '__main__':
    t1 = time.time()
    client_program()
    t2 = time.time()
    print( f"Time ellapsed : {t2 - t1}" )
