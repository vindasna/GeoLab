#!${PYTHON_BINARY}
import sys, os, subprocess, shutil

import time

import csv


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

Convert fsl transform to ants using c3d_affine_tool ( c3d must be installed )

Command example:
----------------
python3 fslToAnts.py \
-i fsl.mat \
-o ants.tfm \
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
        prog="python3 fslToAnts.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input",
        type=is_file, required=True, metavar="<path>",
        help="Path to the fsl .mat transform" )

    required.add_argument(
        "-o", "--output",
        type=str, required=True, metavar="<path>",
        help="Path where to save ants .tfm transform" )

    required.add_argument(
        "-r", "--reference",
        type=is_file, required=True, metavar="<path>",
        help="Reference image for fsl" )

    required.add_argument(
        "-m", "--moving",
        type=is_file, required=True, metavar="<path>",
        help="Moving image for fsl" )

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

    fsl_transform_path = inputs[ "input" ]
    ants_transform_path = inputs[ "output" ]
    reference_image = inputs[ "reference" ]
    moving_image = inputs[ "moving" ]


    # Converting
    print( "Converting... ", end ="" )
    command = [ "c3d_affine_tool",
                "-ref", reference_image,
                "-src", moving_image,
                fsl_transform_path,
                "-fsl2ras",
                "-oitk", ants_transform_path ]
    run_sh_process( cmd = command, shell = True )
    print( "Done" )

if __name__ == "__main__":
    main()
