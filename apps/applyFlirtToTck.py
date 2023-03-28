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

Apply flirt transform to .tck file

Command example:
----------------
python3 normalizeHCPSubjectToMNI.py \
-m in.tck \
-r ref.nii.gz \
-t transform.mat
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
        prog="python3 normalizeHCPSubjectToMNI.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input",
        type=is_file, required=True, metavar="<path>",
        help=( "Input .tck file" ) )

    required.add_argument(
        "-o", "--output",
        type=str, required=True, metavar="<path>",
        help=( "Input .tck file" ) )

    required.add_argument(
        "-m", "--moving",
        type=is_file, required=True, metavar="<path>",
        help=( "Moving .nii/.nii.gz image (where the tractogram is)" ) )

    required.add_argument(
        "-r", "--reference",
        type=is_file, required=True, metavar="<path>",
        help="Reference image (.nii/.nii.gz)" )

    required.add_argument(
        "-t", "--transform",
        type=is_file, required=True, metavar="<path>",
        help="FSL's flirt transform (.mat)" )

    # Optional arguments
    parser.add_argument(
        "-fsl2ants", "--fsl-to-ants",
        type=str,
        help="Path to the fslToAnts.py.")

    parser.add_argument(
        "-f", "--force",
        action="store_true",
        help=( "Force to overwrite files already existing" ) )

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

    input_tck_path = inputs[ "input" ]
    if not input_tck_path.endswith( ".tck" ) :
        print( f"ERROR : -i input must be a .tck file" )

    output_tck_path = inputs[ "output" ]
    if not output_tck_path.endswith( ".tck" ) :
        print( f"ERROR : -o input must be a .tck file" )

    moving_path = inputs[ "moving" ]
    if ( not moving_path.endswith( ".nii.gz" ) and
                                          not moving_path.endswith( ".nii" ) ) :
        print( f"ERROR : -m input must be a .tck file" )

    reference_path = inputs[ "reference" ]
    if ( not reference_path.endswith( ".nii.gz" ) and
                                       not reference_path.endswith( ".nii" ) ) :
        print( f"ERROR : -r input must be a .tck file" )

    transform_path = inputs[ "transform" ]
    if not transform_path.endswith( ".mat" ) :
        print( f"ERROR : -t input must be a .mat file" )

    fsl2ants_command = inputs[ "fsl_to_ants" ]
    if fsl2ants_command :
        if not os.path.isfile( fsl2ants_command ) :
            print( f"ERROR : -fs2ants argument {fsl2ants_command} does not "
                                                                      "exists" )
    else :
        fsl2ants_command = "fslToAnts.py"

    force = inputs[ "force" ]


    tmpDir = os.path.join( os.path.dirname( transform_path ),
                                                       "tmpDirApplyFlirtToTck" )
    if not os.path.isdir( tmpDir ) :
        os.mkdir( tmpDir )
    ################# Converting fsl transforms to ANTs format #################
    if verbose :
        print( "Converting fsl transforms to ANTs format... " )

    transform_to_ants = os.path.join( tmpDir, os.path.basename(
                                         transform_path ).replace( ".mat", ".tfm" ) )
    if not os.path.isfile( transform_to_ants ) or force :
        convert_fsl_to_ants_command = [ fsl2ants_command,
                                        "-i", transform_path,
                                        "-o", transform_to_ants,
                                        "-r", reference_path,
                                        "-m", moving_path ]
        run_sh_process( cmd = convert_fsl_to_ants_command, shell = True )
    else :
        print( f"File {transform_to_ants} exists, use -f to overwrite" )

    print( "Done" )

    ############################## ANTS to MRtrix ##############################
    ants_to_mrtrix_dir = os.path.join( tmpDir, "ANTs_to_MRtrix" )
    if not os.path.isdir( ants_to_mrtrix_dir ) :
        os.mkdir( ants_to_mrtrix_dir )

    # Initialise warp
    if verbose :
        print( "Initialise warp... " )
    inv_identity_warp = os.path.join( ants_to_mrtrix_dir,
                                                     "inv_identity_warp[].nii" )

    is_init_warp = False
    for i in range( 3 ) :
        if os.path.isfile( inv_identity_warp.replace( "[]", f"{i}" ) ) :
            is_init_warp = True
            break

    if not is_init_warp or force :
        warp_init_command = [ "warpinit",
                              reference_path,
                              inv_identity_warp,
                              "-force" ]
        run_sh_process( cmd = warp_init_command, shell = True )
    else :
        print( f"At least one of the {inv_identity_warp} exists, use -f to "
                                                                   "overwrite" )
    print( "Done" )

    # Apply the transformation to identity warp
    if verbose :
        print( "Apply transformation to identity warp..." )

    inv_mrtrix_warp = os.path.join( ants_to_mrtrix_dir,
                                                       "inv_mrtrix_warp[].nii" )

    is_init_warp_mrtrix = False
    for i in range( 3 ) :
        if os.path.isfile( inv_mrtrix_warp.replace( "[]", f"{i}" ) ) :
            is_init_warp_mrtrix = True
            break

    if not is_init_warp_mrtrix or force :
        for i in range( 3 ) :
            apply_command = [ "antsApplyTransforms",
                              "-d", str( 3 ),
                              "-e", str( 0 ),
                              "-i", inv_identity_warp.replace( "[]", f"{i}" ),
                              "-o", inv_mrtrix_warp.replace( "[]", f"{i}" ),
                              "-r", moving_path,
                              "-t", f"[{transform_to_ants}, 1]",
                              "--default-value", str( 2147483647 ) ]
            run_sh_process( cmd = apply_command, shell = True )

    else :
        print( f"At least one of the {inv_mrtrix_warp} exists, use -f to "
                                                                   "overwrite" )

    print( "Done" )

    # Fix warp
    if verbose :
        print( "Fix warp" )

    inv_mrtrix_warp_corrected = os.path.join( ants_to_mrtrix_dir,
                                               "inv_mrtrix_warp_corrected.mif" )

    if not os.path.isfile( inv_mrtrix_warp_corrected ) or force :
        fix_warp_command = [ "warpcorrect",
                             inv_mrtrix_warp,
                             inv_mrtrix_warp_corrected,
                             "-marker", str( 2147483647 ),
                             "-force" ]
        run_sh_process( cmd = fix_warp_command, shell = True )
    else :
        print( f"File {inv_mrtrix_warp_corrected} exists, use -f to overwrite" )

    print( "Done" )


    # Transform tracks file
    if not os.path.isfile( output_tck_path ) or force :
        apply_warp_mrtrix = [ "tcktransform",
                              input_tck_path,
                              inv_mrtrix_warp_corrected,
                              output_tck_path,
                              "-force" ]
        run_sh_process( cmd = apply_warp_mrtrix, shell = True )

    if os.path.isdir( tmpDir ) :
        shutil.rmtree( tmpDir )

if __name__ == "__main__" :
    main()
