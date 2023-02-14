#!/usr/bin/python3
import sys, os, subprocess, shutil

import pickle

import numpy as np

import matplotlib.pyplot as plt

from sklearn import metrics
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error, r2_score
import statsmodels.api as sm
import statsmodels.formula.api as smf

import piecewise_regression

import pwlf

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
measures_options = [ "FA", "MD", "OD", "ICVF", "ISOVF" ]
#-----------------------------------------#


DOC = """
---------------------------

Plot the statistical analysis of a bundle

Command example:
----------------
python3 launchTractogramTGCC.py \
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
        prog="python3 getAgeSubject.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-ukb-dti", "--ukb-dti-dir",
        type=is_dir, required=True, metavar="<path>",
        help=( "UKBiobank directory with dti measures" ) )

    required.add_argument(
        "-t", "--tsv",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the .tsv file with information of the subjects in "
                                                               "ukb-dti-dir" ) )
    required.add_argument(
        "-bc", "--baseline-characteristics",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the .tsv file with baseline characteristics of "
                                                   "subjects in ukb-dti-dir" ) )

    required.add_argument(
        "-m", "--measure",
        type=str, required=True, metavar="<string>",
        help=( "Name of the measure to plot (options : FA, MD, OD, ICVF, "
                                                                    "ISOVF)" ) )

    required.add_argument(
        "-bundle", "--bundle-name",
        type=str, required=True, metavar="<string>",
        help=( "Name of the bundle to plot" ) )

    required.add_argument(
        "-o", "--out",
        type=str, required=True, metavar="<path>",
        help=( "Output directory where to save plots" ) )

    # Optional arguments
    parser.add_argument(
        "-v", "--verbose",
        type=int, choices=[0, 1, 2], default=0,
        help="Increase the verbosity level: 0 silent, 1 verbose." )


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

def readBaselineCharacteristics( path, verbose ) :
    data_dict = {}
    i = 0
    print( f"Reading {path}" )
    with open( path, "r" ) as f :
        for line in f :
            try :
                participant_id, sex, df_34, df_52, df_189, df_21022 = (
                                                            line.split( "\t" ) )
            except :
                nbColumns = len( line.split( "\t" ) )
                print( f"ERROR in readBaselineCharacteristics : expected 6 "
                       f"tabulated columns, got {nbColumns}" )
                sys.exit( 1 )

            data_dict[ participant_id ] = { "sex" : sex, "34" : df_34,
                                            "52" : df_52, "189" : df_189,
                                            "21022" : df_21022 }
    if data_dict[ "participant_id" ][ "sex" ] != "31" :
        print( f"ERROR in readBaselineCharacteristics : problem with title of "
               f"age datafield" )
        sys.exit( 1 )

    return( data_dict )

def short_summary( est ) :
    # return HTML( est.summary().tables[1].as_html() )
    return( est.summary().tables[1] )

def getPValuesStatsModel( fitted_model, exog_names ) :
    pvalues = fitted_model.pvalues
    params = fitted_model.params
    if len( exog_names ) != len( params ) :
        print( "ERROR in getPValuesStatsModel : the length of the params "
               "in fitted_model is different than the length of exog_names" )
        sys.exit( 1 )

    line_to_print = "Name      \tValue     \tp-value  \n"
    for i in range( len( exog_names ) ) :
        line_to_print += ( f"{exog_names[ i ]:10}\t{params[ i ]:10.4e}\t"
                                                     f"{pvalues[ i ]:10.4e}\n" )

    print( line_to_print )

def getPValueInteractionAgeSex( fitted_model, exog_names ) :
    pvalues = fitted_model.pvalues
    if len( exog_names ) != len( pvalues ) :
        print( "ERROR in getPValueInteractionAgeSex : the length of the "
               "pvalues in fitted_model is different than the length of "
                                                                  "exog_names" )
        sys.exit( 1 )

    for i in range( len( exog_names ) ) :
        if exog_names[ i ] == "age:sex" :
            return( pvalues[ i ] )

    return( None )


def main() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    ukb_dti_dir = inputs[ "ukb_dti_dir" ]

    tsv_file = inputs[ "tsv" ]

    baseline_characteristics_path = inputs[ "baseline_characteristics" ]

    bundle_name = inputs[ "bundle_name" ]

    output_dir = inputs[ "out" ]
    if not os.path.isdir( output_dir ) :
        os.mkdir( output_dir )

    measure = inputs[ "measure" ]
    if measure not in measures_options :
        print( "ERROR : the input of -m/--measure must be FA, MD, OD, ICVF or "
                                                                       "ISOVF" )
        sys.exit( 1 )

    subjects = os.listdir( ukb_dti_dir )
    nbSubjects = len( subjects )


    out_data_dict_path = os.path.join( output_dir, f"{measure}AndAge.txt" )
    print( f"Fusing DTI mean of {measure} in {out_data_dict_path}" )
    getFileForPlotMiguelCommand = [ "getFileForPlotsMiguel",
                                    "-ukb-dti", ukb_dti_dir,
                                    "-m", measure,
                                    "-tsv",  tsv_file,
                                    "-o", out_data_dict_path ]
    if ( not os.path.isfile( out_data_dict_path ) ) :
        run_sh_process( getFileForPlotMiguelCommand, shell = True )

    data_dict = readFileWithMeansMeasure( out_data_dict_path )
    nbBundles = len( data_dict )

    baseline_characteristics_dict = readBaselineCharacteristics(
                                        baseline_characteristics_path, verbose )

    # To have reproductible results
    np.random.seed( 100 )

    count_bundle = 1
    print( "Saving plot and models..." )
    maxAbsSlope = 0
    maxAbsSlopeBundle = ""
    nbBundlesWithDti = 0
    bundlesWithSignificantAgeSexInteraction = []
    slope_dict = {}
    for bundle in data_dict :
        if bundle != bundle_name :
            continue


        age_data = data_dict[ bundle ][ 0 ]
        measure_data = data_dict[ bundle ][ 1 ]

        X = np.array( age_data )
        y = np.array( measure_data )

        totalNbSubjectsForBundle = len( y )
        print( f"Total number of subject for bundle : "
                                                 f"{totalNbSubjectsForBundle}" )

        age_min = np.min( X )
        age_max = np.max( X )
        print( f"Age min : {age_min}\t|\tAge max : {age_max}" )

        ########################################################################
        sampling = round( len( X ) / 4 )
        if sampling == 1 : # To avoid division by 0 later
            sampling = 2
        minSamples = X.min()
        maxSamples = X.max()
        samplesX = np.zeros( sampling )
        samplesSizeX = ( ( maxSamples - minSamples ) / ( sampling - 1 ) )
        for i in range( sampling ):
            samplesX[ i ] = ( minSamples + samplesSizeX * i )

        plt.hist( X, samplesX, density = True )
        plt.show()
        plt.clf()
        plt.close()


    print( "\nDone" )


if __name__ == "__main__" :
    main()
