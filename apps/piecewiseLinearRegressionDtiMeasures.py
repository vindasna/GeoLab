#!/usr/bin/python3
import sys, os, subprocess, shutil

import time

import pickle

import numpy as np

import matplotlib.pyplot as plt

from sklearn import metrics
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error, r2_score

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning



from sklearn.pipeline import make_pipeline

import piecewise_regression

import pwlf

from multiprocessing import Process, Lock, Manager, Value

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

#------------------------------ Global variables ------------------------------#

measures_options = [ "FA", "MD", "OD", "ICVF", "ISOVF" ]
slope_dict_piecewise = {}
slope_dict_linear = {}
verbose = 0

# To have reproductible results
_seed = 100
np.random.seed( _seed )

standartOut = sys.stdout
standartErr = sys.stderr

processCounter = 0

onlyLinear = False
onlyPiecewise = False


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
        prog="python3 piecewiseLinearRegressionDtiMeasures.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-ukb-dti", "--ukb-dti-dir",
        type=is_dir, required=True, metavar="<path>",
        help=( "UKBiobank directory with dti measures" ) )

    required.add_argument(
        "-a", "--atlas",
        type=is_dir, required=True, metavar="<path>",
        help=( "Path to the atlas" ) )

    required.add_argument(
        "-t", "--tsv",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the .tsv file with information of the subjects in "
                                                               "ukb-dti-dir" ) )

    required.add_argument(
        "-o", "--out",
        type=str, required=True, metavar="<path>",
        help=( "Output directory where to save plots" ) )

    required.add_argument(
        "-m", "--measure",
        type=str, required=True, metavar="<string>",
        help=( "Name of the measure to plot (options : FA, MD, OD, ICVF, "
                                                                    "ISOVF)" ) )

    required.add_argument(
        "-f", "--format",
        type=str, required=True, metavar="<string>", choices=[ ".bundles",
        ".bundlesdata", ".trk", ".tck" ],
        help=( "Format of the atlas to chose among [ .bundles, .bundlesdata, "
               ".trk, .tck ]" ) )

    # Optional arguments
    parser.add_argument(
        "-bn", "--bundle-name",
        type=str, default = None,
        help=("Name of the bundle to plot, it has to be the name of a "
              "bundle in the atlas directory withouth extension (default : "
              "save plot for all bundles)"))

    parser.add_argument(
        "-ol", "--only-linear",
        action='store_true', default=False,
        help=("Only do linear regression and not piecewise)"))

    parser.add_argument(
        "-opw", "--only-piecewise",
        action='store_true', default=False,
        help=("Only do piecewise regression and not linear)"))

    parser.add_argument(
        "-nbThreads", "--nbThreads",
        type=int, default = 4,
        help=("Chose the number of process to launch at the same time with "
              "multiprocessing (default : 4)"))

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
                        # Age
                        data_dict[ bundleName ][ 0 ].append( int(
                                                            infoSubject[ 0 ] ) )
                        # Measure value
                        data_dict[ bundleName ][ 1 ].append( float(
                                                            infoSubject[ 1 ] ) )
                        # Subject
                        data_dict[ bundleName ][ 2 ].append(
                                          infoSubject[ 2 ].replace( "\n", "" ) )
                        # Increment number of  measures/subjects for this bundle
                        data_dict[ bundleName ][ 3 ] += 1
                    except :
                        print( f"|{line}|" )
                        sys.exit()
            i += 1

    return( data_dict )


def

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


def getPValuesPiecewiseRegression( pw_fit ) :
    _results = pw_fit.get_results()
    p_const = _results[ "estimates" ][ "const" ][ "p_t" ]
    p_alpha1 = _results[ "estimates" ][ "alpha1" ][ "p_t" ]
    p_beta1 = _results[ "estimates" ][ "beta1" ][ "p_t" ]
    p_breakpoint = _results[ "estimates" ][ "breakpoint1" ][ "p_t" ]
    p_alpha2 = _results[ "estimates" ][ "alpha2" ][ "p_t" ]
    return( [ p_const, p_alpha1, p_beta1, p_breakpoint, p_alpha2 ] )

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


def polynomialRegression( data_dict, bundle, output_dir, measure, logFilePath,
                                                                        lock ) :
    global verbose, processCounter, seed, slope_dict_linear, \
                                      measures_options, standartOut, standartErr
    if not os.path.isdir( os.path.join( output_dir, "Linear" ) ) :
        os.mkdir( os.path.join( output_dir, "Linear" ) )
    outPath = os.path.join( output_dir, "Linear", f"{bundle}.png" )

    if os.path.isfile( outPath ) :
        return


    if ( verbose < 2 ) :
        logFile = open( logFilePath, 'a' )
        sys.stderr = logFile
        sys.stdout = logFile


    age_data = data_dict[ bundle ][ 0 ]
    measure_data = data_dict[ bundle ][ 1 ]
    subject_id_data = data_dict[ bundle ][ 2 ]

    x = np.array( age_data )
    y = np.array( measure_data )
    if ( x.shape[ 0 ] != y.shape[ 0 ] ) :
        if ( verbose < 2 ) :
            sys.stderr = standartErr
            sys.stdout = standartOut
            logFile.close()
        return
    else :
        if ( x.shape[ 0 ] < 3 ) :
            if ( verbose < 2 ) :
                sys.stderr = standartErr
                sys.stdout = standartOut
                logFile.close()
            return


    quartiles_y = np.quantile( y, [ 0.0, 0.25, 0.5, 0.75, 1 ] )
    quartiles_x = np.quantile( x, [ 0.0, 0.25, 0.5, 0.75, 1 ] )
    # q25 = quartiles_x[ 1 ]
    # q25 = quartiles_y[ 1 ]
    q25 = quartiles_y[ 0 ]

    # q75 = quartiles_x[ 3 ]
    # q75 = quartiles_y[ 3 ]
    q75 = quartiles_y[ 4 ]

    X = []
    Y = []
    for i in range( y.shape[ 0 ] ) :
        tmpMeasure = y[ i ]
        tmpAge = x[ i ]
        if q25 <= tmpMeasure <= q75 :
        # if q25 <= tmpAge <= q75 :
            Y.append( tmpMeasure )
            X.append( tmpAge )
    X = np.array( X )
    Y = np.array( Y )

    if ( X.shape[ 0 ] != Y.shape[ 0 ] ) :
        if ( verbose < 2 ) :
            sys.stderr = standartErr
            sys.stdout = standartOut
            logFile.close()
        return
    else :
        if ( X.shape[ 0 ] < 3 ) :
            if ( verbose < 2 ) :
                sys.stderr = standartErr
                sys.stdout = standartOut
                logFile.close()
            return

    # include_bias = 0 -> "intercept" = 0 but here we use LinearRegression to
    # take care of this
    # poly = PolynomialFeatures( degree = 2, include_bias = False )
    poly = PolynomialFeatures( degree = 1, include_bias = False )
    poly_features = poly.fit_transform( X.reshape( -1, 1 ) )
    poly_reg_model = LinearRegression()
    poly_reg_model.fit( poly_features, Y )

    _X = np.linspace( np.min( X ), np.max( X ), 500 )
    _poly_features = poly.fit_transform( _X.reshape( -1, 1 ) )
    _y = poly_reg_model.predict( _poly_features )

    plt.plot( X, Y, "bo" )
    plt.plot( _X, _y, "r-" )
    plt.xlabel( "Age" )
    plt.ylabel( measure )
    plt.savefig( outPath, dpi = 300, bbox_inches = "tight" )
    plt.clf()

    with lock :
        # y = a * x + b
        slope_dict_linear[ bundle ] = { "a" : poly_reg_model.coef_,
                                        "b" : poly_reg_model.intercept_ }
        output_slopes = os.path.join( output_dir, f"slopesLinear.pickle" )

        tmpSlopes_dict = dict( slope_dict_linear )
        pickle.dump( tmpSlopes_dict, open( output_slopes, 'wb' ) )

    if ( verbose < 2 ) :
        sys.stderr = standartErr
        sys.stdout = standartOut
        logFile.close()


def testHeteroskedasticity( X, y ) :
    fii = LinearRegression().fit( X.reshape( -1, 1 ), y.reshape( -1, 1 ) )
    print( "Computing residuals..." )

    _y = fii.predict( X.reshape( -1, 1 ) )

    residuals = y - _y.reshape( y.shape )

    fii_residuals = LinearRegression().fit( X.reshape( -1, 1 ),
                                                    residuals.reshape( -1, 1 ) )
    _residuals = fii_residuals.predict( X.reshape( -1, 1 ) )

    plt.plot( X, residuals, "bo" )
    plt.plot( X, _residuals, "ro" )
    plt.show()
    plt.clf()
    plt.close()
    sys.exit( 1 )


def savePiecewiseRegresionPlot( piecewiseModel, X, y, outPath, measure ) :
    # _X = np.linspace( np.min( X ), np.max( X ), np.max( X ) - np.min( X ) )
    _X = np.linspace( np.min( X ), np.max( X ), 500 )
    try :
        const =  piecewiseModel.get_results()[ "estimates" ][ "const" ][
                                                                    "estimate" ]
        alpha1 =  piecewiseModel.get_results()[ "estimates" ][ "alpha1" ][
                                                                    "estimate" ]
        beta1 =  piecewiseModel.get_results()[ "estimates" ][ "beta1" ][
                                                                    "estimate" ]
        alpha2 =  piecewiseModel.get_results()[ "estimates" ][ "alpha2" ][
                                                                    "estimate" ]
        breakpoint1 =  piecewiseModel.get_results()[ "estimates" ][
                                                   "breakpoint1" ][ "estimate" ]
        _y = np.where( _X < breakpoint1, alpha1 * _X + const,
                            const + alpha1 * _X + beta1 * ( _X - breakpoint1 ) )

        # printSummaryPiecewiceRegression( piecewiseModel, X, y )

    except :
        b1 = piecewiseModel.intercepts[ 0 ]
        a1 = piecewiseModel.slopes[ 0 ]
        b2 = piecewiseModel.intercepts[ 1 ]
        a2 = piecewiseModel.slopes[ 1 ]
        breakpoint = piecewiseModel.fit_breaks[ 1 ]

        # printSummaryPwlf( piecewiseModel, X, y )

        _y = np.where( _X < breakpoint, a1 * _X + b1, a2 * _X + b2 )

    plt.plot( X, y, "bo" )
    plt.plot( _X, _y, "r-" )
    plt.xlabel( "Age" )
    plt.ylabel( measure )

    plt.savefig( outPath, dpi = 300, bbox_inches = "tight" )
    plt.clf()



def saveLinearRegressionPlot( sm_model, X, y, outPath, measure ) :
    _X = np.linspace( np.min( X ), np.max( X ), np.max( X ) - np.min( X ) )
    try :
        _y = sm_model.predict( _X )
    except :
        _y = sm_model.predict( _X.reshape( -1, 1 ) )


    plt.plot( X, y, "bo" )
    plt.plot( _X, _y, "r-" )
    plt.xlabel( "Age" )
    plt.ylabel( measure )

    plt.savefig( outPath, dpi = 300, bbox_inches = "tight" )
    plt.clf()

def printSummaryPiecewiceRegression( pw_fit, X, y ) :
    const =  pw_fit.get_results()[ "estimates" ][ "const" ][ "estimate" ]
    alpha1 =  pw_fit.get_results()[ "estimates" ][ "alpha1" ][ "estimate" ]
    beta1 =  pw_fit.get_results()[ "estimates" ][ "beta1" ][ "estimate" ]
    alpha2 =  pw_fit.get_results()[ "estimates" ][ "alpha2" ][ "estimate" ]
    breakpoint1 =  pw_fit.get_results()[ "estimates" ][ "breakpoint1" ][
                                                                "estimate" ]
    y_fit = np.where( X < breakpoint1, alpha1 * X + const,
                     const + alpha1 * X + beta1 * ( X - breakpoint1 ) )
    mse = mean_squared_error( y, y_fit )
    print( f"Mean squared error : {mse}" )
    pw_fit.summary()

def printSummaryPwlf( pwlf_model, X, y ) :
    # predict for the determined points
    y_fit = pwlf_model.predict(X)
    mse = mean_squared_error( y, y_fit )
    print( f"Mean squared error : {mse}" )

    parameters = np.concatenate((pwlf_model.beta, pwlf_model.fit_breaks[1:-1]))

    print( f"R2 : {pwlf_model.r_squared()}" )

    # # p-values and standard error
    p = pwlf_model.p_values(method='non-linear', step_size=1e-4)
    se = pwlf_model.se  # standard errors
    parameters = np.concatenate((pwlf_model.beta,  pwlf_model.fit_breaks[1:-1]))

    header = ['Parameter type', 'Parameter value', 'Standard error', 't',
              'P > np.abs(t) (p-value)']
    print(*header, sep=' | ')
    values = np.zeros((parameters.size, 5), dtype=np.object_)
    values[:, 1] = np.around(parameters, decimals=5)
    values[:, 2] = np.around(se, decimals=5)
    values[:, 3] = np.around(parameters / se, decimals=5)
    values[:, 4] = np.around(p, decimals=5)

    for i, row in enumerate(values):
        if i < pwlf_model.beta.size:
            row[0] = 'Beta'
            print(*row, sep=' | ')
        else:
            row[0] = 'Breakpoint'
            print(*row, sep=' | ')

def getSlopesAndInterceptsPiecewise( piecewiseModel ) :
    const =  piecewiseModel.get_results()[ "estimates" ][ "const" ][
                                                                "estimate" ]
    alpha1 =  piecewiseModel.get_results()[ "estimates" ][ "alpha1" ][
                                                                "estimate" ]
    beta1 =  piecewiseModel.get_results()[ "estimates" ][ "beta1" ][
                                                                "estimate" ]
    alpha2 =  piecewiseModel.get_results()[ "estimates" ][ "alpha2" ][
                                                                "estimate" ]
    breakpoint =  piecewiseModel.get_results()[ "estimates" ][
                                               "breakpoint1" ][ "estimate" ]

    b1 = const
    a1 =  alpha1
    b2 = const + ( alpha1 - alpha2 ) * breakpoint
    a2 = alpha2

    return( b1, a1, b2, a2, breakpoint )


def computeLinearRegresionWithWeights( X, Y, weights ) :
    regr = LinearRegression()
    regr.fit( X.reshape(-1, 1), Y, weights )
    return( regr )

def saveSlopes( slope_dict, outpath ) :
    tmpSlopes_dict = dict( slope_dict )
    pickle.dump( tmpSlopes_dict, open( outpath, 'wb' ) )


def mixedLiearModel( ) :
    pass

def piecewiseRegressionPerBundle( data_dict, bundle, output_dir, measure,
                                                           logFilePath, lock ) :
    global verbose, processCounter, seed, slope_dict_piecewise, \
                          slope_dict_linear, measures_options, standartOut, \
                                          standartErr, onlyLinear, onlyPiecewise
    output_slopes_linear = os.path.join( output_dir, f"slopesLinear.pickle" )
    output_slopes_piecewise = os.path.join( output_dir,
                                                     f"slopesPiecewise.pickle" )

    if not onlyLinear :
        piecewiseDir = os.path.join( output_dir, "PieceWise" )
        if not os.path.isdir( piecewiseDir ) :
            os.mkdir( piecewiseDir )
        outPiecewisePath = os.path.join( piecewiseDir, f"{bundle}.png" )
        if os.path.isfile( outPiecewisePath ) and onlyPiecewise :
            return

    if not onlyPiecewise :
        linearDir = os.path.join( output_dir, "Linear" )
        if not os.path.isdir( linearDir ) :
            os.mkdir( linearDir )
        outLinearPath = os.path.join( linearDir, f"{bundle}.png" )
        if os.path.isfile( outLinearPath ) and onlyLinear :
            return

    if ( verbose < 2 ) :
        logFile = open( logFilePath, 'a' )
        sys.stderr = logFile
        sys.stdout = logFile

    age_data = data_dict[ bundle ][ 0 ]
    measure_data = data_dict[ bundle ][ 1 ]
    subject_id_data = data_dict[ bundle ][ 2 ]

    x = np.array( age_data )
    y = np.array( measure_data )
    if ( x.shape[ 0 ] != y.shape[ 0 ] ) :
        if not onlyPiecewise :
            with lock :
                slope_dict_linear[ bundle ] = -1
                saveSlopes( slope_dict_linear, output_slopes_linear )
        if not onlyLinear :
            with lock :
                slope_dict_piecewise[ bundle ] = -1
                saveSlopes( slope_dict_piecewise, output_slopes_piecewise )
        if ( verbose < 2 ) :
            sys.stderr = standartErr
            sys.stdout = standartOut
            logFile.close()
        return
    else :
        if ( x.shape[ 0 ] < 3 ) :
            if not onlyPiecewise :
                with lock :
                    slope_dict_linear[ bundle ] = -1
                    saveSlopes( slope_dict_linear, output_slopes_linear )
            if not onlyLinear :
                with lock :
                    slope_dict_piecewise[ bundle ] = -1
                    saveSlopes( slope_dict_piecewise, output_slopes_piecewise )
            if ( verbose < 2 ) :
                sys.stderr = standartErr
                sys.stdout = standartOut
                logFile.close()
            return

    # quartiles_y = np.quantile( y, [ 0.0, 0.25, 0.5, 0.75, 1 ] )
    # quartiles_x = np.quantile( x, [ 0.0, 0.25, 0.5, 0.75, 1 ] )
    quartiles_y = np.quantile( y, [ 0.0, 0.25, 0.5, 0.75, 1 ] )
    quartiles_x = np.quantile( x, [ 0.0, 0.25, 0.5, 0.75, 1 ] )
    # q25 = quartiles_x[ 1 ]
    # q25 = quartiles_y[ 1 ]
    q25 = quartiles_y[ 0 ]

    # q75 = quartiles_x[ 3 ]
    # q75 = quartiles_y[ 3 ]
    q75 = quartiles_y[ 4 ]

    X = []
    Y = []
    for i in range( y.shape[ 0 ] ) :
        tmpMeasure = y[ i ]
        tmpAge = x[ i ]
        # if q25 <= tmpAge <= q75 :
        if q25 <= tmpMeasure <= q75 :
            Y.append( tmpMeasure )
            X.append( tmpAge )
    X = np.array( X )
    Y = np.array( Y )

    if ( X.shape[ 0 ] != Y.shape[ 0 ] ) :
        if not onlyPiecewise :
            with lock :
                slope_dict_linear[ bundle ] = -1
                saveSlopes( slope_dict_linear, output_slopes_linear )
        if not onlyLinear :
            with lock :
                slope_dict_piecewise[ bundle ] = -1
                saveSlopes( slope_dict_piecewise, output_slopes_piecewise )
        if ( verbose < 2 ) :
            sys.stderr = standartErr
            sys.stdout = standartOut
            logFile.close()
        return
    else :
        if ( X.shape[ 0 ] < 3 ) :
            if not onlyPiecewise :
                with lock :
                    slope_dict_linear[ bundle ] = -1
                    saveSlopes( slope_dict_linear, output_slopes_linear )
            if not onlyLinear :
                with lock :
                    slope_dict_piecewise[ bundle ] = -1
                    saveSlopes( slope_dict_piecewise, output_slopes_piecewise )
            if ( verbose < 2 ) :
                sys.stderr = standartErr
                sys.stdout = standartOut
                logFile.close()
            return

    nbSubjectsByAge_dict = {}
    for tmpAge in X :
        if tmpAge in nbSubjectsByAge_dict.keys() :
            nbSubjectsByAge_dict[ tmpAge ] += 1
        else :
            nbSubjectsByAge_dict[ tmpAge ] = 1

    weights = []
    for tmpAge in X :
        # weights.append( nbSubjectsByAge_dict[ tmpAge ] )
        weights.append( 1 )
    weights = np.array( weights )
    # weights = 1 / ( weights / weights.max() )
    # weights = ( weights / weights.max() )
    weights = ( weights / np.sum( weights ) )

    if ( X.shape[ 0 ] != weights.shape[ 0 ] ) :
        if not onlyPiecewise :
            with lock :
                slope_dict_linear[ bundle ] = -1
                saveSlopes( slope_dict_linear, output_slopes_linear )
        if not onlyLinear :
            with lock :
                slope_dict_piecewise[ bundle ] = -1
                saveSlopes( slope_dict_piecewise, output_slopes_piecewise )
        if ( verbose < 2 ) :
            sys.stderr = standartErr
            sys.stdout = standartOut
            logFile.close()
        return

    if X.shape[ 0 ] < 500 :
        if not onlyPiecewise :
            with lock :
                slope_dict_linear[ bundle ] = -1
                saveSlopes( slope_dict_linear, output_slopes_linear )
        if not onlyLinear :
            with lock :
                slope_dict_piecewise[ bundle ] = -1
                saveSlopes( slope_dict_piecewise, output_slopes_piecewise )
        if ( verbose < 2 ) :
            sys.stderr = standartErr
            sys.stdout = standartOut
            logFile.close()
        return

    if not onlyLinear :
        if not os.path.isfile( outPiecewisePath )  :
            pw_fit = piecewise_regression.Fit( X, Y, n_breakpoints = 1,
                              max_iterations = 100, tolerance = 1e-5,
                                                                  n_boot = 100 )
            # Saving figure
            p_values = [ 1 ] * 5
            try :
                _p_values = getPValuesPiecewiseRegression( pw_fit )
                print( f"\n{bcolors.OKBLUE}Bundle : {bundle}{bcolors.ENDC}" )
                b1, a1, b2, a2, breakpoint = getSlopesAndInterceptsPiecewise(
                                                                        pw_fit )
                savePiecewiseRegresionPlot( pw_fit, x, y, outPiecewisePath,
                                                                       measure )
                p_values = _p_values
            except :
                my_pwlf = pwlf.PiecewiseLinFit(X, Y, seed = _seed )
                #  my_pwlf = pwlf.PiecewiseLinFit( X, Y )
                res = my_pwlf.fit(2)
                b1 = my_pwlf.intercepts[ 0 ]
                a1 = my_pwlf.slopes[ 0 ]
                b2 = my_pwlf.intercepts[ 1 ]
                a2 = my_pwlf.slopes[ 1 ]
                breakpoint = my_pwlf.fit_breaks[ 1 ]

                _p_values = my_pwlf.p_values(method='non-linear', step_size=1e-4)

                print( f"\n{bcolors.OKBLUE}Bundle : {bundle}{bcolors.ENDC}" )
                savePiecewiseRegresionPlot( my_pwlf, x, y, outPiecewisePath,
                                                                       measure )
                p_values = _p_values


            _change_model = False
            for _p in p_values :
                if _p == "-" :
                    _p = 0
                if float( _p ) > 0.05 :
                    _change_model = True

                if _change_model :
                    my_pwlf_continous = pwlf.PiecewiseLinFit( X, Y,
                                                                  seed = _seed )
                    res = my_pwlf_continous.fit(1)
                    b1 = b2 = my_pwlf_continous.intercepts[ 0 ]
                    a1 = a2 = my_pwlf_continous.slopes[ 0 ]
                    breakpoint = np.min( X )

                    _X = np.linspace( np.min( X ), np.max( X ), 500 )
                    _y = my_pwlf_continous.predict( _X )
                    printSummaryPwlf( my_pwlf_continous, X, Y )
                    plt.plot( x, y, "bo" )
                    plt.plot( _X, _y, "r-" )
                    plt.xlabel( "Age" )
                    plt.ylabel( measure )

                    plt.savefig( outPiecewisePath, dpi = 300,
                                                         bbox_inches = "tight" )
                    plt.clf()
                    break

            with lock :
                slope_dict_piecewise[ bundle ] = { "b1" : b1, "a1" : a1,
                               "b2" : b2, "a2" : a2, "breakpoint" : breakpoint }
                saveSlopes( slope_dict_piecewise, output_slopes_piecewise )


        if not onlyPiecewise :
            if not os.path.isfile( outLinearPath )  :
                my_pwlf_continous = pwlf.PiecewiseLinFit(X, Y, seed = _seed )
                res = my_pwlf_continous.fit(1)
                b1 = b2 = my_pwlf_continous.intercepts[ 0 ]
                a1 = a2 = my_pwlf_continous.slopes[ 0 ]
                breakpoint = np.min( X )

                _X = np.linspace( np.min( X ), np.max( X ), 500 )
                _y = my_pwlf_continous.predict( _X )
                printSummaryPwlf( my_pwlf_continous, X, Y )
                plt.plot( x, y, "bo" )
                plt.plot( _X, _y, "r-" )
                plt.xlabel( "Age" )
                plt.ylabel( measure )

                plt.savefig( outLinearPath, dpi = 300, bbox_inches = "tight" )
                plt.clf()

                with lock :
                    slope_dict_linear[ bundle ] = { "b1" : b1, "a1" : a1,
                               "b2" : b2, "a2" : a2, "breakpoint" : breakpoint }
                    saveSlopes( slope_dict_linear, output_slopes_linear )

        if ( verbose < 2 ) :
            sys.stderr = standartErr
            sys.stdout = standartOut
            logFile.close()

        return

    else :
        # my_pwlf_continous = pwlf.PiecewiseLinFit(X, Y, seed = _seed )
        # res = my_pwlf_continous.fit(1)
        # b1 = b2 = my_pwlf_continous.intercepts[ 0 ]
        # a1 = a2 = my_pwlf_continous.slopes[ 0 ]
        # breakpoint = np.min( X )
        #
        # _X = np.linspace( np.min( X ), np.max( X ), 500 )
        # _y = my_pwlf_continous.predict( _X )
        # printSummaryPwlf( my_pwlf_continous, X, Y )

        linearFittedModel = computeLinearRegresionWithWeights( X, Y, weights )
        _X = np.linspace( np.min( X ), np.max( X ), 500 )
        _y = linearFittedModel.predict( _X.reshape( -1, 1 ) )

        b1 = b2 = linearFittedModel.intercept_
        a1 = a2 = linearFittedModel.coef_[ 0 ]
        breakpoint = np.min( X )


        plt.plot( x, y, "bo" )
        plt.plot( _X, _y, "r-" )
        plt.xlabel( "Age" )
        plt.ylabel( measure )

        plt.savefig( outLinearPath, dpi = 300, bbox_inches = "tight" )
        plt.clf()

        with lock :
            slope_dict_linear[ bundle ] = { "b1" : b1, "a1" : a1, "b2" : b2,
                                          "a2" : a2, "breakpoint" : breakpoint }
            saveSlopes( slope_dict_linear, output_slopes_linear )

        if ( verbose < 2 ) :
            sys.stderr = standartErr
            sys.stdout = standartOut
            logFile.close()

        return


def readSlopes( path ) :
    with open( path, "rb" ) as f :
        slope_dict = pickle.load( f )
    return( slope_dict )



def main() :
    global verbose, processCounter, slope_dict_piecewise, slope_dict_linear, \
           measures_options, standartOut, standartErr, onlyLinear, onlyPiecewise
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    ukb_dti_dir = inputs[ "ukb_dti_dir" ]

    atlas_dir = inputs[ "atlas" ]

    tsv_file = inputs[ "tsv" ]

    output_dir = inputs[ "out" ]
    if not os.path.isdir( output_dir ) :
        os.mkdir( output_dir )

    measure = inputs[ "measure" ]
    if measure not in measures_options :
        print( "ERROR : the input of -m/--measure must be FA, MD, OD, ICVF or "
                                                                       "ISOVF" )
        sys.exit( 1 )


    format = inputs[ "format" ]

    nbThreads = inputs[ "nbThreads" ]
    if ( nbThreads < 1 ) :
        print( f"ERROR : argument -nbThreads must be greater than 0 but got "
               f"{nbThreads}" )
        exit( 1 )

    bundleName = inputs[ "bundle_name" ]

    onlyLinear = inputs[ "only_linear" ]

    onlyPiecewise = inputs[ "only_piecewise" ]


    ############################################################################
    lock = Lock() # Create lock for critical regions for multiprocessing
    manager = Manager() # Multiprocessing shared variable
    # processCounter = Value( 'i', 0 ) # Multiprocessing shared variable
    ############################################################################
    logFilePath = os.path.join( output_dir, "log.txt" )
    if os.path.isfile( logFilePath ) :
        os.remove( logFilePath )

    subjects = os.listdir( ukb_dti_dir )
    nbSubjects = len( subjects )

    tmp = os.listdir( atlas_dir )
    bundles_options = []
    for _tmpFile in tmp :
        if _tmpFile.endswith( format ) :
            bundles_options.append( _tmpFile.replace( format, "" ) )

    tmp = None # Free memory

    if ( len( bundles_options ) == 0 ) :
        print( f"ERROR : no bundles with format {format} where found in "
                                                                f"{atlas_dir}" )
        exit( 1 )


    out_data_dict_path = os.path.join( output_dir, f"{measure}AndAge.txt" )
    print( f"Fusing DTI mean of {measure} in {out_data_dict_path}" )

    getFileForPlotMiguelCommand = [ "getFileForPlotsRegression",
                                    "-ukb-dti", ukb_dti_dir,
                                    "-m", measure,
                                    "-tsv",  tsv_file,
                                    "-o", out_data_dict_path ]
    if not os.path.isfile( out_data_dict_path ) :
        run_sh_process( getFileForPlotMiguelCommand, shell = True )

    data_dict = readFileWithMeansMeasure( out_data_dict_path )
    nbBundles = len( data_dict )

    count_bundle = 1
    print( "Computing plot and models..." )
    maxAbsSlope = 0
    maxAbsSlopeBundle = ""
    nbBundlesWithDti = 0
    bundlesWithSignificantAgeSexInteraction = []


    output_slopes_piecewise = os.path.join( output_dir,
                                                     f"slopesPiecewise.pickle" )
    if os.path.isfile( output_slopes_piecewise ) :
        slope_dict_piecewise = manager.dict( readSlopes(
                                                     output_slopes_piecewise ) )
    else :
        slope_dict_piecewise = manager.dict()

    output_slopes_linear = os.path.join( output_dir, f"slopesLinear.pickle" )
    if os.path.isfile( output_slopes_linear ) :
        slope_dict_linear = manager.dict( readSlopes( output_slopes_linear ) )
    else :
        slope_dict_linear = manager.dict()

    if ( bundleName ) :
        if ( not bundleName in data_dict.keys() ) :
            print( f"ERROR : argument of --bundle-name -> {bundleName} is not "
                   "in keys of data_dict" )
            exit( 1 )
        piecewiseRegressionPerBundle( data_dict, bundleName, output_dir,
                                                    measure, logFilePath, lock )
        return

    processList = []
    # for bundle in data_dict :
    for bundle in bundles_options :
        p1 = Process( target = piecewiseRegressionPerBundle, args = [ data_dict,
               bundle, output_dir, measure, logFilePath, lock ], daemon = True )
        p1.start()
        processList.append( p1 )

        # p2 = Process( target = polynomialRegression, args = [ data_dict,
        #        bundle, output_dir, measure, logFilePath, lock ], daemon = True )
        # p2.start()
        # processList.append( p2 )

        printMessage = True
        while ( len( processList ) >= nbThreads ) :
            if ( printMessage ) :
                print( f"Waiting for ressources, max number of threads reached "
                       f": {len( processList )} | Bundle : [{count_bundle}/"
                       f"{nbBundles}]", end = "\r", flush = True )
                printMessage = False

            processIndex = 0
            for processTmp in processList :
                if not processTmp.is_alive() :
                    processTmp.terminate()
                    del processList[ processIndex ]
                    # processIndex is not increase because we removed an element
                else :
                    processIndex += 1

        print( f"Bundle : [{count_bundle}/{nbBundles}]" + " " * 75, end = "\r" )


        count_bundle += 1

    for processTmp in processList:
        processTmp.join()


    if not onlyLinear :
        tmpSlopes_dict_piecewise = dict( slope_dict_piecewise )
        pickle.dump( tmpSlopes_dict_piecewise, open( output_slopes_piecewise,
                                                                        'wb' ) )

    if not onlyPiecewise :
        tmpSlopes_dict_linear = dict( slope_dict_linear )
        pickle.dump( tmpSlopes_dict_linear, open( output_slopes_linear, 'wb' ) )


    print( "\nDone" )


if __name__ == "__main__" :
    t1 = time.time()
    main()
    duration = time.time() - t1
    print( f"Duration : {duration} s" )
