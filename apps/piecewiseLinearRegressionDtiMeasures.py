#!/usr/bin/python3
import sys, os, subprocess, shutil

import pickle

import numpy as np

import matplotlib.pyplot as plt

from sklearn import metrics
from sklearn.linear_model import LinearRegression, Ridge
# from sklearn.preprocessing import PolynomialFeatures, SplineTransformer
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error, r2_score
import statsmodels.api as sm
import statsmodels.formula.api as smf


from sklearn.pipeline import make_pipeline

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

# bundles_options = ['atlas_lh_IP-SP_1', 'atlas_rh_Op-SF_0', 'atlas_lh_Tr-Ins_0',
# 'atlas_lh_Op-SF_0', 'atlas_lh_PoC-PrC_0', 'atlas_lh_PrC-Ins_0',
# 'atlas_rh_RMF-SF_0', 'atlas_rh_IT-MT_2', 'atlas_lh_PoCi-SF_0',
# 'atlas_rh_PrC-SP_0', 'atlas_rh_PoC-SM_0', 'atlas_rh_PoC-SP_1',
# 'atlas_rh_Fu-LO_1', 'atlas_rh_Cu-Li_0', 'atlas_lh_MOF-ST_0',
# 'atlas_rh_MOF-ST_0', 'atlas_rh_CMF-SF_1', 'atlas_lh_ST-Ins_0',
# 'atlas_lh_PoC-PrC_3', 'atlas_lh_IP-LO_1', 'atlas_lh_Or-Ins_0',
# 'atlas_lh_Op-Ins_0', 'atlas_rh_SP-SM_0', 'atlas_lh_IP-IT_0',
# 'atlas_lh_Tr-SF_0', 'atlas_rh_SM-Ins_0', 'atlas_lh_CMF-PrC_0',
# 'atlas_lh_IP-SM_0', 'atlas_rh_Or-Ins_0', 'atlas_rh_RMF-SF_1',
# 'atlas_rh_Op-PrC_0', 'atlas_rh_PrC-SM_0','atlas_rh_PoCi-PrCu_2',
# 'atlas_lh_PoC-PrC_1','atlas_lh_LOF-Or_0', 'atlas_lh_RMF-SF_1',
# 'atlas_lh_CMF-SF_0', 'atlas_rh_MT-SM_0','atlas_lh_RAC-SF_1',
# 'atlas_rh_IP-SM_0','atlas_rh_MT-ST_0', 'atlas_rh_PoC-SP_0',
# 'atlas_rh_IP-IT_0', 'atlas_rh_PoC-PrC_2','atlas_lh_PrC-SF_0',
# 'atlas_lh_LOF-RMF_0','atlas_rh_IP-LO_0', 'atlas_lh_PrC-SM_0',
# 'atlas_lh_PoC-Ins_0', 'atlas_rh_CMF-SF_0','atlas_rh_Tr-Ins_0',
# 'atlas_lh_PoC-SM_0','atlas_rh_IC-PrCu_0', 'atlas_lh_MT-SM_0',
# 'atlas_lh_IT-MT_0', 'atlas_rh_LOF-RMF_1','atlas_lh_IP-SP_0',
# 'atlas_lh_ST-TT_0','atlas_lh_LOF-ST_0', 'atlas_lh_Op-PrC_0',
# 'atlas_lh_CMF-Op_0', 'atlas_lh_PoCi-PrCu_1','atlas_rh_RAC-SF_0',
# 'atlas_rh_Op-Tr_0','atlas_lh_SP-SM_0', 'atlas_lh_RMF-SF_0',
# 'atlas_lh_PoCi-RAC_0', 'atlas_rh_Tr-SF_0','atlas_lh_IP-MT_0',
# 'atlas_rh_LOF-MOF_0','atlas_lh_CMF-PrC_1', 'atlas_rh_IP-MT_0',
# 'atlas_rh_CMF-RMF_0', 'atlas_rh_CMF-PrC_1','atlas_rh_PoCi-PrCu_1',
# 'atlas_lh_LOF-RMF_1','atlas_rh_LOF-ST_0', 'atlas_rh_PrC-Ins_0',
# 'atlas_lh_CMF-RMF_0', 'atlas_lh_SM-Ins_0','atlas_rh_ST-TT_0',
# 'atlas_lh_PoC-SM_1','atlas_rh_CMF-PrC_0', 'atlas_lh_PoC-PrC_2',
# 'atlas_rh_CAC-PrCu_0', 'atlas_rh_PoCi-RAC_0','atlas_rh_IP-SP_0',
# 'atlas_lh_CAC-PrCu_0','atlas_rh_IT-MT_1', 'atlas_lh_CMF-PoC_0',
# 'atlas_rh_LOF-RMF_0', 'atlas_lh_Fu-LO_0','atlas_rh_Op-Ins_0',
# 'atlas_lh_IC-PrCu_0','atlas_rh_PoC-PrC_1', 'atlas_rh_CAC-PoCi_0',
# 'atlas_lh_PoCi-PrCu_0', 'atlas_rh_PoC-PrC_0','atlas_lh_MT-ST_0',
# 'atlas_rh_LO-SP_0' ]
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
        "-bc", "--baseline-characteristics",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the .tsv file with baseline characteristics of "
                                                   "subjects in ukb-dti-dir" ) )

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
        help=( "Name of the measure to plot (options : FA, MD, OD, ICVF, "
                                                                    "ISOVF)" ) )

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


def polynomialRegression( X, y ) :
    # include_bias = 0 -> "intercept" = 0 but here we use LinearRegression to
    # take care of this
    poly = PolynomialFeatures( degree = 2, include_bias = False )
    poly_features = poly.fit_transform( X.reshape( -1, 1) )
    poly_reg_model = LinearRegression()
    poly_reg_model.fit( poly_features, y )

    _X = np.linspace( np.min( X ), np.max( X ), 500 )
    _poly_features = poly.fit_transform( _X.reshape( -1, 1 ) )
    _y = poly_reg_model.predict( _poly_features )

    plt.plot( X, y, "bo" )
    plt.plot( _X, _y, "r-" )
    plt.show()
    plt.clf()
    plt.close()
    sys.exit( 1 )

# def bSplineRegression( X, y ) :
#     # Knots ~ breakpoints
#     # model = make_pipeline(SplineTransformer(n_knots = 3, degree = 1 ),
#     #                                                      Ridge( alpha = 1e-3 ) )
#     model = make_pipeline(SplineTransformer(n_knots = 3, degree = 1 ),
#                                                             LinearRegression() )
#     model.fit( X.reshape( -1, 1 ), y.reshape( -1, 1 ) )
#
#     # _y = model.predict( X.reshape( -1, 1 ) )
#
#     _X = np.linspace( np.min( X ), np.max( X ), 500 )
#     _y = model.predict( _X.reshape( -1, 1 ) )
#
#
#
# #     splt = SplineTransformer( n_knots=3, degree=1 ).fit( X.reshape( -1, 1 ) )
# #     plt.plot( X, splt.transform( X.reshape( -1, 1 ) ), "o" )
# #     print( "x"*100 )
# #     print( splt.n_features_out_ )
# #     # plt.legend(axes[1].lines, [f"spline {n}" for n in range(3)])
# #     # plt.set_title("SplineTransformer")
# #
# # #   plot knots of spline
# #     knots = splt.bsplines_[0].t
# #     plt .vlines(knots[1:-1], ymin=0, ymax=0.8, linestyles="dashed")
# #     plt.show()
# #     sys.exit( 0 )
#
#     plt.plot( X, y, "bo" )
#     # plt.plot( X, _y, "r-" )
#     plt.plot( _X, _y, "r-" )
#     plt.show()
#     plt.clf()
#     plt.close()
#     sys.exit( 1 )


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

        printSummaryPiecewiceRegression( piecewiseModel, X, y )

    except :
        b1 = piecewiseModel.intercepts[ 0 ]
        a1 = piecewiseModel.slopes[ 0 ]
        b2 = piecewiseModel.intercepts[ 1 ]
        a2 = piecewiseModel.slopes[ 1 ]
        breakpoint = piecewiseModel.fit_breaks[ 1 ]

        printSummaryPwlf( piecewiseModel, X, y )

        _y = np.where( _X < breakpoint, a1 * _X + b1, a2 * _X + b2 )

    plt.plot( X, y, "bo" )
    plt.plot( _X, _y, "r-" )
    plt.xlabel( "Age" )
    plt.ylabel( measure )

    plt.savefig( outPath, dpi = 300, bbox_inches = "tight" )
    plt.clf()
    # plt.show()
    # plt.clf()
    # plt.close()
    # sys.exit( 0 )


def saveLinearRegressionPlot( sm_model, X, y, outPath, measure ) :
    _X = np.linspace( np.min( X ), np.max( X ), np.max( X ) - np.min( X ) )
    try :
        _y = sm_model.predict( _X )
        # _y = sm_model.predict( X )
    except :
        _y = sm_model.predict( _X.reshape( -1, 1 ) )
        # _y = sm_model.predict( X.reshape( -1, 1 ) )


    plt.plot( X, y, "bo" )
    plt.plot( _X, _y, "r-" )
    plt.xlabel( "Age" )
    plt.ylabel( measure )

    plt.savefig( outPath, dpi = 300, bbox_inches = "tight" )
    plt.clf()
    # plt.show()
    # sys.exit( 0 )

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

    # _X = np.linspace( np.min( X ), np.max( X ), 500 )
    # _y = np.where( _X < breakpoint1, a1 * _X + b1, a2 * _X + b2 )
    #
    # plt.plot( X, y, "bo" )
    # plt.plot( _X, _y, "r-" )
    # plt.show()
    # plt.clf()
    # plt.close()

def main() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    ukb_dti_dir = inputs[ "ukb_dti_dir" ]

    atlas_dir = inputs[ "atlas" ]

    tsv_file = inputs[ "tsv" ]

    baseline_characteristics_path = inputs[ "baseline_characteristics" ]

    output_dir = inputs[ "out" ]
    if not os.path.isdir( output_dir ) :
        os.mkdir( output_dir )

    measure = inputs[ "measure" ]
    if measure not in measures_options :
        print( "ERROR : the input of -m/--measure must be FA, MD, OD, ICVF or "
                                                                       "ISOVF" )
        sys.exit( 1 )


    format = inputs[ "format" ]


    ############################################################################
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
    run_sh_process( getFileForPlotMiguelCommand, shell = True )

    data_dict = readFileWithMeansMeasure( out_data_dict_path )
    nbBundles = len( data_dict )

    baseline_characteristics_dict = readBaselineCharacteristics(
                                        baseline_characteristics_path, verbose )

    # To have reproductible results
    _seed = 100
    np.random.seed( _seed )

    count_bundle = 1
    print( "Saving plot and models..." )
    maxAbsSlope = 0
    maxAbsSlopeBundle = ""
    nbBundlesWithDti = 0
    bundlesWithSignificantAgeSexInteraction = []
    slope_dict = {}
    for bundle in data_dict :
        if verbose == 1 :
            print( f"Bundle : [{count_bundle}/{nbBundles}]", end = "\r" )
        else :
            print( f"Bundle : [{count_bundle}/{nbBundles}]" )
            # pass

        ########################################################################
        # if bundle != "atlas_lh_CMF-Op_0" :
        # if bundle != "atlas_rh_PoCi-PrCu_1" :
        # if bundle != "atlas_lh_Fu-LO_0" :
        # if bundle != "atlas_rh_RAC-SF_0" :
        # if bundle != "atlas_rh_LOF-MOF_0" and bundle != "atlas_rh_PrC-Ins_0" and bundle != "atlas_lh_LOF-RMF_0" and bundle != "atlas_rh_PrC-Ins_0" :
        # if bundle != "atlas_rh_LOF-MOF_0" :
        #     count_bundle += 1
        #     continue
        ########################################################################

        age_data = data_dict[ bundle ][ 0 ]
        measure_data = data_dict[ bundle ][ 1 ]
        subject_id_data = data_dict[ bundle ][ 2 ]
        sex_data = [ float( baseline_characteristics_dict[ subject_id.replace(
                   "sub-", "" ) ][ "sex" ] ) for subject_id in subject_id_data ]

        X = np.array( age_data )
        y = np.array( measure_data )

        # polynomialRegression( X, y )
        # bSplineRegression( X, y )
        # testHeteroskedasticity( X, y )


        ############### Select best number of breakking points #################
        # max_breakpoints = 5
        # pw_fit_models = piecewise_regression.model_selection.ModelSelection(
        #                                 X, y, max_breakpoints = max_breakpoints,
        #                                 max_iterations = 100, tolerance = 1e-05,
        #                                 min_distance_between_breakpoints = 0.01,
        #                                 min_distance_to_edge = 0.02,
        #                                 verbose = True )
        # _i = 0
        # for _summary in pw_fit_models.model_summaries :
        #     print( f"{bcolors.OKBLUE}Model with {_i} breakpoints"
        #                                                      f"{bcolors.ENDC}" )
        #     for _key in _summary :
        #         if _key == "estimates" :
        #             try :
        #                 print( f"{_key}")
        #                 for _key2 in _summary[ _key ] :
        #                     print( f"{_key2}" )
        #                     for _key3 in _summary[ _key ][ _key2 ] :
        #                         print( f"\t\t{_key3} : {_summary[ _key ][ _key2 ][ _key3 ]}" )
        #             except :
        #                 print( f"{_key}\t{_summary[ _key ]} " )
        #     # print( _summary )
        #     _i += 1
        # # for _i in range( max_breakpoints ) :
        # #     print( f"{bcolors.OKBLUE}Model with {_i} breakpoints "
        # #                                                      f"{bcolors.ENDC}" )
        # #     _model = pw_fit_models.models[ _i ]
        # #     _model.summary()
        # sys.exit( 1 )
        ########################################################################

        pw_fit = piecewise_regression.Fit( X, y, n_breakpoints = 1,
                          max_iterations = 100, tolerance = 1e-5, n_boot = 100 )
        # Saving figure
        outPath = os.path.join( output_dir, f"{bundle}.png" )
        p_values = [ 1 ] * 5
        try :
            _p_values = getPValuesPiecewiseRegression( pw_fit )
            print( f"\n{bcolors.OKBLUE}Bundle : {bundle}{bcolors.ENDC}" )
            b1, a1, b2, a2, breakpoint = getSlopesAndInterceptsPiecewise(
                                                                        pw_fit )
            savePiecewiseRegresionPlot( pw_fit, X, y, outPath, measure )
            p_values = _p_values
            # sys.exit( 1 )
        except :
            # mod = sm.OLS( y, X )
            # my_pwlf = mod.fit()


            my_pwlf = pwlf.PiecewiseLinFit(X, y, seed = _seed )
            #  my_pwlf = pwlf.PiecewiseLinFit( X, y )
            res = my_pwlf.fit(2)
            b1 = my_pwlf.intercepts[ 0 ]
            a1 = my_pwlf.slopes[ 0 ]
            b2 = my_pwlf.intercepts[ 1 ]
            a2 = my_pwlf.slopes[ 1 ]
            breakpoint = my_pwlf.fit_breaks[ 1 ]

            _p_values = my_pwlf.p_values(method='non-linear', step_size=1e-4)

            print( f"\n{bcolors.OKBLUE}Bundle : {bundle}{bcolors.ENDC}" )
            savePiecewiseRegresionPlot( my_pwlf, X, y,outPath, measure )
            p_values = _p_values


        _change_model = False
        for _p in p_values :
            if _p == "-" :
                _p = 0
            if float( _p ) > 0.05 :
                _change_model = True

            if _change_model :
                my_pwlf_continous = pwlf.PiecewiseLinFit(X, y, seed = _seed )
                res = my_pwlf_continous.fit(1)
                b1 = b2 = my_pwlf_continous.intercepts[ 0 ]
                a1 = a2 = my_pwlf_continous.slopes[ 0 ]
                breakpoint = np.min( X )

                _X = np.linspace( np.min( X ), np.max( X ), 500 )
                _y = my_pwlf_continous.predict( _X )
                printSummaryPwlf( my_pwlf_continous, X, y )
                plt.plot( X, y, "bo" )
                plt.plot( _X, _y, "r-" )
                plt.xlabel( "Age" )
                plt.ylabel( measure )

                plt.savefig( outPath, dpi = 300, bbox_inches = "tight" )
                plt.clf()
                # plt.show()
                # plt.clf()
                # plt.close()
                break

        slope_dict[ bundle ] = { "b1" : b1, "a1" : a1, "b2" : b2, "a2" : a2,
                                                     "breakpoint" : breakpoint }

        """
        # Plot the data, fit, breakpoints and confidence intervals
        pw_fit.plot_data(color="grey", s=20)
        # Pass in standard matplotlib keywords to control any of the plots
        pw_fit.plot_fit(color="red", linewidth=4)
        pw_fit.plot_breakpoints()
        pw_fit.plot_breakpoint_confidence_intervals()
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()
        plt.clf()
        """


        count_bundle += 1


    output_slopes = os.path.join( output_dir, f"slopes.pickle" )
    pickle.dump(slope_dict, open( output_slopes, 'wb' ) )

    print( "\nDone" )


if __name__ == "__main__" :
    main()
