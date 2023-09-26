#!${PYTHON_BINARY}

import sys, os, subprocess, shutil

import numpy as np
import scipy.stats
import scipy.special
import matplotlib.pyplot as plt

import argparse
import textwrap
from argparse import RawTextHelpFormatter

import nibabel as nib

import multiprocessing

from fitter import Fitter, get_common_distributions, get_distributions

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
python3 analyseAtlasBundle.py \
-i input.ana \
-v 1 \

"""

def is_file_or_directory(filepath):
    """ Check file's existence - argparse 'type' argument.
    """
    if not os.path.isfile(filepath):
        if not os.path.isdir(filepath):
            raise argparse.ArgumentError( "File or directory does not exist: %s"
                                                                     % filepath)
    return filepath

def is_file(filepath):
    """ Check file's existence - argparse 'type' argument.
    """
    if not os.path.isfile(filepath):
        raise argparse.ArgumentError( "File does not exist: %s" % filepath )

    return filepath


def get_cmd_line_args():
    """
    Create a command line argument parser and return a dict mapping
    <argument name> -> <argument value>.
    """
    parser = argparse.ArgumentParser(
        prog="python3 analyseAtlasBundle.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group( "required arguments" )
    required.add_argument(
        "-i", "--input",
        type=is_file_or_directory, required=True, metavar="<path>",
        help="Input directoy/file with the bundles to analyse" )

    required.add_argument(
        "-f", "--format", choices=[".bundles", ".trk", ".tck"],
        type=str, required=True, metavar="<str>",
        help="Format of bundles in atlas" )

    # Optional arguments
    parser.add_argument(
        "-o", "--output",
        type=str, metavar="<path>",
        help="Output directory ( default = input_dir/Analysis ) " )
    parser.add_argument(
        "-r", "--reference",
        type=is_file_or_directory, metavar="<path>",
        help="Reference image from which extract resolution and volume size" )
    parser.add_argument(
        "-ac", "--analyse-command-file",
        type=is_file, metavar="<path>",
        help="Path to the analyseAtlasBundle command" )
    parser.add_argument(
        "-pai", "--process-atlas-information",
        type=is_file, metavar="<path>",
        help="Path to the processAtlasInformation command" )
    parser.add_argument(
        "-show", "--interactive",
        action='store_true', default=False,
        help="Show the plot instead of saving it as a .png image" )
    parser.add_argument(
        "-saveFig", "--saveFigure",
        type=str, choices=["False", "false", "True", "true" ], default="True",
        help=("If \"true\", save plots with distribution of parameters for each"
              " bundle (default = true)") )
    parser.add_argument(
        "-useMDF", "--useMDF",
        action='store_true', default=False,
        help="Use MDF distance for analysis" )
    parser.add_argument(
        "-useMeanMDAD", "--useMeanMDAD",
        action='store_true', default=False,
        help="Use mean instead of max for MDAD" )
    parser.add_argument(
        "-force", "--overwrite",
        action='store_true', default=False,
        help="Overwrite output file even if it already exists" )
    parser.add_argument(
        "-parallel", "--multiprocessing",
        action='store_true', default=False,
        help="Use multiprocessing for computations" )
    parser.add_argument(
        "-v", "--verbose",
        type=int, choices=[0, 1, 2], default=1,
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
        print( f"{bcolors.OKBLUE}Output:{bcolors.ENDC} \n{output_cmd}\n"
               f"{bcolors.OKBLUE}Error:{bcolors.ENDC} \n{error_cmd}" )
        sys.exit( 1 )
        return( 0 )
    return( output_cmd )

def read_bundle( bundle_filename, verbose ) :

    if ( verbose > 1 ):
        print ( f'Reading {bundle_filename} ...' )

    # Checking extension bundle_filename
    if ( ( not bundle_filename.endswith( ".bundlesdata" ) ) and
                              ( not bundle_filename.endswith( ".bundles" ) ) ) :
        print( "ERROR : Wrong format for the input of readBundlesFile, got "
              f"{bundle_filename} which is not .bundles/bundlesdata " )
        sys.exit()
    else:
        if bundle_filename.endswith( ".bundlesdata" ):
            bundle_filename = bundle_filename.replace( ".bundlesdata",
                                                                    ".bundles" )

    attributes = [ 'binary',
                   'bundles',
                   'byte_order',
                   'curves_count',
                   'radio',
                   'length',
                   'nSubjects',
                   'data_file_name',
                   'format',
                   'io_mode',
                   'item_count',
                   'label_type',
                   'labels',
                   'object_type',
                   'averageRadius',
                   'minRadius',
                   'maxRadius',
                   'averageAngle',
                   'minAngle',
                   'maxAngle',
                   'averageDirectionAngle',
                   'minDirectionAngle',
                   'maxDirectionAngle',
                   'averageShapeAngle',
                   'minShapeAngle',
                   'maxShapeAngle',
                   'averageLength',
                   'minLength',
                   'maxLength',
                   'averageDisimilarity',
                   'minDisimilarity',
                   'maxDisimilarity',
                   'density',
                   'centerBundleX',
                   'centerBundleY',
                   'centerBundleZ',
                   'resolutionX',
                   'resolutionY',
                   'resolutionZ',
                   'sizeX',
                   'sizeY',
                   'sizeZ',
                   'space_dimension' ]

    minf_dict = dict()
    for i in range( len( attributes ) ) :
        minf_dict[ attributes[ i ] ] = -1


    with open( bundle_filename, "r" ) as f:
        for line in f:
            words = line.split( ":" )
            attribute_name = words[ 0 ].replace( "'", "").replace( " ", "" )
            if attribute_name in minf_dict :
                minf_dict[ attribute_name ] = words[ 1 ].replace( ",\n", "" )

    return minf_dict

def write_bundle( bundle_filename, minf_dict ) :
    if minf_dict[ "io_mode" ] != " 'binary'" :
        minf_dict[ "io_mode" ] = " 'binary'"
    if minf_dict[ "item_count" ] != " 1" :
        minf_dict[ "item_count" ] = " 1"
    if minf_dict[ "label_type" ] != " 'std_string'" :
        minf_dict[ "label_type" ] = " 'std_string'"
    if minf_dict[ "labels" ] != " [ '255' ]":
        minf_dict[ "labels" ] = " [ '255' ]"
    if minf_dict[ "object_type" ] != " 'BundleMap'":
        minf_dict[ "object_type" ] = " 'BundleMap'"

    minf = "attributes = {\n"
    for i in minf_dict :
        name = i
        try :
            value = minf_dict[ i ].replace( " ", "" )
        except :
            value = minf_dict[ i ]
        minf += f"    '{name}' : {value},\n"

    minf += "}"

    f = open( bundle_filename, 'w' )
    f.write( minf )
    f.close()


def read_minf( bundle_filename, verbose ) :
    # Checking extension bundle_filename
    if ( ( not bundle_filename.endswith( ".trk" ) ) and
                              ( not bundle_filename.endswith( ".tck" ) ) and
                              ( not bundle_filename.endswith( ".bundles" ) ) ) :
        print( "ERROR : Wrong format for the input of read_minf, got "
              f"{bundle_filename} which is not .trk/.tck/.bundles/"
               ".bundlesdata " )
        sys.exit()

    if bundle_filename.endswith( ".trk" ) :
        bundle_filename = bundle_filename.replace( ".trk", ".minf" )
    elif bundle_filename.endswith( ".tck" ) :
        bundle_filename = bundle_filename.replace( ".tck", ".minf" )
    else:
        if bundle_filename.endswith( ".bundlesdata" ):
            bundle_filename = bundle_filename.replace( ".bundlesdata",
                                                                    ".bundles" )

    if ( verbose > 1 ):
        print ( f'Reading {bundle_filename} ...' )

    attributes = [ 'averageRadius',
                   'minRadius',
                   'maxRadius',
                   'averageAngle',
                   'minAngle',
                   'maxAngle',
                   'averageDirectionAngle',
                   'minDirectionAngle',
                   'maxDirectionAngle',
                   'averageShapeAngle',
                   'minShapeAngle',
                   'maxShapeAngle',
                   'averageLength',
                   'minLength',
                   'maxLength',
                   'averageDisimilarity',
                   'minDisimilarity',
                   'maxDisimilarity',
                   'density',
                   'centerBundleX',
                   'centerBundleY',
                   'centerBundleZ',
                   'averageDistanceBetweenMedialPoints',
                   'minDistanceBetweenMedialPoints',
                   'maxDistanceBetweenMedialPoints',
                   'disimilarityWithAtlas',
                   'coverageWithAtlas',
                   'overlapWithAtlas',
                   'adjacencyWithAtlas' ]

    minf_dict = dict()
    for i in range( len( attributes ) ) :
        minf_dict[ attributes[ i ] ] = -1


    with open( bundle_filename, "r" ) as f:
        for line in f :
            words = line.split( ":" )
            attribute_name = words[ 0 ].replace( "'", "").replace( " ", "" )
            if attribute_name in minf_dict :
                minf_dict[ attribute_name ] = words[ 1 ].replace( ",\n", "" )

    return minf_dict

def write_minf( bundle_filename, minf_dict ) :
    minf = "attributes = {\n"
    for i in minf_dict :
        name = i
        try :
            value = minf_dict[ i ].replace( " ", "" )
        except :
            value = minf_dict[ i ]
        minf += f"    '{name}' : {value},\n"

    minf += "}"

    f = open( bundle_filename, 'w' )
    f.write( minf )
    f.close()

def getCurvesCount( bundle_filename, verbose ):

    if ( verbose ):
        print ( f'Reading ... ... {bundle_filename}' )

    # Checking extension bundle_filename
    if ( ( not bundle_filename.endswith( ".bundlesdata" ) ) and
                              ( not bundle_filename.endswith( ".bundles" ) ) ) :
        print( "ERROR : Wrong format for the input of readBundlesFile, got "
              f"{bundle_filename} which is not .bundles/bundlesdata " )
        sys.exit()
    else:
        if bundle_filename.endswith( ".bundlesdata" ):
            bundle_filename = bundle_filename.replace( ".bundlesdata",
                                                                    ".bundles" )

    with open( bundle_filename, "r" ) as f:
        for line in f:
            words = line.split( ":" )
            if words[ 0 ].replace( "curves_count", "" ) != words[ 0 ] :
                curves_count = words[ 1 ].replace( ",\n", "" )

    return int( curves_count )

def estimateNormalDistribution( samples, x_values, verbose = 0 ) :
    n_samples = len( samples )
    null_values = []
    for i in range( n_samples ) :
        if samples[ i ] == 0 :
            null_values.append( i )
    if null_values :
        samples = np.delete( samples, null_values )


    loc, scale = scipy.stats.norm.fit( samples )
    estimated_pdf = scipy.stats.norm.pdf( x_values, loc, scale )


    mean, variance = scipy.stats.norm.stats( loc = loc, scale = scale,
                                                            moments = 'mv' )
    median = scipy.stats.norm.median( loc = loc, scale = scale )

    confidence_interval = scipy.stats.norm.interval( 0.90, loc = loc,
                                                                 scale = scale )

    if ( verbose > 1 ) :
        print( f"mean = {mean}   |   median = {median}   |   "
                                                      f"variance = {variance}" )

    return estimated_pdf, mean, median, variance, confidence_interval

def estimateLogNormalDistribution( samples, x_values, verbose = 0 ) :
    n_samples = len( samples )
    null_values = []
    for i in range( n_samples ) :
        if samples[ i ] == 0 :
            null_values.append( i )
    if null_values :
        samples = np.delete( samples, null_values )

    sigma, loc, scale = scipy.stats.lognorm.fit( samples )
    estimated_pdf = scipy.stats.lognorm.pdf( x_values, sigma, loc, scale )


    mean, variance = scipy.stats.lognorm.stats( sigma, loc = loc, scale = scale,
                                                            moments = 'mv' )
    median = scipy.stats.lognorm.median( sigma, loc = loc, scale = scale )

    confidence_interval = scipy.stats.lognorm.interval( 0.90, sigma, loc = loc,
                                                                 scale = scale )

    if ( verbose > 1 ) :
        print( f"mean = {mean}   |   median = {median}   |   "
                                                      f"variance = {variance}" )

    return estimated_pdf, mean, median, variance, confidence_interval


def estimateGammaDistribution( samples, x_values, verbose = 0 ) :
    n_samples = len( samples )
    null_values = []
    for i in range( n_samples ) :
        if samples[ i ] == 0 :
            null_values.append( i )
    if null_values :
        samples = np.delete( samples, null_values )


    alpha, loc, scale = scipy.stats.gamma.fit( samples )
    estimated_pdf = scipy.stats.gamma.pdf( x_values, alpha, loc, scale )

    mean, variance = scipy.stats.gamma.stats( alpha , loc = loc, scale = scale,
                                                                moments = 'mv' )
    median = scipy.stats.gamma.median( alpha, loc = loc, scale = scale )

    confidence_interval = scipy.stats.gamma.interval( 0.90, alpha, loc = loc,
                                                                 scale = scale )

    if ( verbose > 1 ) :
        print( f"mean = {mean}   |   median = {median}   |   "
                                                      f"variance = {variance}" )

    return estimated_pdf, mean, median, variance, confidence_interval

def getBestFittingDistribution( samples, samplesX, distributions = None ) :
    if not distributions :
        distributions = [ 'gamma', 'lognorm', 'beta', 'burr', 'norm' ]

    fitterDistribution = Fitter( samples, distributions = distributions,
                                                                 timeout = 200 )
    fitterDistribution.fit()
    bestDistribution = fitterDistribution.get_best( method='sumsquare_error' )
    bestDistributionName = list( bestDistribution.keys() )[ 0 ]
    bestDistributionParameters = bestDistribution[ bestDistributionName ]

    if bestDistributionName == "gamma" :
        alpha = bestDistributionParameters[ 'a' ]
        loc = bestDistributionParameters[ 'loc' ]
        scale = bestDistributionParameters[ 'scale' ]
        estimated_pdf = scipy.stats.gamma.pdf( samplesX, alpha, loc, scale )
        mean, variance = scipy.stats.gamma.stats( alpha , loc = loc,
                                                 scale = scale, moments = 'mv' )
        median = scipy.stats.gamma.median( alpha, loc = loc, scale = scale )

        confidence_interval = scipy.stats.gamma.interval( 0.90, alpha,
                                                      loc = loc, scale = scale )
    elif bestDistributionName == "lognorm" :
        sigma = bestDistributionParameters[ 's' ]
        loc = bestDistributionParameters[ 'loc' ]
        scale = bestDistributionParameters[ 'scale' ]
        estimated_pdf = scipy.stats.lognorm.pdf( samplesX, sigma, loc, scale )

        mean, variance = scipy.stats.lognorm.stats( sigma, loc = loc,
                                                 scale = scale, moments = 'mv' )
        median = scipy.stats.lognorm.median( sigma, loc = loc, scale = scale )

        confidence_interval = scipy.stats.lognorm.interval( 0.90, sigma,
                                                      loc = loc, scale = scale )
    elif bestDistributionName == "beta" :
        a = bestDistributionParameters[ 'a' ]
        b = bestDistributionParameters[ 'b' ]
        loc = bestDistributionParameters[ 'loc' ]
        scale = bestDistributionParameters[ 'scale' ]
        estimated_pdf = scipy.stats.beta.pdf( samplesX, a, b, loc, scale )

        mean, variance = scipy.stats.beta.stats( a, b, loc = loc,
                                                 scale = scale, moments = 'mv' )
        median = scipy.stats.beta.median( a, b, loc = loc, scale = scale )

        confidence_interval = scipy.stats.beta.interval( 0.90, a, b,
                                                      loc = loc, scale = scale )

    elif bestDistributionName == "burr" :
        c = bestDistributionParameters[ 'c' ]
        d = bestDistributionParameters[ 'd' ]
        loc = bestDistributionParameters[ 'loc' ]
        scale = bestDistributionParameters[ 'scale' ]

        estimated_pdf = scipy.stats.burr.pdf( samplesX, c, d, loc, scale )

        mean, variance = scipy.stats.burr.stats( c, d, loc = loc,
                                                 scale = scale, moments = 'mv' )
        median = scipy.stats.burr.median( c, d, loc = loc, scale = scale )

        confidence_interval = scipy.stats.burr.interval( 0.90, c, d,
                                                      loc = loc, scale = scale )

    elif bestDistributionName == "norm" :
        loc = bestDistributionParameters[ 'loc' ]
        scale = bestDistributionParameters[ 'scale' ]
        estimated_pdf = scipy.stats.norm.pdf( samplesX, loc, scale )


        mean, variance = scipy.stats.norm.stats( loc = loc, scale = scale,
                                                                moments = 'mv' )
        median = scipy.stats.norm.median( loc = loc, scale = scale )

        confidence_interval = scipy.stats.norm.interval( 0.90, loc = loc,
                                                                 scale = scale )

    else :
        print( f"ERROR : In getBestFittingDistribution the only distribution "
               f"available are Gamma, Normal, LogNormal, beta and burr, "
               f"got {bestDistributionName}..." )
        sys.exit( 1 )

    return( estimated_pdf, mean, median, variance, confidence_interval )

# def prepareForHist( samples, sampling, distribution = None, verbose = 0 ) :
def prepareForHist( samples, sampling, distribution = None, verbose = 0 ) :
    minSamples = samples.min()
    maxSamples = samples.max()
    samplesX = np.zeros( sampling )
    samplesSizeX = ( ( maxSamples - minSamples ) / ( sampling - 1 ) )
    for i in range( sampling ):
        samplesX[ i ] = ( minSamples + samplesSizeX * i )

    try :
        mean = np.mean( samples )
        median = np.median( samples )
        variance = np.var( samples )

    except :
        mean = 0
        median = 0
        variance = 0


    ( estimated_pdf, mean_dist, median_dist, variance_dist,
         confidence_interval ) = getBestFittingDistribution( samples, samplesX )

    # if distribution == "Gamma" :
    #     ( estimated_pdf, mean_dist, median_dist, variance_dist,
    #          confidence_interval ) = ( estimateGammaDistribution( samples,
    #                                              samplesX, verbose = verbose ) )
    # elif distribution == "Normal" :
    #     ( estimated_pdf, mean_dist, median_dist, variance_dist,
    #         confidence_interval ) = ( estimateNormalDistribution( samples,
    #                                              samplesX, verbose = verbose ) )
    # elif distribution == "LogNormal" :
    #     ( estimated_pdf, mean_dist, median_dist, variance_dist,
    #      confidence_interval ) = ( estimateLogNormalDistribution( samples,
    #                                              samplesX, verbose = verbose ) )
    # else :
    #     print( f"WARNING : In prepareForHist the only distribution "
    #            f"available are Gamma, Normal and LogNormal, "
    #            f"got {distribution}... Not estimating distribution" )
    #     estimated_pdf = np.zeros( sampling )
    #     mean_dist = 0
    #     median_dist = 0
    #     variance_dist = 0
    #     confidence_interval = [ 0, 1 ]

    return ( samplesX, estimated_pdf, mean, median, variance, mean_dist,
                               median_dist, variance_dist, confidence_interval )

def plotAnalysis( inputs ):
    input = inputs[ 0 ]
    bundle_filename = inputs[ 1 ]
    inputFormat = inputs[ 2 ]
    output = inputs[ 3 ]
    reference = inputs[ 4 ]
    interactive = inputs[ 5 ]
    isSaveFigure = inputs[ 6 ]
    force = inputs[ 7 ]
    verbose = inputs[ 8 ]
    ############################################################################
    f = open( input,'rb' )
    nbFibersAtlasBundle = np.fromfile( f, dtype = np.int32, count = 1 )
    numberDisimilaritiesAtlasBundle = np.fromfile( f,
                                                   dtype = np.int32, count = 1 )
    numberDistancesBetweenMedialPoints = np.fromfile( f,
                                                   dtype = np.int32, count = 1 )

    distancesToCenterAtlasBundle = np.fromfile( f, dtype = np.float32,
                                                count = nbFibersAtlasBundle[0] )
    lengthsAtlasBundle = np.fromfile( f, dtype = np.float32,
                                                count = nbFibersAtlasBundle[0] )
    disimilaritiesAtlasBundle = np.fromfile( f, dtype = np.float32,
                                    count = numberDisimilaritiesAtlasBundle[0] )
    anglesAtlasBundle = np.fromfile( f, dtype = np.float32,
                                    count = numberDisimilaritiesAtlasBundle[0] )
    directionAnglesAtlasBundle = np.fromfile( f, dtype = np.float32,
                                    count = numberDisimilaritiesAtlasBundle[0] )
    shapeAnglesAtlasBundle = np.fromfile( f, dtype = np.float32,
                                                count = nbFibersAtlasBundle[0] )
    distancesBetweenMedialPointsBundle = np.fromfile( f, dtype = np.float32,
                                 count = numberDistancesBetweenMedialPoints[0] )
    f.close()

    # bundle_filename = input.replace( ".ana", ".bundles" )
    if ( inputFormat == ".bundles" ) :
        attributes = read_bundle( bundle_filename, verbose )
    elif ( inputFormat == ".trk" or inputFormat == ".tck" ) :
        attributes = read_minf( bundle_filename, verbose )
    else :
        print( f"ERROR : in plotAnalysis(), bundle_filename is not in "
               f".bundles, .tck or .trk format, got {bundle_filename}" )
        sys.exit( 1 )

    ########################################################################
    if ( isSaveFigure ) :
        fig, axs = plt.subplots( 3, 3 )
    sampling_1 = round( nbFibersAtlasBundle[ 0 ] / 4 )
    if sampling_1 == 1 : # To avoid division by 0 later
        sampling_1 = 2
    sampling_2 = round( numberDisimilaritiesAtlasBundle[ 0 ] / 50 )
    if sampling_2 == 1 : # To avoid division by 0 later
        sampling_2 = 2

    fontsize = 7

    if ( nbFibersAtlasBundle[ 0 ] != 0 and
                                numberDisimilaritiesAtlasBundle[ 0 ] != 0 ):
        ####################################################################
        if verbose > 1 :
            print( "\nDistance to center estimation statistics : " )

        n_samples = len( distancesToCenterAtlasBundle )
        null_values = []
        for i in range( n_samples ) :
            if distancesToCenterAtlasBundle[ i ] == 0 :
                null_values.append( i )
        if null_values :
            distancesToCenterAtlasBundle = np.delete(
                                     distancesToCenterAtlasBundle, null_values )

        distancesToCenterAtlasBundle = distancesToCenterAtlasBundle[
                                   np.isfinite( distancesToCenterAtlasBundle ) ]


        ( distancesToCenterX, estimated_pdf, mean, median, variance, mean_dist,
             median_dist, variance_dist, confidence_interval ) = prepareForHist(
                                                   distancesToCenterAtlasBundle,
                                                   sampling_1,
                                                   distribution = "Gamma",
                                                   verbose = verbose )

        confidence_interval = list( confidence_interval )
        if confidence_interval[ 0 ] < 0 :
            confidence_interval[ 0 ] = 0

        if ( isSaveFigure ) :
            n, bins, patches = axs[ 0, 0 ].hist( distancesToCenterAtlasBundle,
                                            distancesToCenterX, density = True )

            axs[ 0, 0 ].set( ylabel = 'Number of fibers' )
            axs[ 0, 0 ].set_title('Distance to center')

            axs[ 0, 0 ].text( 0.23, 0.77, f'mean = {mean:.3f}   |   '
                              f'mean_dist = {mean_dist:.3f} \n'
                              f'median = {median:.3f}   |   '
                              f'median_dist = {median_dist:.3f} \n'
                              f'variance = {variance:.3f}   |   '
                              f'variance_dist = {variance_dist:.3f}\n'
                              f'CI : [ {confidence_interval[ 0 ]:.3f}, '
                              f'{confidence_interval[ 1 ]:.3f} ]',
                              fontsize = fontsize,
                              transform = plt.gcf().transFigure )


            kde = scipy.stats.gaussian_kde( distancesToCenterAtlasBundle )
            axs[ 0, 0 ].plot( distancesToCenterX, kde( distancesToCenterX ),
                                                                          "k-" )

            if estimated_pdf.any() :
                axs[ 0, 0 ].plot( distancesToCenterX, estimated_pdf )
                axs[ 0, 0 ].axvline( confidence_interval[ 0 ], color = 'r' )
                axs[ 0, 0 ].axvline( confidence_interval[ 1 ], color = 'r' )

        attributes[ "averageRadius" ] = mean_dist
        attributes[ "minRadius" ] = confidence_interval[ 0 ]
        attributes[ "maxRadius" ] = confidence_interval[ 1 ]

        ####################################################################
        if verbose > 1 :
            print( "\nLenght estimation statistics : " )

        n_samples = len( lengthsAtlasBundle )
        null_values = []
        for i in range( n_samples ) :
            if lengthsAtlasBundle[ i ] == 0 :
                null_values.append( i )
        if null_values :
            lengthsAtlasBundle = np.delete( lengthsAtlasBundle, null_values )

        lengthsAtlasBundle = lengthsAtlasBundle[
                                             np.isfinite( lengthsAtlasBundle ) ]

        ( lengthsX, estimated_pdf, mean, median, variance, mean_dist,
             median_dist, variance_dist, confidence_interval ) = prepareForHist(
                                                    lengthsAtlasBundle,
                                                    sampling_1,
                                                    distribution = "Gamma",
                                                    verbose = verbose )

        confidence_interval = list( confidence_interval )
        if ( isSaveFigure ) :
            if confidence_interval[ 0 ] < 0 :
                confidence_interval[ 0 ] = 0

            n, bins, patches = axs[ 0, 1 ].hist( lengthsAtlasBundle, lengthsX,
                                                                density = True )

            axs[ 0, 1 ].set_title( 'Lenghts' )


            axs[ 0, 1 ].text( 0.48, 0.77, f'mean = {mean:.3f}   |   '
                              f'mean_dist = {mean_dist:.3f} \n'
                              f'median = {median:.3f}   |   '
                              f'median_dist = {median_dist:.3f} \n'
                              f'variance = {variance:.3f}   |   '
                              f'variance_dist = {variance_dist:.3f}\n'
                              f'CI : [ {confidence_interval[ 0 ]:.3f}, '
                              f'{confidence_interval[ 1 ]:.3f} ]',
                              fontsize = fontsize,
                              transform = plt.gcf().transFigure )

            kde = scipy.stats.gaussian_kde( lengthsAtlasBundle )
            axs[ 0, 1 ].plot( lengthsX, kde( lengthsX ) , "k-" )

            if estimated_pdf.any() :
                axs[ 0, 1 ].plot( lengthsX, estimated_pdf )
                axs[ 0, 1 ].axvline( confidence_interval[ 0 ], color = 'r' )
                axs[ 0, 1 ].axvline( confidence_interval[ 1 ], color = 'r' )

        attributes[ "averageLength" ] = mean_dist
        attributes[ "minLength" ] = confidence_interval[ 0 ]
        attributes[ "maxLength" ] = confidence_interval[ 1 ]

        ####################################################################
        if verbose > 1 :
            print( "\nDisimilarity estimation statistics : " )

        n_samples = len( disimilaritiesAtlasBundle )
        null_values = []
        for i in range( n_samples ) :
            if disimilaritiesAtlasBundle[ i ] == 0 :
                null_values.append( i )
        if null_values :
            disimilaritiesAtlasBundle = np.delete(
                                        disimilaritiesAtlasBundle, null_values )

        disimilaritiesAtlasBundle = disimilaritiesAtlasBundle[
                                      np.isfinite( disimilaritiesAtlasBundle ) ]

        ( disimilaritiesX, estimated_pdf, mean, median, variance, mean_dist,
             median_dist, variance_dist, confidence_interval ) = prepareForHist(
                                                      disimilaritiesAtlasBundle,
                                                      sampling_2,
                                                      distribution = "Gamma",
                                                      verbose = verbose )

        confidence_interval = list( confidence_interval )
        if confidence_interval[ 0 ] < 0 :
            confidence_interval[ 0 ] = 0

        if ( isSaveFigure ) :
            n, bins, patches = axs[ 0, 2 ].hist( disimilaritiesAtlasBundle,
                                               disimilaritiesX, density = True )

            axs[ 0, 2 ].set_title('Disimilarities')

            axs[ 0, 2 ].text( 0.77, 0.77, f'mean = {mean:.3f}   |   '
                              f'mean_dist = {mean_dist:.3f} \n'
                              f'median = {median:.3f}   |   '
                              f'median_dist = {median_dist:.3f} \n'
                              f'variance = {variance:.3f}   |   '
                              f'variance_dist = {variance_dist:.3f}\n'
                              f'CI : [ {confidence_interval[ 0 ]:.3f}, '
                              f'{confidence_interval[ 1 ]:.3f} ]',
                              fontsize = fontsize,
                              transform = plt.gcf().transFigure )


            kde = scipy.stats.gaussian_kde( disimilaritiesAtlasBundle )
            axs[ 0, 2 ].plot( disimilaritiesX, kde( disimilaritiesX ) , "k-" )

            if estimated_pdf.any() :
                axs[ 0, 2 ].plot( disimilaritiesX, estimated_pdf )
                axs[ 0, 2 ].axvline( confidence_interval[ 0 ], color = 'r' )
                axs[ 0, 2 ].axvline( confidence_interval[ 1 ], color = 'r' )

        attributes[ "averageDisimilarity" ] = mean_dist
        # attributes[ "minDisimilarity" ] = confidence_interval[ 0 ]
        attributes[ "minDisimilarity" ] = 0
        attributes[ "maxDisimilarity" ] = confidence_interval[ 1 ]

        ####################################################################
        if verbose > 1 :
            print( "\nAngle estimation statistics : " )

        n_samples = len( anglesAtlasBundle )
        null_values = []
        for i in range( n_samples ) :
            if anglesAtlasBundle[ i ] == 0 :
                null_values.append( i )
        if null_values :
            anglesAtlasBundle = np.delete( anglesAtlasBundle, null_values )

        anglesAtlasBundle = anglesAtlasBundle[
                                              np.isfinite( anglesAtlasBundle ) ]

        ( anglesX, estimated_pdf, mean, median, variance, mean_dist,
             median_dist, variance_dist, confidence_interval ) = prepareForHist(
                                                        anglesAtlasBundle,
                                                        sampling_2,
                                                        distribution = "Gamma",
                                                        verbose = verbose )

        confidence_interval = list( confidence_interval )
        if confidence_interval[ 0 ] < 0 or confidence_interval[ 0 ] > 180 :
            confidence_interval[ 0 ] = mean_dist - 3 * np.sqrt( variance_dist )
            if confidence_interval[ 0 ] < 0 :
                confidence_interval[ 0 ] = 0
        if confidence_interval[ 1 ] < 0 or confidence_interval[ 1 ] > 180 :
            confidence_interval[ 1 ] = mean_dist + 3 * np.sqrt( variance_dist )
            if confidence_interval[ 1 ] > 180 :
                confidence_interval[ 1 ] = 90

        if ( isSaveFigure ) :
            n, bins, patches = axs[ 1, 0 ].hist( anglesAtlasBundle, anglesX,
                                                                density = True )

            axs[ 1, 0 ].set_title('Angles')

            axs[ 1, 0 ].text( 0.23, 0.48, f'mean = {mean:.3f}   |   '
                              f'mean_dist = {mean_dist:.3f} \n'
                              f'median = {median:.3f}   |   '
                              f'median_dist = {median_dist:.3f} \n'
                              f'variance = {variance:.3f}   |   '
                              f'variance_dist = {variance_dist:.3f}\n'
                              f'CI : [ {confidence_interval[ 0 ]:.3f}, '
                              f'{confidence_interval[ 1 ]:.3f} ]',
                              fontsize = fontsize,
                              transform = plt.gcf().transFigure )


            kde = scipy.stats.gaussian_kde( anglesAtlasBundle )
            axs[ 1, 0 ].plot( anglesX, kde( anglesX ) , "k-" )

            if estimated_pdf.any() :
                axs[ 1, 0 ].plot( anglesX, estimated_pdf )
                axs[ 1, 0 ].axvline( confidence_interval[ 0 ], color = 'r' )
                axs[ 1, 0 ].axvline( confidence_interval[ 1 ], color = 'r' )


        attributes[ "averageAngle" ] = mean_dist
        attributes[ "minAngle" ] = confidence_interval[ 0 ]
        attributes[ "maxAngle" ] = confidence_interval[ 1 ]


        ####################################################################
        if verbose > 1 :
            print( "\nDirection angle estimation statistics : " )

        n_samples = len( directionAnglesAtlasBundle )
        null_values = []
        for i in range( n_samples ) :
            if directionAnglesAtlasBundle[ i ] <= 0 :
                null_values.append( i )
        if null_values :
            directionAnglesAtlasBundle = np.delete(
                                       directionAnglesAtlasBundle, null_values )

        directionAnglesAtlasBundle = directionAnglesAtlasBundle[
                                     np.isfinite( directionAnglesAtlasBundle ) ]

        ( directionAnglesX, estimated_pdf, mean, median, variance, mean_dist,
             median_dist, variance_dist, confidence_interval ) = prepareForHist(
                                                  directionAnglesAtlasBundle,
                                                  sampling_2,
                                                  distribution = "Gamma",
                                                  verbose = verbose )

        confidence_interval = list( confidence_interval )
        if confidence_interval[ 0 ] < 0 or confidence_interval[ 0 ] > 180 :
            confidence_interval[ 0 ] = mean_dist - 3 * np.sqrt( variance_dist )
            if confidence_interval[ 0 ] < 0 :
                confidence_interval[ 0 ] = 0
        if confidence_interval[ 1 ] < 0 or confidence_interval[ 1 ] > 180 :
            confidence_interval[ 1 ] = mean_dist + 3 * np.sqrt( variance_dist )
            if confidence_interval[ 1 ] > 180 :
                confidence_interval[ 1 ] = 180

        if ( isSaveFigure ) :
            n, bins, patches = axs[ 1, 1 ].hist( directionAnglesAtlasBundle,
                                              directionAnglesX, density = True )

            axs[ 1, 1 ].set_title( 'Direction angles' )

            axs[ 1, 1 ].text( 0.48, 0.48, f'mean = {mean:.3f}   |   '
                              f'mean_dist = {mean_dist:.3f} \n'
                              f'median = {median:.3f}   |   '
                              f'median_dist = {median_dist:.3f} \n'
                              f'variance = {variance:.3f}   |   '
                              f'variance_dist = {variance_dist:.3f}\n'
                              f'CI : [ {confidence_interval[ 0 ]:.3f}, '
                              f'{confidence_interval[ 1 ]:.3f} ]',
                              fontsize = fontsize,
                              transform = plt.gcf().transFigure )

            kde = scipy.stats.gaussian_kde( directionAnglesAtlasBundle )
            axs[ 1, 1 ].plot( anglesX, kde( anglesX ) , "k-" )

            if estimated_pdf.any() :
                axs[ 1, 1 ].plot( anglesX, estimated_pdf )
                axs[ 1, 1 ].axvline( confidence_interval[ 0 ], color = 'r' )
                axs[ 1, 1 ].axvline( confidence_interval[ 1 ], color = 'r' )

        attributes[ "averageDirectionAngle" ] = mean_dist
        attributes[ "minDirectionAngle" ] = confidence_interval[ 0 ]
        attributes[ "maxDirectionAngle" ] = confidence_interval[ 1 ]

        ####################################################################
        if verbose > 1 :
            print( "\nShape angle estimation statistics : " )

        n_samples = len( shapeAnglesAtlasBundle )
        null_values = []
        for i in range( n_samples ) :
            if shapeAnglesAtlasBundle[ i ] == 0 :
                null_values.append( i )
        if null_values :
            shapeAnglesAtlasBundle = np.delete(
                                           shapeAnglesAtlasBundle, null_values )

        shapeAnglesAtlasBundle = shapeAnglesAtlasBundle[
                                     np.isfinite( shapeAnglesAtlasBundle ) ]

        ( shapeAnglesX, estimated_pdf, mean, median, variance, mean_dist,
             median_dist, variance_dist, confidence_interval ) = prepareForHist(
                                                    shapeAnglesAtlasBundle,
                                                    sampling_1,
                                                    distribution = "Gamma",
                                                    verbose = verbose )

        confidence_interval = list( confidence_interval )
        if confidence_interval[ 0 ] < 0 or confidence_interval[ 0 ] > 180 :
            confidence_interval[ 0 ] = mean_dist - 3 * np.sqrt( variance_dist )
            if confidence_interval[ 0 ] < 0 :
                confidence_interval[ 0 ] = 0
        if confidence_interval[ 1 ] < 0 or confidence_interval[ 1 ] > 180 :
            confidence_interval[ 1 ] = mean_dist + 3 * np.sqrt( variance_dist )
            if confidence_interval[ 1 ] > 180 :
                confidence_interval[ 1 ] = 180

        if ( isSaveFigure ) :
            n, bins, patches = axs[ 1, 2 ].hist( shapeAnglesAtlasBundle,
                                                 shapeAnglesX, density = True )

            axs[ 1, 2 ].set_title( 'Shape angles' )

            axs[ 1, 2 ].text( 0.77, 0.48, f'mean = {mean:.3f}   |   '
                              f'mean_dist = {mean_dist:.3f} \n'
                              f'median = {median:.3f}   |   '
                              f'median_dist = {median_dist:.3f} \n'
                              f'variance = {variance:.3f}   |   '
                              f'variance_dist = {variance_dist:.3f}\n'
                              f'CI : [ {confidence_interval[ 0 ]:.3f}, '
                              f'{confidence_interval[ 1 ]:.3f} ]',
                              fontsize = fontsize,
                              transform = plt.gcf().transFigure )

            kde = scipy.stats.gaussian_kde( shapeAnglesAtlasBundle )
            axs[ 1, 2 ].plot( shapeAnglesX, kde( shapeAnglesX ) , "k-" )

            if estimated_pdf.any() :
                axs[ 1, 2 ].plot( shapeAnglesX, estimated_pdf )
                axs[ 1, 2 ].axvline( confidence_interval[ 0 ], color = 'r' )
                axs[ 1, 2 ].axvline( confidence_interval[ 1 ], color = 'r' )

        attributes[ "averageShapeAngle" ] = mean_dist
        attributes[ "minShapeAngle" ] = confidence_interval[ 0 ]
        attributes[ "maxShapeAngle" ] = confidence_interval[ 1 ]

        ####################################################################
        if verbose > 1 :
            print( "\nDistance between medial points estimation statistics : " )

        n_samples = len( distancesBetweenMedialPointsBundle )
        null_values = []
        for i in range( n_samples ) :
            if distancesBetweenMedialPointsBundle[ i ] == 0 :
                null_values.append( i )
        if null_values :
            distancesBetweenMedialPointsBundle = np.delete(
                               distancesBetweenMedialPointsBundle, null_values )

        distancesBetweenMedialPointsBundle = distancesBetweenMedialPointsBundle[
                             np.isfinite( distancesBetweenMedialPointsBundle ) ]

        ( distancesBetweenMedialPointsBundleX, estimated_pdf, mean, median,
                        variance, mean_dist, median_dist, variance_dist,
                        confidence_interval ) = prepareForHist(
                                             distancesBetweenMedialPointsBundle,
                                             sampling_1,
                                             distribution = "Gamma",
                                             verbose = verbose )

        confidence_interval = list( confidence_interval )
        if confidence_interval[ 0 ] < 0 or confidence_interval[ 0 ] > 180 :
            confidence_interval[ 0 ] = mean_dist - 3 * np.sqrt( variance_dist )
            if confidence_interval[ 0 ] < 0 :
                confidence_interval[ 0 ] = 0
        if confidence_interval[ 1 ] < 0 or confidence_interval[ 1 ] > 180 :
            confidence_interval[ 1 ] = mean_dist + 3 * np.sqrt( variance_dist )
            if confidence_interval[ 1 ] > 180 :
                confidence_interval[ 1 ] = 180

        if ( isSaveFigure ) :
            n, bins, patches = axs[ 2, 0 ].hist(
                                            distancesBetweenMedialPointsBundle,
                                            distancesBetweenMedialPointsBundleX,
                                                                density = True )

            axs[ 2, 0 ].set_title( 'Distance medial points' )

            axs[ 2, 0 ].text( 0.77, 0.23, f'mean = {mean:.3f}   |   '
                              f'mean_dist = {mean_dist:.3f} \n'
                              f'median = {median:.3f}   |   '
                              f'median_dist = {median_dist:.3f} \n'
                              f'variance = {variance:.3f}   |   '
                              f'variance_dist = {variance_dist:.3f}\n'
                              f'CI : [ {confidence_interval[ 0 ]:.3f}, '
                              f'{confidence_interval[ 1 ]:.3f} ]',
                              fontsize = fontsize,
                              transform = plt.gcf().transFigure )

            kde = scipy.stats.gaussian_kde( distancesBetweenMedialPointsBundle )
            axs[ 2, 0 ].plot( distancesBetweenMedialPointsBundleX,
                             kde( distancesBetweenMedialPointsBundleX ) , "k-" )

            if estimated_pdf.any() :
                axs[ 2, 0 ].plot( distancesBetweenMedialPointsBundleX,
                                                                 estimated_pdf )
                axs[ 2, 0 ].axvline( confidence_interval[ 0 ], color = 'r' )
                axs[ 2, 0 ].axvline( confidence_interval[ 1 ], color = 'r' )

        attributes[ "averageDistanceBetweenMedialPoints" ] = mean_dist
        attributes[ "minDistanceBetweenMedialPoints" ] = confidence_interval[
                                                                             0 ]
        attributes[ "maxDistanceBetweenMedialPoints" ] = confidence_interval[
                                                                             1 ]

        ####################################################################
        if ( isSaveFigure ) :
            for i in range( 3 ) :
                for j in range( 3 ):
                    axs[ i, j ].tick_params( axis = 'x', labelsize = 4 )
                    axs[ i, j ].tick_params( axis = 'y', labelsize = 4 )


        if ( inputFormat == ".bundles" ) :
            if ( ( float( attributes[ "resolutionX" ] ) == -1 or
                   float( attributes[ "resolutionY" ] ) == -1 or
                   float( attributes[ "resolutionZ" ] ) == -1 ) and reference and
                   inputFormat == ".bundles" ) :
                attributes[ "resolutionX" ] = reference[ "pixdim" ][ 1 ]
                attributes[ "resolutionY" ] = reference[ "pixdim" ][ 2 ]
                attributes[ "resolutionZ" ] = reference[ "pixdim" ][ 3 ]

            if ( ( float( attributes[ "sizeX" ] ) == -1 or
                   float( attributes[ "sizeY" ] ) == -1 or
                   float( attributes[ "sizeZ" ] ) == -1 ) and reference and
                   inputFormat == ".bundles" ) :
                attributes[ "sizeX" ] = reference[ "dim" ][ 1 ]
                attributes[ "sizeY" ] = reference[ "dim" ][ 2 ]
                attributes[ "sizeZ" ] = reference[ "dim" ][ 3 ]
            write_bundle( bundle_filename, attributes )
        elif ( inputFormat == ".trk" ) :
            bundleInfoFilename = bundle_filename.replace( ".trk", ".minf" )
            write_minf( bundleInfoFilename, attributes )
        elif ( inputFormat == ".tck" ) :
            bundleInfoFilename = bundle_filename.replace( ".tck", ".minf" )
            write_minf( bundleInfoFilename, attributes )
        else :
            print( f"ERROR : in plotAnalysis(), bundle_filename is not in "
                   f".bundles, .tck or .trk format, got {bundle_filename}" )
            sys.exit( 1 )

        if ( interactive ) :
            plt.show()
        elif ( isSaveFigure ) :
            if ( os.path.isfile( output ) and not force ):
                if verbose :
                    print( f'File {output} already exists')
            else:
                # plt.savefig( output, dpi = 1500, bbox_inches="tight" )
                plt.savefig( output, dpi = 300, bbox_inches = "tight" )


    else:
        if verbose > 1 :
            print( f'WARNING : Not enough fibers in {input} to do analysis' )

    if os.path.isfile( input ):
        os.remove( input )

    plt.close( 'all' )

def main():
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()
    input_path = inputs[ "input" ]
    output_path = inputs[ "output" ]
    inputFormat = inputs[ "format" ]
    if os.path.isdir( input_path ) :
        if not output_path :
            output_path = os.path.join( input_path, "Analysis" )
        if not os.path.isdir( output_path ) :
            os.mkdir( output_path )
    if os.path.isfile( input_path ) :
        if not input_path.endswith( inputFormat ) :
            print( f"ERROR : input file is not the same format as input format" )
        input_dir = os.path.dirname( input_path )
        output_path = os.path.join( input_dir, "Analysis" )
        if not os.path.isdir( output_path ) :
            os.mkdir( output_path )

    reference_filename = inputs[ "reference" ]
    if not ( os.path.isfile( reference_filename ) and (
                                    reference_filename.endswith( ".nii.gz") or
                                     reference_filename.endswith( ".nii" ) ) ) :
        print( f"WARNING : Problem reading reference image : {reference}, the "
                "file might not exists or it is not a nifti file. Ignoring "
                "reference image ... " )
        reference = None
    else :
        reference = nib.load( reference_filename ).header


    interactive = inputs[ "interactive" ]

    isSaveFigure = inputs[ "saveFigure" ]
    if isSaveFigure == "False" or isSaveFigure == "false" :
        isSaveFigure = False
    else :
        isSaveFigure = True

    useMDF = inputs[ "useMDF" ]

    useMeanMDAD = inputs[ "useMeanMDAD" ]

    force = inputs[ "overwrite" ]

    parallel = inputs[ "multiprocessing" ]

    analyse_atlas_bundles_command = inputs[ "analyse_command_file" ]
    if not analyse_atlas_bundles_command :
        analyse_atlas_bundles_command = ( "analyseAtlasBundle" )

    process_atlas_information_command = inputs[ "process_atlas_information" ]
    if not process_atlas_information_command :
        process_atlas_information_command = ( "processAtlasInformation" )


    ############################################################################
    tmp_dir = ""
    if os.path.isfile( input_path ) :
        tmp_dir = os.path.join( output_path, "tmpDir" )
        if not os.path.isdir( tmp_dir ) :
            os.mkdir( tmp_dir )
        shutil.copy( input_path, tmp_dir )
        if input_path.endswith( ".bundles" ) :
            shutil.copy( input_path.replace( ".bundles", ".bundlesdata" ),
                                                                       tmp_dir )
        elif input_path.endswith( ".bundlesdata" ) :
            shutil.copy( input_path.replace( ".bundlesdata", ".bundles" ),
                                                                       tmp_dir )
        elif input_path.endswith( ".trk" ) or input_path.endswith( ".tck" ) :
            shutil.copy( input_path.replace( inputFormat, ".minf" ), tmp_dir )
        input_path = tmp_dir
        output_path = tmp_dir

    bundles_path = input_path



    print( "Processing atlas information ...   ", end = "", flush = True )
    if os.path.isfile( reference_filename ) :
        process_atlas_command = [ process_atlas_information_command,
                                  "-a ", input_path,
                                  "-r ", reference_filename,
                                  "-f ", inputFormat,
                                  "-v ", str( 1 ) ]
        if useMDF :
            process_atlas_command.append( "-useMDF" )
            process_atlas_command.append( "true" )
        if useMeanMDAD :
            process_atlas_command.append( "-useMeanMDAD" )
            process_atlas_command.append( "true" )
    else :
        process_atlas_command = [ process_atlas_information_command,
                                  "-a ", input_path,
                                  "-f ", inputFormat,
                                  "-v ", str( 1 ) ]
        if useMDF :
            process_atlas_command.append( "-useMDF" )
            process_atlas_command.append( "true" )
        if useMeanMDAD :
            process_atlas_command.append( "-useMeanMDAD" )
            process_atlas_command.append( "true" )

    run_sh_process( cmd = process_atlas_command, shell = True )
    print( "Done", flush = True )

    print( "Analysing bundles ...   ", end = "", flush = True )
    analyse_command = [ analyse_atlas_bundles_command,
                            "-a ", input_path,
                            "-f ", inputFormat,
                            "-o ", output_path,
                            "-v ", str( 1 ) ]
    if useMDF :
            analyse_command.append( "-useMDF" )
            analyse_command.append( "true" )
    if useMeanMDAD :
        analyse_command.append( "-useMeanMDAD" )
        analyse_command.append( "true" )
    
    run_sh_process( cmd = analyse_command, shell = True )
    print( "Done", flush = True )

    input_filenames = []
    output_filenames = []
    bundle_filenames = []
    for filename in os.listdir( output_path ):
        input_path_file = os.path.join( output_path, filename )
        if ( filename.endswith( ".ana" ) ):
            input_filenames.append( input_path_file )
            output_path_file = os.path.join( output_path, filename.replace(
                                                          ".ana", ".png" ) )
            output_filenames.append( output_path_file )
            if os.path.isdir( bundles_path ) :
                bundle_filenames.append( os.path.join( bundles_path,
                                     filename.replace( ".ana", inputFormat ) ) )
            elif os.path.isfile( bundles_path ) :
                bundle_filenames.append( bundles_path )
            else :
                print( f"ERROR :Input {bundles_path} is not a directory nor a "
                                                                       "file " )
                sys.exit()
        else:
            if verbose > 2 :
                print( f'WARNING : Impossible to process {input_path_file},'
                        ' wrong file format (must be .ana)' )

    nb_input_filenames = len( input_filenames )
    if nb_input_filenames:
        if ( interactive and nb_input_filenames > 1 ):
            if verbose:
                print( "Interactive mode not compatible when input "
                       "is a directory with more than one file.")
            interactive = 0

        if not parallel:
            counterFileProcessed = 1
            for file_nb in range( nb_input_filenames ) :
                print( f"Processing file : [{counterFileProcessed}/"
                       f"{nb_input_filenames}]", end = "\r" )
                input = input_filenames[ file_nb ]
                bundle_filename = bundle_filenames[ file_nb ]
                output = output_filenames[ file_nb ]
                inputs = [ input, bundle_filename, inputFormat, output,
                          reference, interactive, isSaveFigure, force, verbose ]
                bundle_name = bundle_filename.replace( inputFormat, "" )
                print( f"Processing files : [ {file_nb+1} / "
                       f"{nb_input_filenames} ]   |   "
                       f"name : {os.path.basename( bundle_name ):30}",
                       end = '\r' )
                plotAnalysis( inputs )

            print( "\nDone" )
        else:
            if ( verbose > 1 ) :
                print( "Using multiprocessing" )
            print( "Processing..." )
            args = []
            for file_nb in range( nb_input_filenames ) :
                input = input_filenames[ file_nb ]
                bundle_filename = bundle_filenames[ file_nb ]
                output = output_filenames[ file_nb ]
                inputs = [ input, bundle_filename, inputFormat, output,
                          reference, interactive, isSaveFigure, force, verbose ]
                args.append( inputs )

            pool_obj = multiprocessing.Pool()
            pool_obj.map( plotAnalysis, args )

            print( "\nDone" )
    else:
        print( f'Error : no .ana files found in input directory {input_path}' )
        sys.exit(1)

    if tmp_dir :
        if os.path.isdir( tmp_dir ):
            shutil.rmtree( tmp_dir )

if __name__ == "__main__":
    main()
