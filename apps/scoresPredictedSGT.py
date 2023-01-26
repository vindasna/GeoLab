#!/usr/bin/python3

import os, sys, shutil

import numpy as np

from sklearn.metrics import f1_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import confusion_matrix


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

Plot the statistical analysis of a bundle

Command example:
----------------
python3 scorePredictionSGT.py \
-pl predictedLabels.txt \
-pd predictedLabels.dict \
-td trueLabels.dict \
-td trueLabels.dict \

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

def is_dir(dirpath):
    """ Check file's existence - argparse 'type' argument.
    """
    if not os.path.isdir(dirpath):
        raise argparse.ArgumentError( "Directory does not exist: %s" % dirpath )

    return dirpath


def get_cmd_line_args():
    """
    Create a command line argument parser and return a dict mapping
    <argument name> -> <argument value>.
    """
    parser = argparse.ArgumentParser(
        prog="python3 scorePredictionSGT.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group( "required arguments" )
    required.add_argument(
        "-pl", "--predicted-labels",
        type=is_file, required=True, metavar="<path>",
        help="Path to the predicted labels .txt (or .bin for SupWMA)" )

    required.add_argument(
        "-pd", "--predicted-dictionary",
        type=is_file, required=True, metavar="<path>",
        help=("Path to the predicted labels dictionary .dict (or .bin for "
                                                                    "SupWMA)") )

    required.add_argument(
        "-tl", "--true-labels",
        type=is_file, required=True, metavar="<path>",
        help="Path to the true labels .txt (or .bin for SupWMA)" )

    required.add_argument(
        "-td", "--true-dictionary",
        type=is_file, required=True, metavar="<path>",
        help=("Path to the true labels dictionary .dict (or .bin for SupWMA)") )

    required.add_argument(
        "-o", "--output",
        type=str, required=True, metavar="<path>",
        help="Output directory " )


    # Optional arguments
    parser.add_argument(
        "-force", "--overwrite",
        action='store_true', default=False,
        help="Overwrite output file even if it already exists" )
    parser.add_argument(
        "-supWMA", "--supWMA",
        action='store_true', default=False,
        help="Use flag if labels come from SupWMA" )
    parser.add_argument(
        "-v", "--verbose",
        type=int, choices=[0, 1, 2], default=1,
        help="Increase the verbosity level: 0 silent, 1 verbose." )


    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)
    verbose = kwargs.pop("verbose")

    return kwargs, verbose


def readLabels( path ) :
    nbFibers = 0
    with open( path, "r" ) as f:
        for line in f:
            words = line.split( ":" )
            if int( words[ 0 ] ) > nbFibers :
                nbFibers = int( words[ 0 ] )

    nbFibers += 1

    labels = {}
    with open( path, "r" ) as f:
        for line in f:
            words = line.split( ":" )
            if ( int( words[ 0 ] ) not in labels.keys() ) :
                labels[ int( words[ 0 ] ) ] = [ int( words[ 1 ] ) ]
            else :
                labels[ int( words[ 0 ] ) ].append( int( words[ 1 ] ) )

    return( labels )

def readLabelsSupWMA( path ) :
    with open( path, "rb" ) as f:
        x = np.fromfile( f, dtype = np.int16 )

    labels = {}
    for i in range( x.shape[ 0 ] ) :
        labels[ i ] = [ x[ i ] ]

    return( labels )


def readDict( path ) :
    _dict = {}
    with open( path, "r" ) as f:
        for line in f:
            words = line.split( ":" )
            _dict[ int( words[ 1 ] ) ] = words[ 0 ].replace( " ", "" )

    return( _dict )

def saveDict( inDict, path ) :
    with open( path, 'w' ) as f :
        for _key in inDict :
            f.write( f"{inDict[ _key ]} : {_key}\n" )


def printScoresPerBundles( bundlesNames, scores, scoreName = None ) :
    nbBundlesDict = len( bundlesNames )
    nbBundlesScores = len( scores )
    if nbBundlesDict != nbBundlesScores :
        print( f"ERROR in printScoresPerBundles() : budlesDict and scores must"
               f" have the same length, got {nbBundlesDict} and "
               f"{nbBundlesScores} respectivly" )
        return

    if scoreName :
        print( f"Bundle Name\t|\t{scoreName}" )
    for label in range( nbBundlesDict ) :
        bundleName = bundlesNames[ label ]
        print( f"{bundleName}\t:\t{scores[ label ]}" )

def getLabelFromBundleName( bundlesDict, bundleName ) :
    for label in bundlesDict.keys() :
        if bundlesDict[ label ] == bundleName :
            return( label )
    print( f"Error in getLabelFromBundleName() : \'{bundleName}\' not found in "
           f"bundlesDict" )
    print( f"Bundles names in bundlesDict : " )
    for label in bundlesDict.keys() :
        print( f"\'{bundlesDict[ label ]}\'" )
    sys.exit( 1 )

def saveScoresPerBundle( accuracies_per_bundle,
                         f1_score_per_bundle_model,
                         precision_score_per_bundle_model,
                         recall_score_per_bundle_model,
                         labelsNames,
                         out_scores_per_bundle_path ) :
    accLen = len( accuracies_per_bundle )
    f1Len = len( f1_score_per_bundle_model )
    preLen = len( precision_score_per_bundle_model )
    recLen = len( recall_score_per_bundle_model )
    lnLen = len( labelsNames )
    if ( accLen != f1Len ) :
        print( "ERROR in saveScoresPerBundle() : F1-scores must be the same "
              f"dimension as accuracies, got {f1Len} and {accLen} respectivly" )
        sys.exit( 1 )
    if ( accLen != preLen ) :
        print( "ERROR in saveScoresPerBundle() : precisions must be the same "
             f"dimension as accuracies, got {preLen} and {accLen} respectivly" )
        sys.exit( 1 )
    if ( accLen != recLen ) :
        print( "ERROR in saveScoresPerBundle() : recalls must be the same "
             f"dimension as accuracies, got {recLen} and {accLen} respectivly" )
        sys.exit( 1 )
    if ( accLen != lnLen ) :
        print( "ERROR in saveScoresPerBundle() : labels names must be the same "
              f"dimension as accuracies, got {lnLen} and {accLen} respectivly" )
        sys.exit( 1 )

    nbBundles = len( accuracies_per_bundle )

    with open( out_scores_per_bundle_path, 'w' ) as f :
        f.write( "BundleName\tAccuracy\tF1-score\tPrecision\tRecall\n" )
        for i in range( nbBundles ) :
            labelName = labelsNames[ i ]
            accuracy_bundle = accuracies_per_bundle[ i ]
            f1_bundle = f1_score_per_bundle_model[ i ]
            precision_bundle = precision_score_per_bundle_model[ i ]
            recall_bundle = recall_score_per_bundle_model[ i ]
            f.write( f"{labelName}\t{accuracy_bundle}\t{f1_bundle}\t"
                     f"{precision_bundle}\t{recall_bundle}\n" )

    #--------------------------------------------------------------------------#
    mean_acc = np.mean( accuracies_per_bundle )
    std_acc = np.std( accuracies_per_bundle )

    mean_f1 = np.mean( f1_score_per_bundle_model )
    std_f1 = np.std( f1_score_per_bundle_model )

    mean_pre = np.mean( precision_score_per_bundle_model )
    std_pre = np.std( precision_score_per_bundle_model )

    mean_rec = np.mean( recall_score_per_bundle_model )
    std_rec = np.std( recall_score_per_bundle_model )

    print( f"Scores per bundles ( mean +- std ) :\n"
           f" Accuracy : {mean_acc} +- {std_acc} \n"
           f" F1-Score : {mean_f1} +- {std_f1} \n"
           f" Precision : {mean_pre} +- {std_pre} \n"
           f" Recall : {mean_rec} +- {std_rec}" )




def saveScoresPerBundle2( sensitivities,
                          specificities,
                          accuracies,
                          precisions,
                          jaccards,
                          labelsNames,
                          out_scores_per_bundle_path ) :

    # Sorting dictionary in alphabetical order
    _tmpDict = {}
    for _key, _value in labelsNames.items() :
        _tmpDict[ _key ] = _value.lower()
    _tmpDict = { k: v for k, v in sorted( _tmpDict.items(),
                                                key = lambda item: item[ 1 ] ) }

    sortedLabelsNames = {}
    for _key in _tmpDict :
        sortedLabelsNames[ _key ] = labelsNames[ _key ]


    senLen = len( sensitivities )
    speLen = len( specificities )
    accLen = len( accuracies )
    preLen = len( precisions )
    jacLen = len( jaccards )
    lnLen = len( sortedLabelsNames )
    if ( senLen != speLen ) :
        print( f"ERROR in saveScoresPerBundle2() : Specificities must be the "
               f"same size as sensitivities, got {speLen} and {senLen} "
               f"respectivly" )
        sys.exit( 1 )
    if ( senLen != accLen ) :
        print( f"ERROR in saveScoresPerBundle2() : accuracies must be the same "
               f"size as sensitivities, got {accLen} and {senLen} respectivly" )
        sys.exit( 1 )
    if ( senLen != preLen ) :
        print( f"ERROR in saveScoresPerBundle2() : precisions must be the same "
               f"size as sensitivities, got {preLen} and {senLen} respectivly" )
        sys.exit( 1 )
    if ( senLen != jacLen ) :
        print( f"ERROR in saveScoresPerBundle2() : jaccards must be the same "
               f"size as sensitivities, got {jacLen} and {senLen} respectivly" )
        sys.exit( 1 )
    if ( senLen != lnLen ) :
        print( "ERROR in saveScoresPerBundle2() : labels names must be the same"
               f" size as sensitivities, got {lnLen} and {senLen} respectivly" )
        sys.exit( 1 )

    nbBundles = len( sensitivities )

    with open( out_scores_per_bundle_path, 'w' ) as f :
        f.write( "BundleName\tSensitivity\tSpecificity\tAccuracy\tPrecision\t"
                 "Jaccard\n" )
        # for i in range( nbBundles ) :
        for i in sortedLabelsNames :
            _sensitivity = sensitivities[ i ]
            _specificity = specificities[ i ]
            _accuracy = accuracies[ i ]
            _precision = precisions[ i ]
            _jaccard = jaccards[ i ]
            _lableName = sortedLabelsNames[ i ]
            f.write( f"{_lableName}\t{_sensitivity}\t{_specificity}\t"
                     f"{_accuracy}\t{_precision}\t{_jaccard}\n" )

    #--------------------------------------------------------------------------#
    mean_sen = np.mean( sensitivities )
    std_sen = np.std( sensitivities )

    mean_spe = np.mean( specificities )
    std_spe = np.std( specificities )

    mean_acc = np.mean( accuracies )
    std_acc = np.std( accuracies )

    mean_pre = np.mean( precisions )
    std_pre = np.std( precisions )

    mean_jac = np.mean( jaccards )
    std_jac = np.std( jaccards )

    print( f"Scores per bundles ( mean +- std ) :\n"
           f" Sensitivity : {mean_sen} +- {std_sen} \n"
           f" Specificity : {mean_spe} +- {std_spe} \n"
           f" Accuracy : {mean_acc} +- {std_acc} \n"
           f" Precision : {mean_pre} +- {std_pre} \n"
           f" Jaccard : {mean_jac} +- {std_jac}" )

def saveConfusionMatrix( confusion_matrix_model, path ) :
    with open( path, 'w' ) as f :
        sizeX = confusion_matrix_model.shape[ 0 ]
        sizeY = confusion_matrix_model.shape[ 1 ]
        for i in range( sizeX ) :
            for j in range( sizeY ) :
                f.write( f"{confusion_matrix_model[ i, j ]}\t" )
            f.write( "\n" )

def readConfusionMatrix( path ) :
    confusion_matrix_model = []
    with open( path, 'r' ) as f :
        for line in f :
            words = line.split( "\t" )
            lineConfusionMatrix = [ int( x ) for x in words if x != "\n" ]
            confusion_matrix_model.append( lineConfusionMatrix )

    return( np.array( confusion_matrix_model ) )


def comparePredictionToTrue( realLabelsPath, realDictPath, predictedLabelsPath,
                          predictedDictPath, isSupWMA = False, force = False ) :
    confusionMatrixSavingPath = os.path.dirname( predictedLabelsPath )
    _tmpBasename = os.path.basename( os.path.dirname( realLabelsPath ) )
    confusionMatrixSavingPath = os.path.join( confusionMatrixSavingPath,
                                         f"confusionMatrix_{_tmpBasename}.tsv" )

    confusionMatrixDictSavingPath = os.path.dirname( predictedLabelsPath )
    _tmpBasename = os.path.basename( os.path.dirname( realLabelsPath ) )
    confusionMatrixDictSavingPath = os.path.join( confusionMatrixDictSavingPath,
                                        f"confusionMatrix_{_tmpBasename}.dict" )

    out_scores_per_bundle_path = os.path.join(
                                     os.path.dirname( predictedLabelsPath ),
                                     f"scoresPerBundle_{_tmpBasename}.tsv" )

    if ( not os.path.isfile( confusionMatrixSavingPath ) or
                not os.path.isfile( confusionMatrixDictSavingPath ) or force ) :
        print( "Reading true labels... ", end = "" )
        realLabels = readLabels( realLabelsPath )
        print( "Done\nReading true labels dictionary... ", end = "" )
        realDict = readDict( realDictPath )

        print( "Done\nReading predicted labels... ", end = "" )
        if ( isSupWMA ) :
            predictedLabels = readLabelsSupWMA( predictedLabelsPath )
        else :
            predictedLabels = readLabels( predictedLabelsPath )
        print( "Done\nReading predicted labels dictionary... ", end = "" )
        predictedDict = readDict( predictedDictPath )
        print( "Done" )
        nbRealLabels = len( realLabels )
        nbPredictedLabels = len( predictedLabels )


        if nbRealLabels != nbPredictedLabels :
            print( f"ERROR : the number of predicted labels has to be the same "
                   f"as the number of true labels, got {nbRealLabels} true  "
                   f"labels and {nbPredictedLabels} predicted labels" )
            sys.exit( 1 )


        labelSupWMAUnlabelled = None
        if isSupWMA :
            labelSupWMAUnlabelled = getLabelFromBundleName( predictedDict,
                                                             "unlabeledFibers" )
            for realLabel in realDict.keys() :
                if realDict[ realLabel ] == "unlabeledFibers" :
                    labelSupWMAUnlabelled = realLabel
                    break

            if not labelSupWMAUnlabelled :
                print( f"WARNING : In labels dictionary of SupWMA did not find "
                       f"\"unlabeledFibers\" label, setting to -1" )
                labelSupWMAUnlabelled = -1
        if isSupWMA :
            print( f"Unlabeled in SupWMA dict : {labelSupWMAUnlabelled}" )


        # Putting missing labels in dictionaries
        if len( realDict.keys() ) > len( predictedDict.keys() ) :
            for _key in realDict :
                if realDict[ _key ] not in predictedDict.values() :
                    predictedDict[ len( predictedDict.keys() ) ] = realDict[
                                                                          _key ]
        else :
            for _key in predictedDict :
                if predictedDict[ _key ] == "unlabeledFibers" :
                    continue
                else :
                    if predictedDict[ _key ] not in realDict.values() :
                        realDict[ len( realDict.keys() ) ] = predictedDict[
                                                                          _key ]

        yTrue = []
        yPred = []
        for i in range( nbPredictedLabels ) :
            if i%1000 == 0 or i == nbPredictedLabels - 1 :
                print( f"Processing : [{i+1}/{nbRealLabels}]", end = "\r" )
            for tmPredictedLabel in predictedLabels[ i ] :
                if isSupWMA and tmPredictedLabel == labelSupWMAUnlabelled :
                    predictedLabel = -1
                else :
                    predictedLabel = tmPredictedLabel
                yPred.append( predictedLabel )
                isPredictedUnlabeled = False
                if predictedLabel == -1 :
                    isPredictedUnlabeled = True
                    if -1 in realLabels[ i ] :
                        yTrue.append( -1 )
                    else :
                        _tmpLabel =  getLabelFromBundleName( predictedDict,
                                              realDict[ realLabels[ i ][ 0 ] ] )
                        yTrue.append( _tmpLabel )


                if isPredictedUnlabeled :
                    continue


                isPredictedEqualsReal = False
                for realLabel in realLabels[ i ] :
                    if realLabel == -1 :
                        continue
                    if ( predictedDict[ predictedLabel ] ==
                                                       realDict[ realLabel ] ) :
                        yTrue.append( predictedLabel )
                        isPredictedEqualsReal = True
                        break
                if not isPredictedEqualsReal :
                    if realLabels[ i ][ 0 ] == -1 :
                        _tmpLabel = -1
                    else :
                        _tmpLabel =  getLabelFromBundleName( predictedDict,
                                              realDict[ realLabels[ i ][ 0 ] ] )
                    yTrue.append( _tmpLabel )


        print( "\nGetting predicted labels" )
        testPredictedLabels = np.unique( yPred )
        print( f"Number of labels in yPred : {len(testPredictedLabels)}" )

        print( "Getting real labels" )
        testRealLabels = np.unique( yTrue )
        print( f"Number of labels in yTrue : {len(testRealLabels)}" )



        confusion_matrix_model = confusion_matrix( yTrue, yPred )
        saveConfusionMatrix( confusion_matrix_model, confusionMatrixSavingPath )
        # accuracy_model = accuracy_score( yTrue, yPred )
        # print( f"Accuracy : {accuracy_model}" )
        #
        # balanced_accuracy_model = balanced_accuracy_score( yTrue, yPred )
        # print( f"Balanced accuracy : {balanced_accuracy_model}" )
        #
        # f1_score_model = f1_score( yTrue, yPred, average = "macro" )
        # # f1_score_model = f1_score( yTrue, yPred, average = "weighted" )
        # print( f"F1 Score : {f1_score_model}" )
        #
        #
        # f1_score_per_bundle_model = f1_score( yTrue, yPred, average = None )
        # precision_score_per_bundle_model = precision_score( yTrue, yPred,
        #                                                             average = None )
        # recall_score_per_bundle_model = recall_score( yTrue, yPred, average = None )
        # cm = np.copy( confusion_matrix_model )
        # cm = ( cm.astype('float') / cm.sum(axis=1)[:, np.newaxis] )
        # accuracies_per_bundle = cm.diagonal()

        testPredictedLabels = list( testPredictedLabels )
        testRealLabels = list( testRealLabels )
        for _tmp in testPredictedLabels :
            if _tmp not in testRealLabels :
                testRealLabels.append( _tmp )
        testRealLabels.sort()
        newPredictedDict = {}
        for _i in range( len( testRealLabels ) ) :
            _realLabel = testRealLabels[ _i ]
            if _realLabel == -1 :
                newPredictedDict[ _i ] = "unlabeledFibers"
            else :
                newPredictedDict[ _i ] = predictedDict[ _realLabel ]
        # if len( testRealLabels ) <= len( testPredictedLabels ) :
        #     for _i in range( len( testRealLabels ) ) :
        #         _realLabel = testRealLabels[ _i ]
        #         if _realLabel == -1 :
        #             newPredictedDict[ _i ] = "unlabeledFibers"
        #         else :
        #             newPredictedDict[ _i ] = predictedDict[ _realLabel ]
        # else :
        #     for _i in range( len( testPredictedLabels ) )  :
        #         _predLabel = testPredictedLabels[ _i ]
        #         if _predLabel == -1 :
        #             newPredictedDict[ _i ] = "unlabeledFibers"
        #         else :
        #             newPredictedDict[ _i ] = testRealLabels[ _predLabel ]

        saveDict( newPredictedDict, confusionMatrixDictSavingPath )
    else :
        confusion_matrix_model = readConfusionMatrix(
                                                     confusionMatrixSavingPath )
        newPredictedDict = readDict( confusionMatrixDictSavingPath )




    sensitivities = []
    specificities = []
    accuracies = []
    precisions = []
    jaccards = []
    nbLabelsNewDict = len( newPredictedDict.keys() )
    print( "Computing scores... " )
    for _label in newPredictedDict :
        print( f"Processing : [{_label+1}/{nbLabelsNewDict}]", end = "\r" )
        _sensitivity = computeSensitivity( confusion_matrix_model, _label )
        _specificity = computeSpecificity( confusion_matrix_model, _label )
        _accuracy = computeAccuracy( confusion_matrix_model, _label )
        _precision = computePrecision( confusion_matrix_model, _label )
        _jaccard = computeJaccard( confusion_matrix_model, _label )

        sensitivities.append( _sensitivity )
        specificities.append( _specificity )
        accuracies.append( _accuracy )
        precisions.append( _precision )
        jaccards.append( _jaccard )

    print( "\nDone" )

    saveScoresPerBundle2( sensitivities,
                          specificities,
                          accuracies,
                          precisions,
                          jaccards,
                          newPredictedDict,
                          out_scores_per_bundle_path )



    #--------------------------------------------------------------------------#

    # saveScoresPerBundle( accuracies_per_bundle,
                         # f1_score_per_bundle_model,
                         # precision_score_per_bundle_model,
                         # recall_score_per_bundle_model,
                         # newPredictedDict,
                         # out_scores_per_bundle_path )



def computeSensitivity( confusion_matrix_model, index_class ) :
    _tp = confusion_matrix_model[ index_class, index_class ]

    # Don't know why computing element_not_diag with confusion_matrix_model -
    # confusion_matrix_model.diagonal() is not working
    element_not_diag = np.copy( confusion_matrix_model )
    for i in range( element_not_diag.shape[ 0 ] ) :
        element_not_diag[ i, i ] = 0


    _fn = np.sum( element_not_diag[ index_class, : ] )

    if _tp == 0 :
        sensitivity = 0
    else :
        sensitivity = _tp / ( _tp + _fn )

    return( sensitivity )

def computeSpecificity( confusion_matrix_model, index_class ) :
    _tn = np.sum( confusion_matrix_model.diagonal() ) - confusion_matrix_model[
                                                      index_class, index_class ]
    # Don't know why computing element_not_diag with confusion_matrix_model -
    # confusion_matrix_model.diagonal() is not working
    element_not_diag = np.copy( confusion_matrix_model )
    for i in range( element_not_diag.shape[ 0 ] ) :
        element_not_diag[ i, i ] = 0

    _fn = np.sum( element_not_diag[ index_class, : ] )

    if _tn == 0 :
        specificity = 0
    else :
        specificity = _tn / ( _tn + _fn )

    return( specificity )

def computeAccuracy( confusion_matrix_model, index_class ) :
    _tp = confusion_matrix_model[ index_class, index_class ]
    _tn = np.sum( confusion_matrix_model.diagonal() ) - confusion_matrix_model[
                                                      index_class, index_class ]
    # Don't know why computing element_not_diag with confusion_matrix_model -
    # confusion_matrix_model.diagonal() is not working
    element_not_diag = np.copy( confusion_matrix_model )
    for i in range( element_not_diag.shape[ 0 ] ) :
        element_not_diag[ i, i ] = 0

    _fp = np.sum( element_not_diag[ :, index_class ] )
    _fn = np.sum( element_not_diag[ index_class, : ] )

    if ( _tp + _tn ) == 0 :
        accuracy = 0
    else :
        accuracy = ( _tp + _tn ) / ( _tp + _tn + _fp + _fn )

    return( accuracy )

def computePrecision( confusion_matrix_model, index_class ) :
    _tp = confusion_matrix_model[ index_class, index_class ]

    # Don't know why computing element_not_diag with confusion_matrix_model -
    # confusion_matrix_model.diagonal() is not working
    element_not_diag = np.copy( confusion_matrix_model )
    for i in range( element_not_diag.shape[ 0 ] ) :
        element_not_diag[ i, i ] = 0

    _fp = np.sum( element_not_diag[ :, index_class ] )

    if _tp == 0 :
        precision = 0
    else :
        precision = ( _tp ) / ( _tp + _fp )

    return( precision )

def computeJaccard( confusion_matrix_model, index_class ) :
    _tp = confusion_matrix_model[ index_class, index_class ]

    # Don't know why computing element_not_diag with confusion_matrix_model -
    # confusion_matrix_model.diagonal() is not working
    element_not_diag = np.copy( confusion_matrix_model )
    for i in range( element_not_diag.shape[ 0 ] ) :
        element_not_diag[ i, i ] = 0

    _fp = np.sum( element_not_diag[ :, index_class ] )
    _fn = np.sum( element_not_diag[ index_class, : ] )

    if _tp == 0 :
        jaccard = 0
    else :
        jaccard = ( _tp ) / ( _tp + _fn + _fp )

    return( jaccard )


"""
sensitivity (TP/(TP + FN))
specificity (TN/(TN + FN))
accuracy (TP + TN)/(TP + TN + FP + FN)
precision TP/(TP + FP)
Jaccard index (TP/(TP + FN + FP))
"""



################################################################################
################################################################################
################################################################################
def main() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    force = inputs[ "overwrite" ]

    predictedLabelsPath = inputs[ "predicted_labels" ]
    predictedDictPath = inputs[ "predicted_dictionary" ]
    trueLabelsPath = inputs[ "true_labels" ]
    trueDictPath = inputs[ "true_dictionary" ]
    isSupWMA = inputs[ "supWMA" ]

    outDir = inputs[ "output" ]
    if not os.path.isdir( outDir ) :
        os.mkdir( outDir )

    outConfusionMatrixPath = os.path.join( "confusionMatrix.tsv" )
    if ( os.path.isfile( outConfusionMatrixPath ) and not force ) :
        print( f"ERROR : confusion matrix file {outConfusionMatrixPath} already"
               f" exists and the force flag was not used" )
        sys.exit( 1 )

    outScoresPath = os.path.join( "scores.tsv" )
    if ( os.path.isfile( outScoresPath ) and not force ) :
        print( f"ERROR : confusion matrix file {outScoresPath} already exists "
               f" and the force flag was not used" )
        sys.exit( 1 )


    comparePredictionToTrue( trueLabelsPath, trueDictPath, predictedLabelsPath,
                                            predictedDictPath, isSupWMA, force )


if __name__ == "__main__":
    main()

################################################################################
################################################################################
################################################################################
# disco_darted_fused_dir = "/home/nv264568/Bureau/Doctorat/ISBI2023/fused_2_groups"
# dir_to_analyse = disco_darted_fused_dir
#
# """
# Good combinations
# tmpNum = "9"
# _num = 3
#
# """
#
# # tmpNum = "5"
# # tmpNum = "7"
# # tmpNum = "17"
#
# # tmpNum = "15"
# # tmpNum = "13"
# # tmpNum = "14"
# # tmpNum = "16"
# # realLabelsPath = f"{dir_to_analyse}/newLabels{tmpNum}/labels.txt"
# # realDictPath = f"{dir_to_analyse}/newLabels{tmpNum}/labels.dict"
#
# # realLabelsPath = f"{dir_to_analyse}/semiAutomaticLabels/labels.txt"
# # realDictPath = f"{dir_to_analyse}/semiAutomaticLabels/labels.dict"
#
# realLabelsPath = f"{dir_to_analyse}/semiAutomaticLabelsWithoutAtlasFinal/labels.txt"
# realDictPath = f"{dir_to_analyse}/semiAutomaticLabelsWithoutAtlasFinal/labels.dict"
# #-----------------------------------------------------------------------------#
#
# # _num = 10
# # _num = 6
# # _num = "12"
# _num = "2"
# # predictedLabelsPath = f"{dir_to_analyse}/testFullISBI{_num}/labels.txt"
# # predictedDictPath = f"{dir_to_analyse}/testFullISBI{_num}/labels.dict"
# predictedLabelsPath = f"/home/nv264568/Bureau/Doctorat/ISBI2023/fused_2_groups/tests/testClassicFull{_num}/labels.txt"
# predictedDictPath = f"/home/nv264568/Bureau/Doctorat/ISBI2023/fused_2_groups/tests/testClassicFull{_num}/labels.dict"
# #-----------------------------------------------------------------------------#
# # _numReco = "5"
# # _numReco = "8"
#
# # _numReco = "3"
# _numReco = "24"
#
# # predictedLabelsPathReco = f"{dir_to_analyse}/testFullReco{_numReco}/labels.txt"
# # predictedDictPathReco = f"{dir_to_analyse}/testFullReco{_numReco}/labels.dict"
#
# # predictedLabelsPathReco = f"{dir_to_analyse}/tests/testProjectAtlasRecobundles2/labels.txt"
# # predictedDictPathReco = f"{dir_to_analyse}/tests/testProjectAtlasRecobundles2/labels.dict"
#
# _numReco = 11
# predictedLabelsPathReco = f"/home/nv264568/Bureau/Doctorat/ISBI2023/fused_2_groups/tests/testProjectAtlasRecobundlesFull{_numReco}/labels.txt"
# predictedDictPathReco = f"/home/nv264568/Bureau/Doctorat/ISBI2023/fused_2_groups/tests/testProjectAtlasRecobundlesFull{_numReco}/labels.dict"
#
# # predictedLabelsPathReco = f"/home/nv264568/Bureau/Doctorat/ISBI2023/fused_2_groups/TGCC/testFullReco0.5/labels.txt"
# # predictedDictPathReco = f"/home/nv264568/Bureau/Doctorat/ISBI2023/fused_2_groups/TGCC/testFullReco0.5/labels.dict"
#
# # predictedLabelsPathReco = f"{dir_to_analyse}/testCombReco/labels.txt"
# # predictedDictPathReco = f"{dir_to_analyse}/testCombReco/labels.dict"
#
# # predictedLabelsPathReco = f"{dir_to_analyse}/testFullReco2/labels.txt"
# # predictedDictPathReco = f"{dir_to_analyse}/testFullReco2/labels.dict"
# #-----------------------------------------------------------------------------#
#
# # predictedLabelsSupWMA = f"{dir_to_analyse}/testSupWMA/predicted_labels.bin"
# # predictedDictSupWMA = f"/home/nv264568/Bureau/Doctorat/ISBI2023/labelsDictSupWMA.txt"
# # predictedLabelsSupWMA = f"{dir_to_analyse}/tests/testSupWMA/predicted_labels.bin"
# # predictedDictSupWMA = f"/home/nv264568/Bureau/Doctorat/ISBI2023/labelsDictSupWMA.txt"
# predictedLabelsSupWMA = f"/home/nv264568/Bureau/Doctorat/ISBI2023/fused_2_groups/tests/testSupWMA2/predicted_labels.bin"
# predictedDictSupWMA = f"/home/nv264568/Bureau/Doctorat/ISBI2023/fused_2_groups/tests/testSupWMA2/sf_clusters_train_labels_dictionary.txt"
#
#
# #-----------------------------------------------------------------------------#
# #-----------------------------------------------------------------------------#
# #-----------------------------------------------------------------------------#
#
# print( f"True labels directory : {os.path.dirname( realLabelsPath )}" )
# print( f"Predicted classic directory : {os.path.dirname( predictedLabelsPath )}" )
# print( f"Predicted reco directory : {os.path.dirname( predictedLabelsPathReco )}" )
#
# # print( "##################### Classic #####################" )
# # forceClassic = False
# # comparePredictionToTrue( realLabelsPath, realDictPath, predictedLabelsPath,
# #                      predictedDictPath, isSupWMA = False, force = forceClassic )
#
# # print( "##################### Geolab #####################" )
# # forceGeoLab = False
# # comparePredictionToTrue( realLabelsPath, realDictPath, predictedLabelsPathReco,
# #                   predictedDictPathReco, isSupWMA = False, force = forceGeoLab )
#
# print( "\n##################### SupWMA #####################" )
# forceSupWMA = False
# comparePredictionToTrue( realLabelsPath, realDictPath, predictedLabelsSupWMA,
#
# predictedDictSupWMA, isSupWMA = True, force = forceSupWMA )
#
#
#
# """
# -thrSim 0.05 -tolP 0.0 -tolThr 0.3 -tolMaxAng 0.0 -tolMaxDirAng 0.0 -tolMinShapeAng 0.9
# -tolMaxShapeAng 0.9 -tolLenght 0.0 -tolDBMP 0.0 -cp false -cb -sp false -v 1"
#
# """
