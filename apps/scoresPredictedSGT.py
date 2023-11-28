#!${PYTHON_BINARY}

import os, sys, shutil

import json

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
        "-pd", "--predicted-dictionary",
        type=str, metavar="<path>", default="",
        help=("Path to the predicted labels dictionary .dict (MANDATORY "
              "for GeoLab and SupWMA)" ) )
    parser.add_argument(
        "-force", "--overwrite",
        action='store_true', default=False,
        help="Overwrite output file even if it already exists" )
    parser.add_argument(
        "-supWMA", "--supWMA",
        action='store_true', default=False,
        help="Use flag if labels come from SupWMA" )
    parser.add_argument(
        "-rbx", "--rbx",
        action='store_true', default=False,
        help="Use flag if labels come from RecoBundlesX" )
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

def saveLabels( in_list, path ) :
    nbStreamlines = len( in_list )
    with open( path, 'w' ) as f :
        for tmpFiberIndex in range( nbStreamlines ) :
            tmpCounter = 0
            for tmpLabel in in_list[ tmpFiberIndex ] :
                f.write( f"{tmpFiberIndex} : {tmpLabel}" )
                if tmpFiberIndex != ( nbStreamlines - 1 ) or tmpCounter != ( len( in_list[ tmpFiberIndex ] ) - 1 ) :
                    f.write( "\n" )

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


def readRecoBundlesX( path, nbStreamlines ) :
    out_labels_dict = {}
    out_labels = {}
    for tmpFiber in range( nbStreamlines ) :
        out_labels[ tmpFiber ] = [ -1 ]

    with open( path, 'r' ) as f :
        tmpDict = json.load( f )
    
    counter = 0
    for tmpFile in tmpDict.keys() :
        tmpBundleName = os.path.basename( tmpFile )
        tmpBundleName = tmpBundleName.replace( ".trk", "" )
        tmpLabel = counter
        if tmpLabel not in out_labels_dict.keys() :
            out_labels_dict[ counter ] = tmpBundleName
            counter += 1
        
        for tmpIndex in tmpDict[ tmpFile ][ "indices" ] :
            if len( out_labels[ tmpIndex ] ) == 1 and out_labels[ tmpIndex ][ 0 ] == -1 :
                out_labels[ tmpIndex ] = [ tmpLabel ]
            else :
                out_labels[ tmpIndex ].append( tmpLabel )
    
    """
    out_dict_path = "/volatile/articles_Nabil/GeoLab/RecoBundlesX/labels.dict"
    saveDict( out_labels_dict, out_dict_path )
    out_labels_path = "/volatile/articles_Nabil/GeoLab/RecoBundlesX/labels.txt"
    saveLabels( out_labels, out_labels_path )
    """
    
    return( out_labels_dict, out_labels )


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
    # print( f"Error in getLabelFromBundleName() : \'{bundleName}\' not found in "
    #        f"bundlesDict" )
    # print( f"Bundles names in bundlesDict : " )
    # for label in bundlesDict.keys() :
    #     print( f"\'{bundlesDict[ label ]}\'" )
    sys.exit( 1 )


def saveScoresPerBundle( sensitivities,
                         specificities,
                         accuracies,
                         precisions,
                         jaccards,
                         f1_scores,
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
    f1Len = len( f1_scores )
    lnLen = len( sortedLabelsNames )
    if ( senLen != speLen ) :
        print( f"ERROR in saveScoresPerBundle() : Specificities must be the "
               f"same size as sensitivities, got {speLen} and {senLen} "
               f"respectivly" )
        sys.exit( 1 )
    if ( senLen != accLen ) :
        print( f"ERROR in saveScoresPerBundle() : accuracies must be the same "
               f"size as sensitivities, got {accLen} and {senLen} respectivly" )
        sys.exit( 1 )
    if ( senLen != preLen ) :
        print( f"ERROR in saveScoresPerBundle() : precisions must be the same "
               f"size as sensitivities, got {preLen} and {senLen} respectivly" )
        sys.exit( 1 )
    if ( senLen != jacLen ) :
        print( f"ERROR in saveScoresPerBundle() : jaccards must be the same "
               f"size as sensitivities, got {jacLen} and {senLen} respectivly" )
        sys.exit( 1 )
    if ( senLen != f1Len ) :
        print( "ERROR in saveScoresPerBundle() : f1-scores must be the same"
               f" size as sensitivities, got {f1Len} and {senLen} respectivly" )
        sys.exit( 1 )
    if ( senLen != lnLen ) :
        print( "ERROR in saveScoresPerBundle() : labels names must be the same"
               f" size as sensitivities, got {lnLen} and {senLen} respectivly" )
        sys.exit( 1 )

    nbBundles = len( sensitivities )

    with open( out_scores_per_bundle_path, 'w' ) as f :
        f.write( "BundleName\tSensitivity\tSpecificity\tAccuracy\tPrecision\t"
                 "Jaccard\tF1-score\n" )
        # for i in range( nbBundles ) :
        for i in sortedLabelsNames :
            _sensitivity = sensitivities[ i ]
            _specificity = specificities[ i ]
            _accuracy = accuracies[ i ]
            _precision = precisions[ i ]
            _jaccard = jaccards[ i ]
            _f1_score = f1_scores[ i ]
            _lableName = sortedLabelsNames[ i ]
            f.write( f"{_lableName}\t{_sensitivity}\t{_specificity}\t"
                     f"{_accuracy}\t{_precision}\t{_jaccard}\t{_f1_score}\n" )

    #--------------------------------------------------------------------------#
    mean_sen = np.mean( sensitivities )
    median_sen = np.median( sensitivities )
    std_sen = np.std( sensitivities )

    mean_spe = np.mean( specificities )
    median_spe = np.median( specificities )
    std_spe = np.std( specificities )

    mean_acc = np.mean( accuracies )
    median_acc = np.median( accuracies )
    std_acc = np.std( accuracies )

    mean_pre = np.mean( precisions )
    median_pre = np.median( precisions )
    std_pre = np.std( precisions )

    mean_jac = np.mean( jaccards )
    median_jac = np.median( jaccards )
    std_jac = np.std( jaccards )

    mean_f1 = np.mean( f1_scores )
    median_f1 = np.median( f1_scores )
    std_f1 = np.std( f1_scores )

    print( f"Scores per bundles ( mean (median) +- std ) :\n"
           f" Sensitivity : {mean_sen} ({median_sen}) +- {std_sen} \n"
           f" Specificity : {mean_spe} ({median_spe}) +- {std_spe} \n"
           f" Accuracy : {mean_acc} ({median_acc}) +- {std_acc} \n"
           f" Precision : {mean_pre} ({median_pre}) +- {std_pre} \n"
           f" Jaccard : {mean_jac} ({median_jac}) +- {std_jac}\n"
           f" F1-score : {mean_f1} ({median_f1}) +- {std_f1}" )

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
                  predictedDictPath, outDir, isSupWMA = False, isRBX = False, 
                                                               force = False ) :
    if isSupWMA and isRBX :
        print( f"ERROR : isSupWMA and isRBX both True at the same time" )
        sys.exit( 1 )
    confusionMatrixSavingPath = os.path.join( outDir, f"confusionMatrix.tsv" )

    confusionMatrixDictSavingPath = os.path.join( outDir,
                                                       f"confusionMatrix.dict" )

    out_scores_per_bundle_path = os.path.join( outDir, f"scoresPerBundle.tsv" )

    if ( not os.path.isfile( confusionMatrixSavingPath ) or
                not os.path.isfile( confusionMatrixDictSavingPath ) or force ) :
        print( "Reading true labels... ", end = "" )
        realLabels = readLabels( realLabelsPath )
        print( "Done\nReading true labels dictionary... ", end = "" )
        realDict = readDict( realDictPath )

        print( "Done\nReading predicted labels... ", end = "" )
        if ( isSupWMA ) :
            predictedLabels = readLabelsSupWMA( predictedLabelsPath )
        elif isRBX :
            nbStreamlines = len( realLabels )
            predictedDict, predictedLabels = readRecoBundlesX( predictedLabelsPath, nbStreamlines )
        else :
            predictedLabels = readLabels( predictedLabelsPath )
        print( "Done\nReading predicted labels dictionary... ", end = "" )
        if not isRBX :
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


        ##################################################################################################################
        # Putting missing labels in dictionaries
        """
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
        """
        
        # Getting only labels in SGT
        tmpNewPredictedDict = {}
        for tmpKey in predictedDict.keys() :
            tmpPredictedBundleName = predictedDict[ tmpKey ]
            tmpRealLabel = getLabelFromBundleName( realDict, tmpPredictedBundleName )
            if tmpRealLabel :
                tmpNewPredictedDict[ tmpKey ] = predictedDict[ tmpKey ]
        
        # Rebuild labels and dict to start bundle label at 0 and end at the number of labels bundles
        print( "Rebuilding predicted dictionnary... " )
        newPredictedLabelsDict = {}
        correspondenceDict = {}
        counter = 0
        tmpNbKeys = len( tmpNewPredictedDict.keys() )
        for tmpLabel in tmpNewPredictedDict.keys() :
            print( f"Rebuilding dictionnary... Processing : [{counter+1}/{tmpNbKeys}]", end = "\r" )
            newPredictedLabelsDict[ counter ] = tmpNewPredictedDict[ tmpLabel ]
            correspondenceDict[ tmpLabel ] = counter
            counter += 1

        print( "\nDone" )

        nbFibers = len( predictedLabels.keys() )
        counter = 1
        for tmpFiberIndex in predictedLabels.keys() :
            print( f"Rebuilding predicted labels... Processing : [{counter}/{nbFibers}]", end = "\r" )
            tmpI = 0
            while tmpI < len( predictedLabels[ tmpFiberIndex ] ) :
                tmpLabel = predictedLabels[ tmpFiberIndex ][ tmpI ]
                if tmpLabel != -1 :
                    if tmpLabel in tmpNewPredictedDict.keys() :
                        predictedLabels[ tmpFiberIndex ][ tmpI ] = correspondenceDict[ predictedLabels[ tmpFiberIndex ][ tmpI ] ]
                        tmpI += 1
                    else :
                        # Pas d'incrementation de tmpI car on enlève un élement de la liste
                        while tmpLabel in predictedLabels[ tmpFiberIndex ] :
                            predictedLabels[ tmpFiberIndex ].remove( tmpLabel )

                        if len( predictedLabels[ tmpFiberIndex ] ) == 0 :
                            predictedLabels[ tmpFiberIndex ] = [ -1 ]

                else :
                    tmpI += 1

                            
            counter += 1
        

        predictedDict = newPredictedLabelsDict
        

        ##################################################################################################################

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
                        try :
                            _tmpLabel =  getLabelFromBundleName( predictedDict,
                                              realDict[ realLabels[ i ][ 0 ] ] )
                        except :
                            _tmpLabelName = realDict[ realLabels[ i ][ 0 ] ]
                            if _tmpLabelName not in predictedDict.values() :
                                _tmpLabel = max( [ int ( i ) for i in
                                                    predictedDict.keys() ] ) + 1
                                predictedDict[ _tmpLabel ] = _tmpLabelName
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
                        yTrue.append( -1 )
                    else :
                        try :
                            _tmpLabel =  getLabelFromBundleName( predictedDict,
                                              realDict[ realLabels[ i ][ 0 ] ] )
                        except :
                            _tmpLabelName = realDict[ realLabels[ i ][ 0 ] ]
                            if _tmpLabelName not in predictedDict.values() :
                                _tmpLabel = max( [ int ( i ) for i in
                                                    predictedDict.keys() ] ) + 1
                                predictedDict[ _tmpLabel ] = _tmpLabelName
                        yTrue.append( _tmpLabel )


        print( "\nGetting predicted labels" )
        testPredictedLabels = np.unique( yPred )
        print( f"Number of labels in yPred : {len(testPredictedLabels)}" )

        print( "Getting real labels" )
        testRealLabels = np.unique( yTrue )
        print( f"Number of labels in yTrue : {len(testRealLabels)}" )

        tmpLabels = list ( np.unique( list( testPredictedLabels ) + list( testRealLabels ) ) )
        tmpLabelsDict = {}
        newPredictedDict = {}
        for tmpI in range( len( tmpLabels ) ) :
            tmpLabelsDict[ tmpLabels[ tmpI ] ] = tmpI
            if tmpI not in newPredictedDict.keys() :
                if tmpLabels[ tmpI ] == -1 :
                    newPredictedDict[ tmpI ] = "Unlabelled"
                else :
                    newPredictedDict[ tmpI ] = predictedDict[ tmpLabels[ tmpI ] ]

        newYpred = []
        newYtrue = []
        for tmpI in range( len( yPred ) ) :
            tmpYpred = yPred[ tmpI ]
            tmpYtrue = yTrue[ tmpI ]
            newYpred.append( tmpLabelsDict[ tmpYpred ] )
            newYtrue.append( tmpLabelsDict[ tmpYtrue ] )

        confusion_matrix_model = confusion_matrix( newYtrue, newYpred )
        saveConfusionMatrix( confusion_matrix_model, confusionMatrixSavingPath )

        

        # testPredictedLabels = list( testPredictedLabels )
        # testRealLabels = list( testRealLabels )
        # for _tmp in testPredictedLabels :
        #     if _tmp not in testRealLabels :
        #         testRealLabels.append( _tmp )
        # testRealLabels.sort()
        # newPredictedDict = predictedDict

        # newPredictedDict = {}
        # for _i in range( len( testRealLabels ) ) :
        #     _realLabel = testRealLabels[ _i ]
        #     if _realLabel == -1 :
        #         newPredictedDict[ _i ] = "unlabeledFibers"
        #     else :
        #         newPredictedDict[ _i ] = predictedDict[ _realLabel ]


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
    f1_scores = []
    nbLabelsNewDict = len( newPredictedDict.keys() )
    print( "Computing scores... " )
    for _label in newPredictedDict :
        print( f"Processing : [{_label+1}/{nbLabelsNewDict}]", end = "\r" )
        _sensitivity = computeSensitivity( confusion_matrix_model, _label )
        _specificity = computeSpecificity( confusion_matrix_model, _label )
        _accuracy = computeAccuracy( confusion_matrix_model, _label )
        _precision = computePrecision( confusion_matrix_model, _label )
        _jaccard = computeJaccard( confusion_matrix_model, _label )
        if ( _sensitivity + _precision == 0 ) :
            _f1_score = 0
        else :
            _f1_score = ( 2 * ( _sensitivity * _precision ) / ( _sensitivity +
                                                                  _precision ) )

        sensitivities.append( _sensitivity )
        specificities.append( _specificity )
        accuracies.append( _accuracy )
        precisions.append( _precision )
        jaccards.append( _jaccard )
        f1_scores.append( _f1_score )

    print( "\nDone" )

    saveScoresPerBundle( sensitivities,
                         specificities,
                         accuracies,
                         precisions,
                         jaccards,
                         f1_scores,
                         newPredictedDict,
                         out_scores_per_bundle_path )




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
    isRBX = inputs[ "rbx" ]

    if ( isSupWMA or ( not isSupWMA and not isRBX ) ) and not os.path.isfile( predictedDictPath ) :
        print( f"ERROR : -pd must be given for SupWMA or GeoLab" )
        sys.exit( 1 )

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
                                    predictedDictPath, outDir, isSupWMA, isRBX, force )


if __name__ == "__main__":
    main()
