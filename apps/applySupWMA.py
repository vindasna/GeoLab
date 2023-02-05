#!/usr/bin/python3
import sys, os, subprocess, shutil

import time
import numpy as np
import h5py

import argparse
import textwrap
from argparse import RawTextHelpFormatter

import pickle

import torch
import torch.nn.parallel
import torch.utils.data

from utils.model_supcon import PointNet_SupCon, PointNet_Classifier
from utils.model import PointNetCls
from utils.dataset import TestDataset

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
python3 countSuccessfulSubjects.py \
-in inputDir
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
        prog="python3 applyPointNetClassifier.py",
        description=textwrap.dedent(DOC),
        formatter_class=RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-t", "--tractogram",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the input tractogram (.bundles )" ) )

    required.add_argument(
        "-f", "--features",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the features ( points in .h5 )" ) )

    required.add_argument(
        "-ep", "--encoder-params",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the encoder parameters ( .pickle )" ) )

    required.add_argument(
        "-ew", "--encoder-weights",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the encoder weights (.pth )" ) )

    required.add_argument(
        "-cw", "--classifier-weights",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the classifier weights (.pth )" ) )

    required.add_argument(
        "-ln", "--labels-names",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the label names of the trained model ( .h5 ) " ) )

    required.add_argument(
        "-ld", "--labels-dictionary",
        type=is_file, required=True, metavar="<path>",
        help=( "Path to the labels dictionary (.txt) " ) )

    required.add_argument(
        "-o", "--output",
        type=str, required=True, metavar="<path>",
        help=( "Output directory" ) )

    # Optional arguments
    parser.add_argument(
        "-bs", "--batch-size",
        type=int, default=6144,
        help="Batch size ( default = 6144 ) " )

    parser.add_argument(
        "-nw", "--num-workers",
        type=int, default=4,
        help="Number of data loading workers ( default = 4 )" )

    parser.add_argument(
        "-apl", "--applyPredictedLabels",
        type=is_file,
        help="Path to the applyPredictedLabels binary" )

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


def readBundlesFile( bundle_filename, verbose ):

    if verbose > 1 :
        print (f'Reading ... ... {bundle_filename}')

    # Checking extension bundle_filename
    if bundle_filename.replace(".bundles", "") == bundle_filename:
        bundle_filename = f'{bundle_filename}.bundles'
    else:
        if bundle_filename.replace(".bundlesdata","") != bundle_filename:
            bundle_filename = bundle_filename.replace(".bundlesdata",".bundles")

    ns = dict()
    exec(open(bundle_filename).read(), ns)
    curves_count = (ns[ 'attributes' ][ 'curves_count' ])
    labels = (ns[ 'attributes' ][ 'bundles' ])
    try:
        resolution = [ ns[ 'attributes' ][ 'resolutionX' ],
                        ns[ 'attributes' ][ 'resolutionY' ],
                        ns[ 'attributes' ][ 'resolutionZ' ] ]
    except:
        print("Not resolution info in .bundles file...")
        resolution = [-1, -1, -1]
    try:
        size = [ ns[ 'attributes' ][ 'sizeX' ],
                    ns[ 'attributes' ][ 'sizeY' ],
                    ns[ 'attributes' ][ 'sizeZ' ] ]
    except:
        print("Not size info in .bundles file...")
        size = [-1, -1, -1]


    bundle = []
    nPoints = []
    if curves_count > 0:
       f = open(bundle_filename.replace(".bundles", ".bundlesdata"),'rb')
       for curve in range(curves_count):
           p = np.fromfile(f, dtype=np.int32, count=1)
           nPoints.append(p[0])
           bundle.append(np.fromfile(f, dtype=np.float32, count=p[0]*3))
       f.close()

       if int(max(nPoints)) == int(min(nPoints)):
          n_points = int(max(nPoints))
       else:
             print ("Number of points per curve not equal...")
             n_points = nPoints
    else:
          ns = dict()
          exec(open(bundle_filename + '.bundles').read(), ns)
          n_points = (ns[ 'attributes' ][ 'number_of_points' ])

    return bundle, curves_count, n_points, labels, resolution, size


def writeBundlesFormat (save_path, curves_count, labels = ["1"],
                        resolution = [-1, -1, -1], size = [-1, -1, -1]):
    if save_path.replace(".bundles", "") == save_path:
        save_path = f'{save_path}.bundles'
    else:
        if save_path.replace(".bundlesdata","") != save_path:
            save_path = save_path.replace(".bundlesdata", ".bundles")
    save = open(save_path,'w')
    minf = ("attributes = {"
            "\n     'binary' : 1,"
            "\n     'bundles' : %s,"
            "\n     'byte_order' : 'DCBA',"
            "\n     'curves_count' : %s,"
            "\n     'data_file_name' : '*.bundlesdata',"
            "\n     'format' : 'bundles_1.0',"
            "\n     'io_mode' : 'binary',"
            "\n     'object_type' : 'BundleMap',"
            "\n     'space_dimension' : 3,"
            "\n     'resolutionX' : %s,"
            "\n     'resolutionY' : %s,"
            "\n     'resolutionZ' : %s,"
            "\n     'sizeX' : %s,"
            "\n     'sizeY' : %s,"
            "\n     'sizeZ' : %s"
            "\n     }")
    save.write( minf % ( [ str( labels ), 0 ], curves_count, *resolution,
                                                                       *size ) )
    save.close()


def writeBundlesdataFile(file_name, bundle, n_points):
    if file_name.replace(".bundlesdata","") == file_name:
        file_name = f'{file_name}.bundlesdata'

    print (f'Writing ... ... {file_name}')

    curves_count = len(bundle)
    f = open(file_name,'wb')
    if type(n_points) == int:
        print("Same number of points for every fiber.")
        for i in range(curves_count):
            np.asarray(n_points, dtype=np.int32).tofile(f)
            np.asarray(bundle[i], dtype=np.float32).tofile(f)
    elif type(n_points) == list:
        print("Different number of points for fibers.")
        if curves_count == len(n_points):
            for i in range(curves_count):
                np.asarray(n_points[i], dtype=np.int32).tofile(f)
                np.asarray(bundle[i], dtype=np.float32).tofile(f)
        else:
            print("Number of bundles and number of elements in n_points don't"
                    "match.")
            f.close()
            return(1)
    else:
        print("Invalid type of n_points.")
        f.close()
        return(1)

    f.close()
    return(0)

def convertBundleVectorToMatrix( bundle, curves_count, nbPoints ) :
    # bundle_matrix = np.zeros( ( curves_count, nbPoints, 3 ) )
    # for curve in range( curves_count ) :
    #     for point in range( nbPoints ) :
    #         for i in range( 3 ) :
    #             if np.isnan( bundle[ curve ][ 3 * point + i ] ) :
    #                 print( f"\nERROR : NaN value in bundle {bundle_name}" )
    #                 sys.exit( 1 )
    #             bundle_matrix[ curve, point, i ] = bundle[ curve ][
    #                                                              3 * point + i ]

    bundle_matrix = np.array( bundle ).reshape( ( curves_count, nbPoints, 3 ) )
    return( bundle_matrix )

def convertBundleMatrixToVector( bundle_matrix, curves_count, nbPoints ) :
    bundle_vector = bundle_matrix.reshape( curves_count * nbPoints * 3 )
    return( bundle_vector )




def load_test_data( feat_path, label_names_path, test_batch_size, num_workers,
                                                                 script_name ) :
    """Load test data and labels name in model"""
    # Put test data into loader
    test_dataset = TestDataset( feat_path, None, label_names_path, script_name)
    test_loader = torch.utils.data.DataLoader(
                                         test_dataset,
                                         batch_size = test_batch_size,
                                         shuffle = False,
                                         num_workers = int( num_workers ) )
    test_data_size = len( test_dataset )
    print( script_name, 'The test data size is:{}'.format( test_data_size ) )
    num_classes = len( test_dataset.label_names_in_model )
    # load label names
    label_names = test_dataset.label_names_in_model
    print( 'The label names are: {}'.format( str( label_names ) ) )
    print( script_name, 'The number of classes is:{}'.format( num_classes ) )

    return test_loader, label_names, num_classes


def load_model_point_net( num_classes, classifier_weight_path ):
    classifer = PointNetCls( k = num_classes )

    classifer.load_state_dict( torch.load( classifier_weight_path ) )

    return classifer


def load_model_supWMA( encoder_params, encoder_weight_path,
                                         classifier_weight_path, num_classes ) :
    encoder = PointNet_SupCon( head=encoder_params[ 'head_name' ],
                    feat_dim=encoder_params[ 'encoder_feat_num' ] ).to( device )
    # print( '{} use first feature transform for encoder'.format(
    #                      not encoder_params[ 'not_first_feature_transform' ] ) )
    classifer = PointNet_Classifier( num_classes = num_classes ).to( device )

    # load weights
    encoder.load_state_dict( torch.load( encoder_weight_path ) )
    
    classifer.load_state_dict( torch.load( classifier_weight_path ) )
    
    return encoder, classifer


def test_net( output_prediction_mask_path, num_classes, test_data_loader,
              encoder_params, encoder_weight_path, classifier_weight_path,
                                                                 script_name ) :
    """perform predition of multiple clusters"""
    print('')
    print('===================================')
    print('')
    print(script_name, 'Start multi-cluster prediction.')

    # classifer_net = load_model_point_net(  num_classes, classifier_weight_path )
    encoder_net, classifer_net = load_model_supWMA( encoder_params,
                                   encoder_weight_path, classifier_weight_path,
                                                                   num_classes )
    if not os.path.exists( output_prediction_mask_path ) :
        # Load model
        start_time = time.time()
        with torch.no_grad():
            total_test_correct = 0
            test_labels_lst, test_predicted_lst = [], []
            for j, data in ( enumerate( test_data_loader, 0 ) ) :
                points, labels = data
                points = points.transpose( 2 , 1 )
                points = points.to( device )

                # classifer_net = classifer_net.eval()
                # cls
                # features = points
                # pred = classifer_net( features )
                # _, pred_idx = torch.max( pred, dim = 1 )
                # pred_idx = torch.where( pred_idx < num_classes - 1, pred_idx,
                #                   torch.tensor( num_classes - 1 ).to( device ) )
                #
                # # entire data
                # pred_idx = pred_idx.cpu().detach().numpy().tolist()
                # test_predicted_lst.extend( pred_idx )
                #
                # # enc-cls
                # features, _, _, critical_point_idx = encoder_net.encoder(points)

                encoder_net, classifer_net = ( encoder_net.eval(),
                                                          classifer_net.eval() )
                features = encoder_net.encoder( points ) # Added by Nabil
                pred = classifer_net( features )
                _, pred_idx = torch.max(pred, dim=1)
                pred_idx = torch.where( pred_idx < num_classes, pred_idx,
                                      torch.tensor( num_classes ).to( device ) )

                pred_idx = pred_idx.cpu().detach().numpy().tolist()
                test_predicted_lst.extend( pred_idx )

        end_time = time.time()
        print( 'The total time of prediction is:{} s'.format( round( (
                                                end_time - start_time ), 4 ) ) )
        print( 'The test sample size is: ', len( test_predicted_lst ) )
        # test_prediction_lst_h5 = h5py.File( output_prediction_mask_path, "w" )
        # test_prediction_lst_h5[ 'complete_pred_test' ] = test_predicted_lst
        test_predicted_array = np.asarray( test_predicted_lst )

        # Saving binary with labels
        f = open( output_prediction_mask_path, 'wb' )
        nbFibers = test_predicted_array.shape[ 0 ]
        for i in range( nbFibers ) :
            np.asarray( test_predicted_array[ i ], dtype=np.int16).tofile( f )
        f.close()

        return test_predicted_array

    else:
        # print( script_name, 'Loading prediction result.' )
        # test_prediction_h5 = h5py.File( output_prediction_mask_path, "r" )
        # test_predicted_array = np.asarray(
        #                             test_prediction_h5[ 'complete_pred_test' ] )
        print( f"File {output_prediction_mask_path} already exists" )
        return 0




def tractography_parcellation( input_tractogram_path,
                               labels_dictionary_path,
                               output_prediction_mask_path,
                               output_cluster_folder,
                               applyPredictedLabels_command,
                               script_name ) :
    """Generate the tractography parcellation results with the predicted list"""
    print( '' )
    print( '===================================' )
    print( '' )
    print( script_name, 'Output fiber clusters.' )
    # Added by Nabil
    print( f"Input tractogram path : {input_tractogram_path}" )
    print( f"Labels dictionary path: {labels_dictionary_path}" )
    print( f"Output prediction mask path : {output_prediction_mask_path}" )
    print( f"Output cluster folder : {output_cluster_folder}" )
    print( f"applyPredictedLabels command : {applyPredictedLabels_command}" )
    print( f"Script name : {script_name}" )
    # Tractography Parcellation
    tractographyParcellation_command = [ applyPredictedLabels_command,
                                         "-i", input_tractogram_path,
                                         "-l", output_prediction_mask_path,
                                         "-o", output_cluster_folder,
                                         "-d", labels_dictionary_path,
                                         "-supWMA",
                                         "-v", str( 1 ) ]
    run_sh_process( cmd = tractographyParcellation_command, shell = True )


    # Reading dictionary with labels names
    file = open( labels_dictionary_path, "r" )
    lines = file.readlines()
    file.close
    labels_dictionary = {}
    for line in lines :
        labels_dictionary[ line.split()[ -1 ] ] = line.split()[ 0 ]


    # Renaming files
    clusters = os.listdir( output_cluster_folder )
    for cluster in clusters :
        if cluster.endswith( ".bundles" ) :
            key = cluster.replace( "cluster_", "" )
            key = key.replace( ".bundles", "" )
            old_cluster_name = os.path.join( output_cluster_folder, cluster )
            new_cluster_name = os.path.join( output_cluster_folder,
                                         f"{labels_dictionary[ key ]}.bundles" )
            os.rename( old_cluster_name, new_cluster_name )
            os.rename( f"{old_cluster_name}data", f"{new_cluster_name}data" )




    print( script_name, 'Done! Clusters are in :', output_cluster_folder )


def main() :
    """
    Parse the command line.
    """
    inputs, verbose = get_cmd_line_args()

    input_tractogram_path = inputs[ "tractogram" ]

    feat_path = inputs[ "features" ]

    label_names_path = inputs[ "labels_names" ]

    labels_dictionary_path = inputs[ "labels_dictionary" ]

    test_batch_size = inputs[ "batch_size" ]

    num_workers = inputs[ "num_workers" ]

    applyPredictedLabels_command = inputs[ "applyPredictedLabels" ]
    if applyPredictedLabels_command :
        if not os.path.isfile( applyPredictedLabels_command ) :
            applyPredictedLabels_command = "applyPredictedLabels"
    else :
        applyPredictedLabels_command = "applyPredictedLabels"

    encoder_params_path = inputs[ "encoder_params" ]

    encoder_weight_path = inputs[ "encoder_weights" ]

    classifier_weight_path = inputs[ "classifier_weights" ]


    output_dir = inputs[ "output" ]
    if not os.path.isdir( output_dir ) :
        os.mkdir( output_dir )

    output_prediction_mask_path = os.path.join( output_dir,
                                                        "predicted_labels.bin" )

    output_cluster_folder = os.path.join( output_dir, "Clusters" )
    if not os.path.isdir( output_cluster_folder ) :
        os.mkdir( output_cluster_folder )

    script_name = "applyPointNetClassifier"

    # Load test data
    test_data_loader, label_names, num_classes = load_test_data(
                                            feat_path = feat_path,
                                            label_names_path = label_names_path,
                                            test_batch_size = test_batch_size,
                                            num_workers = num_workers,
                                            script_name = script_name )

    # Load encode parameters
    with open( encoder_params_path, 'rb' ) as f:
        encoder_params = pickle.load( f )
        f.close()

    # Generate prediction
    predicted_arr =  test_net(
                      output_prediction_mask_path = output_prediction_mask_path,
                      num_classes = num_classes,
                      test_data_loader = test_data_loader,
                      encoder_params = encoder_params,
                      encoder_weight_path = encoder_weight_path,
                      classifier_weight_path = classifier_weight_path,
                      script_name = script_name )

    # # Process tractography parcellation
    tractography_parcellation( input_tractogram_path,
                               labels_dictionary_path,
                               output_prediction_mask_path,
                               output_cluster_folder,
                               applyPredictedLabels_command,
                               script_name )




if __name__ == "__main__" :
    # Set device
    use_cpu = True
    if use_cpu:
        device = torch.device("cpu")
    else:
        device = torch.device("cuda:0")
        torch.cuda.empty_cache()
    t1 = time.time()
    main()
    elapsed_time = time.time() - t1
    print( f"Elapsed time : {elapsed_time}" )
