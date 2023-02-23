# GeoLab

This repository release the source-code from "Geolab : Geometric-based  tractography parcellation of Superficial white matter" accepted by the ISBI 2023.
It also makes available the altas used in the article.

GeoLab is a tool for superficial white matter parcellation, which improves the RecoBundles framework (made for deep white matter parcellation) to work on superficial white matter. This method outperforms the SOTAs on semi-ground truth (ARCHI dataset) and has a comparable result on the UKBiobank dataset.

![Alt text](https://github.com/vindasna/GeoLab/blob/main/Pipeline.png "Pipeline")


## Licence

The contents of this repository are released under Apache-2.0 license.

## Dependencies

- C++ compiler with full C++17 support
- Ubuntu 18.04.6 LTS
- CMake >= 3.20.2
- Boost >= 1.65.1
- Eigen >= 3.3
- OenMP >= 4.5
- Python >= 3.6.9
- Dipy >= 1.5
- Numpy >= 1.19.5
- Setproctitle >= 1.2.3
- fitter >= 1.5.1
- nibabel >= 3.2.2
- SupWMA if you want to apply the SupWMA trained with the ESBA atlas.

## Quick install

1. Install dependencies.
2. In ./GeoLab/CMakeLists.txt line 19 ("set(PYTHON_BINARY "/usr/bin/python3")") replace "/usr/bin/python3" by your python binary if you are using a virtual environment..
3. Clone Git repository and compile:

    $ git clone https://github.com/vindasna/GeoLab
   
    $ cd GeoLab
   
    $ mkdir build
   
    $ cd build
   
    $ cmake ..
   
    $ make
   

4. Configure PATH :
   Edit the startup ~/.bashrc or /etc/bash.bashrc file manually by adding this line :
   
    $ export PATH=/<edit as appropriate>/GeoLab/build/bin:$PATH
   
    $ export ESBA_DIR=/<edit as appropriate>/GeoLab/Atlas

5. Check installation :
 
    $ ProjectAtlasGeoLab -h

## Usage example on ESBA atlas

To extract the bundles of the ESBA atlas from a subject you first need to compute the tractogram (.tck/.trk/.bundles), register it to MNI space (recommended : image-based with ANTs) and resample it to 15 points per fiber. Then use the ProjectAtlasGeoLab command :

    $ ProjectAtlasGeoLab -i input${format} -o outputDir -nbPoints 15 -nbThreads ${nbThreads}

* Replace ${format} with {.trk, .tck, .bundles} according to your tractogram format.
* input${format} : subject's tractogram (.tck/.trk/.bundles)
* outputDir : directory where to save the results.
* ${nbThreads} : number of threads to use for OpenMP.
 

## Usage example on other atlas

First, you'll need to resample your atlas to a fixed number of points per fiber, if your atlas is in .tck format you can use MRTrix's command tckreample.

You'll need to analyse your atlas to get the bundle-specific thresholds **this step can be done only once** :

    $ analyseAtlasBundle.py -i atlasDir -f ${format} -r referenceImage.nii

* Replace ${format} with {.trk, .tck, .bundles} according to your tractogram format.
* atlasDir : Path to your atlas directory.
* referenceImage.nii : path to the reference .nii where the atlas is.


You'll also need to precompute the full atlas (all bundles in one single file), the atlas neighborhood and the atlas centroids :

    // Compute full atlas
    $ fuseAtlas -i atlasDir -o outDirFullAtlas -f ${format}

    // Compute atlas neighborhood
    $ computeNeighborhood -i outDirFullAtlas/fullAtlas${format} -o outDirNeighborhoodAtlas -r referenceImage.nii

    // Compute atlas centroids
    $ computeCentroids -i outDirNeighborhoodAtlas -o outDirCentroidsAtlas -r referenceImage.nii -nbPoints ${nbPoints} -nbThreads ${nbThreads} -f ${format}
    
* Replace ${format} with {.trk, .tck, .bundles} according to your tractogram format.
* atlasDir : Path to your atlas directory.
* outDirFullAtlas : output directory of command fuseAtlas.
* referenceImage.nii : path to the reference .nii where the atlas is
* outDirNeighborhoodAtlas : output directory of command computeNeighborhood.
* outDirCentroidsAtlas : output directory of command computeCentroids.
* ${nbThreads} : number of threads to use for OpenMP.


Then use the ProjectAtlasGeoLab command :

    $ ProjectAtlasGeoLab -i input${format} -a atlasDir -ref referenceImage.nii -o outputDir -nbPoints 15 -an NeigborhoodAtlas -anc CentroidsAtlas -nbThreads ${nbThreads}


* Replace ${format} with {.trk, .tck, .bundles} according to your tractogram format.
* input${format} : subject's tractogram (.tck/.trk/.bundles)
* atlasDir : Path to your atlas directory (after analysing it with analyseAtlasBundle).
* referenceImage.nii : path to the reference .nii where the atlas is
* outputDir : directory where to save the results.
* NeigborhoodAtlas : outDirNeighborhoodAtlas.
* CentroidsAtlas : outDirCentroidsAtlas.
* ${nbThreads} : number of threads to use for OpenMP.


## To compute the scores of prediction on labelled data

Your labelled data should be in the form of two files :

* .txt -> labels for each fiber in the form of :
* 
    `fiber_index_k : label_i`
    
    `           ...           `
    
    `fiber_index_l : label_j`
        
  With : 
    * label_i, ..., label_j integers.
            
    * fiber_index_l is the index of the fiber in the tractogram used as input for segmentation.
        
  If a fiber has multiple labels you just need to have several lines for that fiber.
        
        
* .dict -> dictionary for the labels in the form of :
    
    `label_name_i : label_i`
    
    `           ...          ` 
    
    `label_name_j : label_j`

  With : 
        
    * label_i, ..., label_j the same integers as in the .txt.
            
    * label_name_i, ..., label_name_j the names of the labels.
    
Once you have those files you can use the following command :
    
    $ scoresPredictedSGT.py -pl labels.txt -pd labels.dict -tl trueLabels.txt -td trueLabels.dict -o outDir

With : 
  * labels.txt : file produced by ProjectAtlasGeoLab, saved in output directory of ProjectAtlasGeoLab.
  * labels.dict : file produced by ProjectAtlasGeoLab, saved in output directory of ProjectAtlasGeoLab.
  * trueLabels.txt : your labelled data as explained above.
  * trueLabels.dict : your labelled data as explained above.
  * outDir : directory where to save scoresPerBundle.tsv file. This file contains the scores per bundle (label).

If you want to reproduce the result of the paper, the semi-ground truth is in the SGT folder.


## Apply SupWMA model trained with ESBA atlas

First you need to extract the features for SupWMA with extract_bundles_feat.py command :

    $ extract_bundles_feat.py -i SGT.bundles -o outSGT.h5 -v 1

With :
  * SGT.bundles : path to your input tractogram in .bundles format (compatibility with other formats will be added soon)
  * outSGT.h5 : path of the output file which must be .h5

Then you can use the applySupWMA command :

    $ applySupWMA.py -t tractogram.bundles -f tractogram.h5 -ep encoderParameters.pickle -ew encoderWeights.pth -cw classifierWeights.pth -ln labelNames.h5 -ld labelsDictSupWMA.txt -spw SupWMA_path -o outDir

With :
  * tractogram.bundles : your input tractogram in .bundles (compatibility with other formats will be added soon)
  * tractogram.h5 : the output of extract_bundles_feat.py command.
  * encoderParameters.pickle : found in GeoLab/TrainedSupWMA
  * encoderWeights.pth : found in GeoLab/TrainedSupWMA
  * classifierWeights.pth : found in GeoLab/TrainedSupWMA
  * labelNames.h5 : found in GeoLab/TrainedSupWMA
  * labelsDictSupWMA.txt :found in GeoLab/TrainedSupWMA
  * SupWMA_path : path to where you cloned the SupWMA repository.
  * outDir : output directory where to save the results.


## For windows

GeoLab will soon be available for windows as a docker container. 
