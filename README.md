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

## Quick install

1. Install dependencies.
2. Clone Git repository and compile:

   `$ git clone https://github.com/vindasna/GeoLab`
   
   `$ cd GeoLab`
   
   `$ mkdir build`
   
   `$ cd build`
   
   `$ cmake ..`
   
   `$ make`
   

3. Configure PATH :
   Edit the startup ~/.bashrc or /etc/bash.bashrc file manually by adding this line :
   
   `$ export PATH=/<edit as appropriate>/GeoLab/build/bin:$PATH`

4. Check installation :
 
   `$ ProjectAtlasGeoLab -h`

5. If you are using a virtual environment for python, change the first line in {GeoLab Path}/build/bin/analyseAtlasBundle.py :
   
   `$ #!/usr/bin/python3 (original line)`    --->    `$ #!{your python binary}`     
 

## Usage example on ESBA atlas

To extract the bundles of the ESBA atlas from a subject you first need to compute the tractogram (.tck/.trk/.bundles), register it to MNI space (recommended : image-based with ANTs) and resample it to 15 points per fiber. If you use .tck/.trk you need to create a ".minf" for the file :

    `$ createMinf -o {path to tractogram without extension}.minf -f ${format}`

Then use the ProjectAtlasGeoLab command :

    `$ ProjectAtlasGeoLab -i input${format} -a atlasDir -ref mni_icbm152_t1_tal_nlin_asym_09c_brain.nii -o outputDir -cc clientComputeCentroids.py --rb clientRegisterBundles.py -ods dipyServer.py -cds clientCloseServer.py -nbPoints 15 -an Neigborhood${format2} -anc Centroids${format2} -nbThreads ${nbThreads}` 

* Replace ${format} with {.trk, .tck, .bundles} according to your tractogram format.
* Replace ${format2} with {Trk, Tck, Bundles} according to your tractogram format.
* mni_icbm152_t1_tal_nlin_asym_09c_brain.nii : path to the reference image mni_icbm152_t1_tal_nlin_asym_09c_brain.nii
* outDir : directory where to save the results.
* clientComputeCentroids.py : found in dipyServiceClient folder.
* clientRegisterBundles.py : found in dipyServiceClient folder.
* dipyServer.py : found in dipyServiceClient folder.
* clientCloseServer.py : found in ./dipyServiceClient folder.
* Neigborhood${format2} : found in ./Atlas/*.zip.
* Centroids${format2} : found in ./Atlas/*.zip.
* ${nbThreads} : number of threads to use for OpenMP.
 

## Usage example on other atlas

If yout atlas is in .tck/.trk format you'll need to create a ".minf" file for each of the ".tck/.trk" files in your atlas. You can use the createMinf command as in the example ***Usage example on ESBA atlas***.

You'll need to analyse your atlas to get the bundle-specific thresholds :

    `$ analyseAtlasBundle.py -i atlasDir -f ${format} -r mni_icbm152_t1_tal_nlin_asym_09c_brain.nii`

* atlasDir : Path to your atlas directory.
* Replace ${format} with {.trk, .tck, .bundles} according to your tractogram format.
* mni_icbm152_t1_tal_nlin_asym_09c_brain.nii : path to the reference image mni_icbm152_t1_tal_nlin_asym_09c_brain.nii





