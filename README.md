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



