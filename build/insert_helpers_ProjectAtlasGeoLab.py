import re

with open("../apps/ProjectAtlasGeoLab.cxx", "r") as f:
    content = f.read()

# Core helper logic from computeNeighborhood.cxx
helper_code = r"""
////////////////////////////////////////////////////////////////////////////////
//////////////////// Save index in tractogram of neighbors /////////////////////
////////////////////////////////////////////////////////////////////////////////
void saveIndexInTractogram( const char* labelsBinaryFilename,
                            const std::vector<int64_t>& labels )
{
  std::ofstream file ;
  file.open( labelsBinaryFilename, std::ios::binary | std::ios::out ) ;
  if ( file.fail() )
  {
    std::cout << "Problem opening file to write : " << labelsBinaryFilename << std::endl ;
    exit( 1 ) ;
  }
  int n_curves = labels.size() ;
  for ( int track = 0 ; track < n_curves ; track++ )
  {
    int64_t label = labels[ track ] ;
    file.write( reinterpret_cast<char*>( &label ), sizeof( int64_t ) ) ;
  }
  file.close() ;
}

////////////////////////////////////////////////////////////////////////////////
//////////////// Native In-Memory computeNeighborhood implementation ///////////
////////////////////////////////////////////////////////////////////////////////
void computeNeighborhoodNative( const std::string& inputTractogramPath,
                                const std::string& atlasDirectory,
                                const std::string& outputDirectory,
                                float minLength,
                                float maxLength,
                                float toleranceThr,
                                bool useDefaultMaxLength,
                                int nbThreadsCN,
                                int verbose )
{
    if (verbose) std::cout << "Native computeNeighborhood on: " << inputTractogramPath << std::endl;
    
    BundlesMinf inputTractogramInfo( inputTractogramPath.c_str() ) ;
    BundlesData inputTractogram( inputTractogramPath.c_str() ) ;

    std::string format = ".bundles";
    if (endswith(inputTractogramPath, ".trk")) format = ".trk";
    if (endswith(inputTractogramPath, ".tck")) format = ".tck";

    AtlasBundles atlasBundles( atlasDirectory.c_str(),
                               inputTractogram.isBundles,
                               inputTractogram.isTrk,
                               inputTractogram.isTck,
                               verbose ) ;

    int nbBundlesAtlas = atlasBundles.bundlesMinf.size() ;
    int nbThreadsUsed = nbThreadsCN > 0 ? nbThreadsCN : omp_get_max_threads();

    #pragma omp parallel for num_threads( nbThreadsUsed ) schedule(dynamic)
    for ( int bundle = 0 ; bundle < nbBundlesAtlas ; bundle++ )
    {
        BundlesMinf& atlasBundleInfo = atlasBundles.bundlesMinf[ bundle ] ;
        std::vector<float> medialPointAtlasBundle = atlasBundleInfo.centerBundle ;

        float thresholdDistanceBundle = atlasBundleInfo.maxRadius * ( 1 + toleranceThr ) ;
        if ( thresholdDistanceBundle == 0 ) thresholdDistanceBundle = 10 ; // fallback

        float maxLengthBundle = maxLength ;
        if ( useDefaultMaxLength ) maxLengthBundle = atlasBundleInfo.maxLength * 1.2 ;

        int64_t nbFibersTractogram = inputTractogram.curves_count ;
        int nbPoints = inputTractogram.pointsPerTrack[ 0 ] ;

        std::vector<int64_t> indexNeighborFibers ;
        int curveCountNeighborhood = 0 ;
        int64_t nbElementsExtractedNeighborhood = 0 ;

        bool isOK = false ;
        while ( !isOK )
        {
            indexNeighborFibers.clear();
            curveCountNeighborhood = 0 ;
            nbElementsExtractedNeighborhood = 0 ;

            for ( int fiberIndex = 0 ; fiberIndex < nbFibersTractogram ; fiberIndex++ )
            {
                std::array<float, 3> medialPointFiberTractogram{0, 0, 0} ;
                inputTractogram.computeMedialPointFiberWithDistance( fiberIndex, medialPointFiberTractogram ) ;

                float distance = 0 ;
                for ( int i = 0 ; i < 3 ; i++ ) {
                    distance += pow( medialPointFiberTractogram[ i ] - medialPointAtlasBundle[ i ], 2 ) ;
                }
                distance = sqrt( distance ) ;

                float lengthFiber = inputTractogram.computeLengthFiber( fiberIndex ) ;

                if ( distance < thresholdDistanceBundle &&
                     lengthFiber > minLength && lengthFiber < maxLengthBundle )
                {
                    indexNeighborFibers.push_back( fiberIndex ) ;
                    nbElementsExtractedNeighborhood += 3 * nbPoints ;
                    curveCountNeighborhood += 1 ;
                }
            }

            if ( curveCountNeighborhood > 20 ) { // matching minNbCurvesNeighborhood default
                isOK = true ;
            } else {
                thresholdDistanceBundle += 5 ;
            }
        }

        std::vector<float> extractedNeighborhoodTotal( nbElementsExtractedNeighborhood, 0 ) ;
        std::vector<int32_t> pointsPerTrackNeighborhood ;
        int32_t fiberIndexNeighborhood = 0 ;

        for ( int selectedFiberIndex = 0 ; selectedFiberIndex < indexNeighborFibers.size() ; selectedFiberIndex++ )
        {
            int fiberIndex = indexNeighborFibers[ selectedFiberIndex ] ;
            int64_t offsetTractogram = 3 * nbPoints * fiberIndex ;
            int offsetExtractedNeighborhood = 3 * nbPoints * fiberIndexNeighborhood ;

            pointsPerTrackNeighborhood.push_back( nbPoints ) ;
            std::copy( inputTractogram.matrixTracks.begin() + offsetTractogram,
                       inputTractogram.matrixTracks.begin() + offsetTractogram + 3 * nbPoints, 
                       extractedNeighborhoodTotal.begin() + offsetExtractedNeighborhood ) ;

            fiberIndexNeighborhood += 1 ;
        }

        BundlesMinf extractedNeighborhoodInfo( inputTractogramPath.c_str() ) ;
        std::string outFormat = extractedNeighborhoodInfo.getFormat() ;
        std::string outputBundlesFilename = outputDirectory + atlasBundles.bundlesNames[ bundle ] + outFormat ;

        extractedNeighborhoodInfo.curves_count = curveCountNeighborhood ;
        if ( outFormat == ".tck" || outFormat == ".trk" ) {
             extractedNeighborhoodInfo.haveMinf = false ;
        }

        std::vector<int64_t> fibersWithNans;
        std::vector<std::vector<float>> tracksScalars, tracksProperties;
        
        // Setup proper IsBundles/Trk/Tck flags
        bool isBundles = (outFormat == ".bundles" || outFormat == ".bundlesdata");
        bool isTrk = (outFormat == ".trk");
        bool isTck = (outFormat == ".tck");

        BundlesData extractedNeighborhoodData( extractedNeighborhoodTotal,
                                              pointsPerTrackNeighborhood,
                                              fibersWithNans,
                                              tracksScalars,
                                              tracksProperties,
                                              curveCountNeighborhood,
                                              isBundles, isTrk, isTck ) ;
                                              
        extractedNeighborhoodData.write( outputBundlesFilename.c_str(), extractedNeighborhoodInfo ) ;

        // Save index bin
        std::string labelsBinaryFilename = replaceExtension( outputBundlesFilename, "Index.bin" ) ;
        saveIndexInTractogram( labelsBinaryFilename.c_str(), indexNeighborFibers ) ;
    }
}
"""

# Insert before function applyGeoLab roughly line 55
insert_pattern = r'(// Function to get flag position when parsing arguments.*?\n\})'
match = re.search(insert_pattern, content, re.DOTALL)
if match:
    # Adding before applyGeoLab
    content_new = content.replace(match.group(1), match.group(1) + "\n" + helper_code)
    with open("../apps/ProjectAtlasGeoLab.cxx", "w") as f:
         f.write(content_new)
    print("Helpers inserted correctly.")
else:
    print("Could not insert helpers!")

