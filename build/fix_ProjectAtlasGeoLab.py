import re

with open("../apps/ProjectAtlasGeoLab.cxx", "r") as f:
    content = f.read()

# Define the pattern to replace inside applyGeoLab
pattern = r'(\s*)// Computing centroids neighborhood atlas.*?\n(\s*)// Projection'

replacement = r"""\1// Native registration using registerFast (GeoLab)
\1std::string neighborhoodRegistered = replaceExtension(
\1                                                  movedTractogramNeighborhood,
\1                                                  "_moved.bundles" ) ;
\1if ( ( format == ".trk" || format == ".tck" ) )
\1{
\1  neighborhoodRegistered = replaceExtension( neighborhoodRegistered, format) ;
\1}

\1if ( is_file( neighborhoodRegistered ) && !force )
\1{
\1    if ( verbose > 1 ) {
\1         std::cout << "Using existing registered file: " << neighborhoodRegistered << std::endl;
\1    }
\1}
\1else
\1{
\1    if ( verbose > 1 ) {
\1         std::cout << "Starting native registerFast for bundle " << bundleName << "..." << std::endl;
\1    }

\1    try {
\1        // Load atlas and moving neighborhoods into memory
\1        BundlesData atlasNeighborhoodBundles( atlasNeighborhoodFile.c_str() ) ;
\1        BundlesData movingNeighborhoodBundles( movedTractogramNeighborhood.c_str() ) ;

\1        movingNeighborhoodBundles.registerFast( atlasNeighborhoodBundles, TransformType::RIGID ) ;

\1        // Save aligned bundle directly
\1        BundlesMinf movingInfo( movedTractogramNeighborhood.c_str() ) ;
\1        movingNeighborhoodBundles.write( neighborhoodRegistered.c_str(), movingInfo ) ;
\1    }
\1    catch (const std::exception& e) {
\1         // Gracefully handle if calculation crashes
\1         std::cerr << "ERROR in native registerFast for bundle " << bundleName << ": " << e.what() << std::endl;
\1         return; 
\1    }
\1}

\1// Projection"""

match = re.search(pattern, content, re.DOTALL)
if match:
    print("Found match! Replacing...")
    content_new = re.sub(pattern, replacement, content, count=1, flags=re.DOTALL)
    with open("../apps/ProjectAtlasGeoLab.cxx", "w") as f:
        f.write(content_new)
    print("Replacement applied successfully.")
else:
    print("Could not find native registration block based on comments.")

