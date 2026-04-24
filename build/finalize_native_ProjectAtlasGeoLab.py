import re

with open("../apps/ProjectAtlasGeoLab.cxx", "r") as f:
    content = f.read()

# 1. Replace Global SLR (approx lines 2336 to 2394)
slr_pattern = r'(\s*)std::ostringstream globalSLROss.*?\n(\s*)exit\( 1 \) ;'
slr_replacement = r"""\1int nationaleSLRfail = 0 ;
\1if ( is_file( movedTractogramTrk ) && !force )
\1{
\1    if ( verbose > 1 ) std::cout << "Using existing global SLR: " << movedTractogramTrk << std::endl;
\1}
\1else
\1{
\1    try {
\1        if (verbose) std::cout << "Starting native global registerFast..." << std::endl;
\1        BundlesData atlasFull( fullAtlasFilename.c_str() ) ;
\1        BundlesData movingFull( inputBundlesFilename.c_str() ) ;
\1        movingFull.registerFast( atlasFull, TransformType::RIGID ) ;
\1        BundlesMinf movingInfo( inputBundlesFilename.c_str() ) ;
\1        movingFull.write( movedTractogramTrk.c_str(), movingInfo ) ;
\1    } catch (const std::exception& e) {
\1        std::cerr << "ERROR in native global registerFast: " << e.what() << std::endl;
\1        nationaleSLRfail = 1;
\1    }
\1}
\1if ( nationaleSLRfail ) exit(1) ;"""

if re.search(slr_pattern, content, re.DOTALL):
    print("Found Global SLR pattern. Replacing...")
    content = re.sub(slr_pattern, slr_replacement, content, count=1, flags=re.DOTALL)

# 2. Replace computeNeighborhood subcalls
# Pattern 1 for Tractogram Neighborhood
neigh_pattern1 = r'(\s*)std::ostringstream computeNeighborhoodCommandOss.*?\n(\s*)std::string computeNeighborhoodCommand.*?\n(\s*)int isNeighborhoodFail = 0 ;.*?\n\s*if.*?isNeighborhoodFail = run_sh_process.*?\n.*?\n(\s*)if \( isNeighborhoodFail \).*?\n\s*exit\( 1 \) ;'
neigh_replacement1 = r"""\1computeNeighborhoodNative( movedTractogram, atlasDirectory, tmpNeighborhoodDir, 10, 200, toleranceThrComputeNeighborhood, true, nbThreadsCN, verbose ) ;"""

if re.search(neigh_pattern1, content, re.DOTALL):
    print("Found Tractogram Neighborhood pattern. Replacing...")
    content = re.sub(neigh_pattern1, neigh_replacement1, content, count=1, flags=re.DOTALL)

# Pattern 2 for Atlas Neighborhood
neigh_pattern2 = r'(\s*)std::ostringstream computeAtlasNeighborhoodCommandOss.*?\n.*?\n(\s*)if \( isAtlasNeighborhoodFail \).*?\n\s*exit\( 1 \) ;'
neigh_replacement2 = r"""\1computeNeighborhoodNative( fullAtlasFilename, atlasDirectory, tmpNeighborhoodAtlasDir, 10, 200, toleranceThrComputeNeighborhood, true, nbThreadsCN, verbose ) ;"""

if re.search(neigh_pattern2, content, re.DOTALL):
    print("Found Atlas Neighborhood pattern. Replacing...")
    content = re.sub(neigh_pattern2, neigh_replacement2, content, count=1, flags=re.DOTALL)

# 3. Clean Dipy Server Start/Stop
dipy_start_pattern = r'(\s*)// Launching dipy service.*?\n(\s*)while \( c\.running\(\) \)\s*\{\s*\n'
dipy_start_replacement = r"""\1// Launching native projection (dipy service removed)
\1{
"""

if re.search(dipy_start_pattern, content, re.DOTALL):
    print("Found Dipy Start pattern. Replacing...")
    content = re.sub(dipy_start_pattern, dipy_start_replacement, content, count=1, flags=re.DOTALL)

dipy_stop_pattern = r'(\s*)// Closing dipy service.*?\n(\s*)closeDipyServer\( portDipyServer \).*?\n.*?\n(\s*)break ;\s*\n\s*\}'
dipy_stop_replacement = r"""\1}"""

if re.search(dipy_stop_pattern, content, re.DOTALL):
    print("Found Dipy Stop pattern. Replacing...")
    content = re.sub(dipy_stop_pattern, dipy_stop_replacement, content, count=1, flags=re.DOTALL)

# 4. Total Duration displays
duration_pattern = r'(\s*)if \( keepTmpFiles \).*?\n(\s*)return\( 0 \) ;'
duration_replacement = r"""\1const std::chrono::duration< double > duration = std::chrono::system_clock::now() - start_time ;
\1std::cout << "\nTotal execution time: " << duration.count() << " seconds" << std::endl ;
\1return(0);"""

if re.search(duration_pattern, content, re.DOTALL):
    print("Found Duration pattern. Replacing...")
    content = re.sub(duration_pattern, duration_replacement, content, count=1, flags=re.DOTALL)

# 5. Clean closeDipyServer Function
close_dipy_pattern = r'(\s*)void closeDipyServer\( int portDipyServer \).*?\n(\s*)\}'
close_dipy_replacement = r""

if re.search(close_dipy_pattern, content, re.DOTALL):
    print("Found closeDipyServer function. Removing...")
    content = re.sub(close_dipy_pattern, close_dipy_replacement, content, count=1, flags=re.DOTALL)

with open("../apps/ProjectAtlasGeoLab.cxx", "w") as f:
    f.write(content)

print("Replacement complete.")

