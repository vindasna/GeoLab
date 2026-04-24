import re

# 1. Update bundlesData.h
with open("../src/bundlesData.h", "r") as f:
    h_content = f.read()

# Replace previous registerFast declaration with the new one
old_declaration = """Eigen::Matrix4f registerFast(BundlesData& target, 
                                  TransformType regType = TransformType::RIGID, 
                                  float qbThreshold = 10.0f, 
                                  int iterations = 10, 
                                  float outlierRejectionPct = 0.20f,
                                  BundlesData* outMovingCentroids = nullptr,
                                  BundlesData* outTargetCentroids = nullptr);"""

new_declaration = """Eigen::Matrix4f registerFast(BundlesData& target, 
                                  TransformType regType = TransformType::RIGID, 
                                  float qbThreshold = 10.0f, 
                                  int iterations = 10, 
                                  float outlierRejectionPct = 0.20f,
                                  BundlesData* outOriginalMovingCentroids = nullptr,
                                  BundlesData* outMovingCentroids = nullptr,
                                  BundlesData* outTargetCentroids = nullptr);"""

if old_declaration in h_content:
    h_content = h_content.replace(old_declaration, new_declaration)
else:
    # Try regex match in case formatting differs slightly
    h_content = re.sub(
        r'Eigen::Matrix4f\s+registerFast\s*\([^)]+\)\s*;',
        new_declaration,
        h_content
    )

with open("../src/bundlesData.h", "w") as f:
    f.write(h_content)


# 2. Update bundlesData.cxx
with open("../src/bundlesData.cxx", "r") as f:
    cxx_content = f.read()

# Replace function signature in .cxx
old_sig = """Eigen::Matrix4f BundlesData::registerFast(BundlesData& target, 
                                          TransformType regType, 
                                          float qbThreshold, 
                                          int iterations, 
                                          float outlierRejectionPct,
                                          BundlesData* outMovingCentroids,
                                          BundlesData* outTargetCentroids)"""

new_sig = """Eigen::Matrix4f BundlesData::registerFast(BundlesData& target, 
                                          TransformType regType, 
                                          float qbThreshold, 
                                          int iterations, 
                                          float outlierRejectionPct,
                                          BundlesData* outOriginalMovingCentroids,
                                          BundlesData* outMovingCentroids,
                                          BundlesData* outTargetCentroids)"""

if old_sig in cxx_content:
    cxx_content = cxx_content.replace(old_sig, new_sig)

# Add copy for original centroids coordinates FLAT
flat_coordinate_step = """    std::vector<float> targetCentroidsFlat(numTargetCentroids * N_POINTS * 3);
    for(int i = 0; i < numTargetCentroids; ++i) {
        std::copy(targetClusters[i].centroid.begin(), targetClusters[i].centroid.end(), targetCentroidsFlat.begin() + i * N_POINTS * 3);
    }"""

flat_coords_add = """    std::vector<float> targetCentroidsFlat(numTargetCentroids * N_POINTS * 3);
    for(int i = 0; i < numTargetCentroids; ++i) {
        std::copy(targetClusters[i].centroid.begin(), targetClusters[i].centroid.end(), targetCentroidsFlat.begin() + i * N_POINTS * 3);
    }

    // Save copy of ORIGINAL moving centroids BEFORE ICP loop
    std::vector<float> originalMovingCentroidsFlat = movingCentroidsFlat;"""

if flat_coordinate_step in cxx_content:
    cxx_content = cxx_content.replace(flat_coordinate_step, flat_coords_add)

# Add export step for original centroids at the end of registerFast
export_step = """    // --- STEP 5: Export Centroids (Optional) ---
    if (outMovingCentroids != nullptr)"""

export_add = """    // --- STEP 5: Export Centroids (Optional) ---
    if (outOriginalMovingCentroids != nullptr)
    {
        outOriginalMovingCentroids->curves_count = numMovingCentroids;
        outOriginalMovingCentroids->pointsPerTrack.assign(numMovingCentroids, N_POINTS);
        outOriginalMovingCentroids->matrixTracks = originalMovingCentroidsFlat; 
    }

    if (outMovingCentroids != nullptr)"""

if export_step in cxx_content:
    cxx_content = cxx_content.replace(export_step, export_add)

with open("../src/bundlesData.cxx", "w") as f:
    f.write(cxx_content)


# 3. Update apps/registerTractogramsFast.cxx
with open("../apps/registerTractogramsFast.cxx", "r") as f:
    app_content = f.read()

# Add originalMovingCentroids declaration
decl_p = """    BundlesData movingCentroids;
    BundlesData targetCentroids;"""
new_decl_p = """    BundlesData originalMovingCentroids;
    BundlesData movingCentroids;
    BundlesData targetCentroids;"""

if decl_p in app_content:
    app_content = app_content.replace(decl_p, new_decl_p)

# Update function call
call_p = """    Eigen::Matrix4f transformMatrix = movingBundles.registerFast( fixedBundles, regType, qbThr, nbIter, outlierRejectionThr, &movingCentroids, &targetCentroids ) ;"""
new_call_p = """    Eigen::Matrix4f transformMatrix = movingBundles.registerFast( fixedBundles, regType, qbThr, nbIter, outlierRejectionThr, &originalMovingCentroids, &movingCentroids, &targetCentroids ) ;"""

if call_p in app_content:
    app_content = app_content.replace(call_p, new_call_p)

# Update centroids saving logic to save ORIGINAL and MOVED supporting types
save_cent = """    // 2. Save Centroids
    std::cout << "Saving Centroids..." << std::endl;
    movingCentroids.isBundles = movingBundles.isBundles;
    movingCentroids.isTrk = movingBundles.isTrk;
    movingCentroids.isTck = movingBundles.isTck;

    targetCentroids.isBundles = fixedBundles.isBundles;
    targetCentroids.isTrk = fixedBundles.isTrk;
    targetCentroids.isTck = fixedBundles.isTck;

    std::string movingCentroidsPath = outputDir + "moving_centroids" + extMoving;
    std::string targetCentroidsPath = outputDir + "target_centroids" + extFixed;

    BundlesMinf movingCentroidsInfo = movingInfo;
    movingCentroidsInfo.curves_count = movingCentroids.curves_count;
    movingCentroids.write(movingCentroidsPath.c_str(), movingCentroidsInfo);

    BundlesMinf targetCentroidsInfo = fixedInfo;
    targetCentroidsInfo.curves_count = targetCentroids.curves_count;
    targetCentroids.write(targetCentroidsPath.c_str(), targetCentroidsInfo);"""

new_save_cent = """    // 2. Save Centroids
    std::cout << "Saving Centroids..." << std::endl;
    originalMovingCentroids.isBundles = movingBundles.isBundles;
    originalMovingCentroids.isTrk = movingBundles.isTrk;
    originalMovingCentroids.isTck = movingBundles.isTck;

    movingCentroids.isBundles = movingBundles.isBundles;
    movingCentroids.isTrk = movingBundles.isTrk;
    movingCentroids.isTck = movingBundles.isTck;

    targetCentroids.isBundles = fixedBundles.isBundles;
    targetCentroids.isTrk = fixedBundles.isTrk;
    targetCentroids.isTck = fixedBundles.isTck;

    std::string originalCentroidsPath = outputDir + "original_moving_centroids" + extMoving;
    std::string movedCentroidsPath = outputDir + "moved_centroids" + extMoving;
    std::string targetCentroidsPath = outputDir + "target_centroids" + extFixed;

    BundlesMinf movingCentroidsInfo = movingInfo;
    movingCentroidsInfo.curves_count = movingCentroids.curves_count;

    originalMovingCentroids.write(originalCentroidsPath.c_str(), movingCentroidsInfo);
    movingCentroids.write(movedCentroidsPath.c_str(), movingCentroidsInfo);

    BundlesMinf targetCentroidsInfo = fixedInfo;
    targetCentroidsInfo.curves_count = targetCentroids.curves_count;
    targetCentroids.write(targetCentroidsPath.c_str(), targetCentroidsInfo);"""

if save_cent in app_content:
    app_content = app_content.replace(save_cent, new_save_cent)

with open("../apps/registerTractogramsFast.cxx", "w") as f:
    f.write(app_content)

