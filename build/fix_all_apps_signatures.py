import re
import glob

files = glob.glob("../apps/*.cxx")

for f_path in files:
    with open(f_path, "r") as f:
        content = f.read()
    content_changed = False

    # 1. Update computeAverageFiberBundle signature
    if "void computeAverageFiberBundle(" in content:
        old_cafb = """void computeAverageFiberBundle(
                        BundlesData& atlasBundleData,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        int nbPoints,
                        std::vector<float>& averageFiber,
                        std::vector<float>& medialPointAtlasBundle )"""
        new_cafb = """void computeAverageFiberBundle(
                        BundlesData& atlasBundleData,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        int nbPoints,
                        std::vector<float>& averageFiber,
                        std::array<float, 3>& medialPointAtlasBundle )"""
        if old_cafb in content:
            content = content.replace(old_cafb, new_cafb)
            content_changed = True

    # 2. Update computeGravityCenterAtlasBundle signature
    if "void computeGravityCenterAtlasBundle(" in content:
        old_cgca = """void computeGravityCenterAtlasBundle(
                                  BundlesData& atlasBundleData,
                                  int nbPoints,
                                  std::vector<float>& gravityCenterAtlasBundle )"""
        new_cgca = """void computeGravityCenterAtlasBundle(
                                  BundlesData& atlasBundleData,
                                  int nbPoints,
                                  std::array<float, 3>& gravityCenterAtlasBundle )"""
        if old_cgca in content:
            content = content.replace(old_cgca, new_cgca)
            content_changed = True

    # 3. Local variables called xVector or similar in getFeaturesFibersAtlas.cxx
    if "xVector" in content and "getFeaturesFibersAtlas" in f_path:
        content = content.replace("std::vector<float> xVector( 3, 0 ) ;", "std::array<float, 3> xVector{0, 0, 0} ;")
        content = content.replace("std::vector<float> yVector( 3, 0 ) ;", "std::array<float, 3> yVector{0, 0, 0} ;")
        content = content.replace("std::vector<float> zVector( 3, 0 ) ;", "std::array<float, 3> zVector{0, 0, 0} ;")
        content_changed = True

    # 4. Also alignFibersInBundleToPlaneXY.cxx local variables that might be declared differently
    if "alignFibersInBundleToPlaneXY" in f_path:
        content = content.replace("std::vector<float> normalVectorReference( 3, 0 ) ;", "std::array<float, 3> normalVectorReference{0, 0, 0} ;")
        content = content.replace("std::vector<float> normalVectorRegistered( 3, 0 ) ;", "std::array<float, 3> normalVectorRegistered{0, 0, 0} ;")
        content = content.replace("std::vector<float> directionVectorReference( 3, 0 ) ;", "std::array<float, 3> directionVectorReference{0, 0, 0} ;")
        content = content.replace("std::vector<float> directionVectorRegistered( 3, 0 ) ;", "std::array<float, 3> directionVectorRegistered{0, 0, 0} ;")
        content_changed = True

    if content_changed:
        print(f"Updated signatures in {f_path}")
        with open(f_path, "w") as f:
            f.write(content)

