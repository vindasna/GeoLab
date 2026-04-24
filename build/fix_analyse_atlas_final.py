import re

with open("../apps/analyseAtlasBundle.cxx", "r") as f:
    content = f.read()

# 1. Update computeAverageFiberBundle signature (Line 94)
content = content.replace(
    """void computeAverageFiberBundle(
                        BundlesData& atlasBundleData,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        int nbPoints,
                        std::vector<float>& averageFiber,
                        std::vector<float>& medialPointAtlasBundle )""",
    """void computeAverageFiberBundle(
                        BundlesData& atlasBundleData,
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        int nbPoints,
                        std::vector<float>& averageFiber,
                        std::array<float, 3>& medialPointAtlasBundle )"""
)

# 2. Update computeGravityCenterAtlasBundle signature (Line 209)
content = content.replace(
    """void computeGravityCenterAtlasBundle(
                                  BundlesData& atlasBundleData,
                                  int nbPoints,
                                  std::vector<float>& gravityCenterAtlasBundle )""",
    """void computeGravityCenterAtlasBundle(
                                  BundlesData& atlasBundleData,
                                  int nbPoints,
                                  std::array<float, 3>& gravityCenterAtlasBundle )"""
)

# 3. Update local newNormalVectorFiber variables that regular expressions might have missed
content = content.replace("std::vector<float> newNormalVectorFiber1( 3, 0 ) ;", "std::array<float, 3> newNormalVectorFiber1{0, 0, 0} ;")
content = content.replace("std::vector<float> newNormalVectorFiber2( 3, 0 ) ;", "std::array<float, 3> newNormalVectorFiber2{0, 0, 0} ;")

with open("../apps/analyseAtlasBundle.cxx", "w") as f:
    f.write(content)

