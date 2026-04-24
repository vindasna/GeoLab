import re
import glob

# 1. Fix alignFibersInBundleToPlaneXY.cxx
with open("../apps/alignFibersInBundleToPlaneXY.cxx", "r") as f:
    content = f.read()

content = content.replace(
    "std::vector<float> normalVectorReference{ 0, 0, 1 } ;",
    "std::array<float, 3> normalVectorReference{ 0, 0, 1 } ;"
)
content = content.replace(
    "std::vector<float> directionVectorReference{ 1, 0, 0 } ;",
    "std::array<float, 3> directionVectorReference{ 1, 0, 0 } ;"
)

with open("../apps/alignFibersInBundleToPlaneXY.cxx", "w") as f:
    f.write(content)

# 2. Fix processAtlasInformation.cxx and computeCenterAtlasBundles.cxx
files = ["../apps/processAtlasInformation.cxx", "../apps/computeCenterAtlasBundles.cxx"]
for f_path in files:
    with open(f_path, "r") as f:
        content = f.read()
    
    # Update computeDistancesToCenterBundle signature
    old_cdtc = """void computeDistancesToCenterBundle(
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        const std::vector<float>& medialPointAtlasBundle,
                        int nbFibersAtlasBundle,
                        std::vector<float>& distancesToCenterAtlasBundle )"""
    new_cdtc = """void computeDistancesToCenterBundle(
                        const std::vector<float>& medialPointsAtlasBundleFibers,
                        const std::array<float, 3>& medialPointAtlasBundle,
                        int nbFibersAtlasBundle,
                        std::vector<float>& distancesToCenterAtlasBundle )"""
    
    if old_cdtc in content:
        content = content.replace(old_cdtc, new_cdtc)

    with open(f_path, "w") as f:
        f.write(content)

# 3. Fix getFeaturesFibersAtlas.cxx xVector that might have extra spaces
with open("../apps/getFeaturesFibersAtlas.cxx", "r") as f:
    content = f.read()

content = re.sub(r'std::vector\s*<\s*float\s*>\s+([xyz]Vector)\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> \1{0, 0, 0} ;', content)

with open("../apps/getFeaturesFibersAtlas.cxx", "w") as f:
    f.write(content)

