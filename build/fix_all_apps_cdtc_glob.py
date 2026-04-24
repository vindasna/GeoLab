import re
import glob

files = glob.glob("../apps/*.cxx")

for f_path in files:
    with open(f_path, "r") as f:
        content = f.read()
    content_changed = False

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
        content_changed = True
        print(f"Updated computeDistancesToCenterBundle signature in {f_path}")

    if content_changed:
        with open(f_path, "w") as f:
            f.write(content)

