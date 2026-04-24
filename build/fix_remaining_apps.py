import re

# 1. Fix alignFibersInBundleToPlaneXY.cxx
with open("../apps/alignFibersInBundleToPlaneXY.cxx", "r") as f:
    content = f.read()

content = content.replace("std::vector<float> normalVectorFiber( 3, 0 ) ;", "std::array<float, 3> normalVectorFiber{0, 0, 0} ;")
content = content.replace("std::vector<float> normalVector1( 3, 0 ) ;", "std::array<float, 3> normalVector1{0, 0, 0} ;")

with open("../apps/alignFibersInBundleToPlaneXY.cxx", "w") as f:
    f.write(content)

# 2. Fix analyseAtlasBundle.cxx Invalid reference passing
with open("../apps/analyseAtlasBundle.cxx", "r") as f:
    content = f.read()

# Lines 1085 were failing on something about invalid initialization
# Let's read the lines around 1085 using regex or view_file to make exact replacement!

# I will just grep for it in the script or run it iteratively.
# Since line 1085 had "gravityCenterAtlasBundle" or similar.
content = content.replace("std::vector<float> gravityCenterAtlasBundle( 3, 0 ) ;", "std::array<float, 3> gravityCenterAtlasBundle{0, 0, 0} ;")
content = content.replace("std::vector<float> medialPointAtlasBundle( 3, 0 ) ;", "std::array<float, 3> medialPointAtlasBundle{0, 0, 0} ;")

with open("../apps/analyseAtlasBundle.cxx", "w") as f:
    f.write(content)

