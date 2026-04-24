
with open("../src/bundlesData.cxx", "r") as f:
    data = f.read()

# 1. Fix computeRotationMatrixFromVectorAndAngle signature in .cxx
old_crm = """void BundlesData::computeRotationMatrixFromVectorAndAngle(
                                      const std::array<float, 3>& vector,
                                      float angle, // in rad
                                      std::vector<float>& rotationMatrix ) const"""

new_crm = """void BundlesData::computeRotationMatrixFromVectorAndAngle(
                                      const std::array<float, 3>& vector,
                                      float angle, // in rad
                                      std::array<float, 9>& rotationMatrix ) const"""

data = data.replace(old_crm, new_crm)

# 2. Fix local rotationMatrix variables in functions
data = data.replace("std::vector<float> rotationMatrixSamePlane( 9, 0 ) ;", "std::array<float, 9> rotationMatrixSamePlane{0} ;")
data = data.replace("std::vector<float> rotationMatrixSamePlaneTmp( 9, 0 ) ;", "std::array<float, 9> rotationMatrixSamePlaneTmp{0} ;")
data = data.replace("std::vector<float> rotationMatrixSameDirection( 9, 0 ) ;", "std::array<float, 9> rotationMatrixSameDirection{0} ;")
data = data.replace("std::vector<float> rotationMatrixSameDirectionTmp( 9, 0 ) ;", "std::array<float, 9> rotationMatrixSameDirectionTmp{0} ;")

# Also other local vector<float> variables of size 9 used for rotation matrices
# we can just use simple regex for local declarations of size 9
import re
data = re.sub(r'std::vector\s*<\s*float\s*>\s+(\w+RotationMatrix\w*)\s*\(\s*9\s*,\s*0\s*\)\s*;', r'std::array<float, 9> \1{0} ;', data)
# Specific for rotationMatrixSamePlane
data = re.sub(r'std::vector\s*<\s*float\s*>\s+rotationMatrixSamePlane\s*\(\s*9\s*,\s*0\s*\)\s*;', r'std::array<float, 9> rotationMatrixSamePlane{0} ;', data)

with open("../src/bundlesData.cxx", "w") as f:
    f.write(data)

