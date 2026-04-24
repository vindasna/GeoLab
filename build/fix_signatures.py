
with open("../src/bundlesData.cxx", "r") as f:
    data = f.read()

# 1. Fix translationVector local declaration
data = data.replace("std::vector<float> translationVector( 3 , 0 ) ;", "std::array<float, 3> translationVector{0, 0, 0} ;")

# 2. Fix applyRotationMatrixToVector signature block
old_armv = """void BundlesData::applyRotationMatrixToVector(
                                       const std::vector<float>& vector,
                                       const std::vector<float>& rotationMatrix,
                                       std::vector<float>& rotatedVector ) const"""

new_armv = """void BundlesData::applyRotationMatrixToVector(
                                       const std::array<float, 3>& vector,
                                       const std::array<float, 9>& rotationMatrix,
                                       std::array<float, 3>& rotatedVector ) const"""

data = data.replace(old_armv, new_armv)

# 3. Fix applyRotationMatrixToFiber overloads having std::vector<float>& rotationMatrix
# Since there are multiple overloads, let's replace all standard signatures
data = data.replace("const std::vector<float>& rotationMatrix", "const std::array<float, 9>& rotationMatrix")

# 4. Fix putFibersInSamePlane signature block (Overload at line 3701)
old_pfisp = """void BundlesData::putFibersInSamePlane(
                               const std::vector<float>& normalVector1,
                               const std::vector<float>& normalVector2,
                               const std::vector<float>& fiber2,
                               int nbPoints,
                               std::vector<float>& fiber2ToPlane1,
                               std::vector<float>& newNormalVectorFiber2 ) const"""

new_pfisp = """void BundlesData::putFibersInSamePlane(
                               const std::array<float, 3>& normalVector1,
                               const std::array<float, 3>& normalVector2,
                               const std::vector<float>& fiber2,
                               int nbPoints,
                               std::vector<float>& fiber2ToPlane1,
                               std::array<float, 3>& newNormalVectorFiber2 ) const"""

data = data.replace(old_pfisp, new_pfisp)

with open("../src/bundlesData.cxx", "w") as f:
    f.write(data)

