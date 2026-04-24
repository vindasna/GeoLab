import re

with open("../apps/alignFibers.cxx", "r") as f:
    content = f.read()

# 1. Update computeFiberWithVectors Overload 1 (Line 59)
content = content.replace(
    """void computeFiberWithVectors( BundlesData& inputFiber,
                              BundlesMinf& inputFiberInfo,
                              std::vector<float> vector1,
                              std::vector<float> vector2,
                              std::vector<float> vector3,
                              int nbPoints,""",
    """void computeFiberWithVectors( BundlesData& inputFiber,
                              BundlesMinf& inputFiberInfo,
                              const std::array<float, 3>& vector1,
                              const std::array<float, 3>& vector2,
                              const std::array<float, 3>& vector3,
                              int nbPoints,"""
)

# 2. Update computeFiberWithVectors Overload 2 (Line 136)
content = content.replace(
    """void computeFiberWithVectors( BundlesData& inputFiber,
                              BundlesMinf& inputFiberInfo,
                              std::vector<float>& vector1,
                              std::vector<float>& vector2,
                              int nbPoints,""",
    """void computeFiberWithVectors( BundlesData& inputFiber,
                              BundlesMinf& inputFiberInfo,
                              const std::array<float, 3>& vector1,
                              const std::array<float, 3>& vector2,
                              int nbPoints,"""
)

# 3. Update local variables in main
content = content.replace("std::vector<float> newNormalVectorFiber1( 3, 0 ) ;", "std::array<float, 3> newNormalVectorFiber1{0, 0, 0} ;")
content = content.replace("std::vector<float> newNormalVectorFiber2( 3, 0 ) ;", "std::array<float, 3> newNormalVectorFiber2{0, 0, 0} ;")

with open("../apps/alignFibers.cxx", "w") as f:
    f.write(content)

