
with open("../src/bundlesData.cxx", "r") as f:
    data = f.read()

# 1. Fix applyRotationMatrixToFiber parameter line 3372
data = data.replace("const std::array<float, 3>& fiber,", "const std::vector<float>& fiber,")

# 2. Fix rotationMatrixSamePlaneTmp local variable at line 3652
data = data.replace("std::vector<float> rotationMatrixSamePlaneTmp( 9, 0 );", "std::array<float, 9> rotationMatrixSamePlaneTmp{0};")

with open("../src/bundlesData.cxx", "w") as f:
    f.write(data)

