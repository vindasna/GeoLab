import re

with open("../src/AtlasBundles.cxx", "r") as f:
    content = f.read()

# Replace declarations of size 3 vectors for newNormalVectorFiber2
content = re.sub(r'std::vector\s*<\s*float\s*>\s+newNormalVectorFiber2\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> newNormalVectorFiber2{0, 0, 0} ;', content)

# Also check for any variables with similar names
content = re.sub(r'std::vector\s*<\s*float\s*>\s+newNormalVectorTmp\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> newNormalVectorTmp{0, 0, 0} ;', content)

with open("../src/AtlasBundles.cxx", "w") as f:
    f.write(content)

