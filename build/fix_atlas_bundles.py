import re

with open("../src/AtlasBundles.cxx", "r") as f:
    content = f.read()

# Replace declarations of size 3 vectors for medial points
content = re.sub(r'std::vector\s*<\s*float\s*>\s+(medialPoint\w+)\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> \1{0, 0, 0} ;', content)

with open("../src/AtlasBundles.cxx", "w") as f:
    f.write(content)

