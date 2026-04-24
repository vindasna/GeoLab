import re

with open("../src/AtlasBundles.cxx", "r") as f:
    content = f.read()

# Replace declarations of size 3 vectors for normalVector and directionVector
content = re.sub(r'std::vector\s*<\s*float\s*>\s+normalVector\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> normalVector{0, 0, 0} ;', content)
content = re.sub(r'std::vector\s*<\s*float\s*>\s+directionVector\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> directionVector{0, 0, 0} ;', content)

# Check for others like normalVectorFiberAtlas and medialPointFiberAtlas which are also used on lines 572-622
content = re.sub(r'std::vector\s*<\s*float\s*>\s+normalVectorFiberAtlas\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> normalVectorFiberAtlas{0, 0, 0} ;', content)
content = re.sub(r'std::vector\s*<\s*float\s*>\s+medialPointFiberAtlas\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> medialPointFiberAtlas{0, 0, 0} ;', content)

with open("../src/AtlasBundles.cxx", "w") as f:
    f.write(content)

