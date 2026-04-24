import re
import glob

files = [
    "../apps/alignFibers.cxx",
    "../apps/analyseAtlasBundle.cxx",
    "../apps/alignFibersInBundle.cxx",
    "../apps/alignFibersInBundleToPlaneXY.cxx"
]

for f_path in files:
    with open(f_path, "r") as f:
        content = f.read()
    
    # 1. replace std::vector<float> declarations of size 3 into std::array<float, 3>
    # Variables usually named normalVector, directionVector, medialPoint, gravityCenter, etc.
    content = re.sub(r'std::vector\s*<\s*float\s*>\s+(normalVector\w*)\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> \1{0, 0, 0} ;', content)
    content = re.sub(r'std::vector\s*<\s*float\s*>\s+(directionVector\w*)\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> \1{0, 0, 0} ;', content)
    content = re.sub(r'std::vector\s*<\s*float\s*>\s+(medialPoint\w*)\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> \1{0, 0, 0} ;', content)
    content = re.sub(r'std::vector\s*<\s*float\s*>\s+(gravityCenter\w*)\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> \1{0, 0, 0} ;', content)
    content = re.sub(r'std::vector\s*<\s*float\s*>\s+(vector\d+)\s*\(\s*3\s*,\s*0\s*\)\s*;', r'std::array<float, 3> \1{0, 0, 0} ;', content)

    # 2. Check for manual corrections like gravityCenterFiber declared as size 3
    content = content.replace("std::vector<float> gravityCenterFiber( 3, 0 ) ;", "std::array<float, 3> gravityCenterFiber{0, 0, 0} ;")

    with open(f_path, "w") as f:
        f.write(content)

