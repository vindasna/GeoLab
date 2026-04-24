import re
import glob

files = glob.glob("../apps/*.cxx")

for f_path in files:
    with open(f_path, "r") as f:
        content = f.read()
    
    # Robust declaration replacement of size 3 vectors
    # Matches: std::vector<float> name( 3 , 0 ) ;   or  (3,0) etc.
    new_content = re.sub(
        r'std::vector\s*<\s*float\s*>\s+([a-zA-Z0-9_\d]+)\s*\(\s*3\s*,\s*0\s*\)\s*;',
        r'std::array<float, 3> \1{0, 0, 0} ;',
        content
    )
    
    if new_content != content:
        print(f"Updated size-3 declarations in {f_path}")
        with open(f_path, "w") as f:
            f.write(new_content)

