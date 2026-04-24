import re

with open("../src/bundlesData.h", "r") as f:
    content = f.read()

# Replace const std::vector<float>& rotationMatrix with robust regex for spaces
content = re.sub(r'const\s+std::vector\s*<\s*float\s*>\s*&\s+rotationMatrix', r'const std::array<float, 9>& rotationMatrix', content)
# For reference: std::vector<float>& rotationMatrix
content = re.sub(r'std::vector\s*<\s*float\s*>\s*&\s+rotationMatrix', r'std::array<float, 9>& rotationMatrix', content)
# Also fix medialPointFiber
content = re.sub(r'const\s+std::vector\s*<\s*float\s*>\s*&\s+medialPointFiber', r'const std::array<float, 3>& medialPointFiber', content)

with open("../src/bundlesData.h", "w") as f:
    f.write(content)

