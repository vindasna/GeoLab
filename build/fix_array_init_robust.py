import re

with open("../src/bundlesData.cxx", "r") as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    orig_line = line
    if "std::array<float, 3>" in line and "( 3, 0 )" in line:
        line = line.replace("( 3, 0 )", "{0, 0, 0}")
    elif "std::array<float, 3>" in line and "(3, 0)" in line:
        line = line.replace("(3, 0)", "{0, 0, 0}")
    elif "std::array<float, 3>" in line and "({ 0, 0, 0 })" in line:
        line = line.replace("({ 0, 0, 0 })", "{0, 0, 0}")
    elif "std::array<float, 9>" in line and "( 9, 0 )" in line:
        line = line.replace("( 9, 0 )", "{0}")
    elif "std::array<float, 9>" in line and "(9, 0)" in line:
        line = line.replace("(9, 0)", "{0}")
    else:
        # Fallback regex for more spaces
        line = re.sub(r'std::array\s*<\s*float\s*,\s*3\s*>\s+(\w+)\s*\(\s*3\s*,\s*0\s*\)', r'std::array<float, 3> \1{0, 0, 0}', line)
        line = re.sub(r'std::array\s*<\s*float\s*,\s*9\s*>\s+(\w+)\s*\(\s*9\s*,\s*0\s*\)', r'std::array<float, 9> \1{0}', line)
    new_lines.append(line)

with open("../src/bundlesData.cxx", "w") as f:
    f.writelines(new_lines)

