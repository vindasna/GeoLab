import re

with open("../src/bundlesData.h", "r") as f:
    content = f.read()

# 1. Add methods to BundlesData class
methods_to_add = """
  void resampleAll(int targetPoints);
  void registerFast(BundlesData& target, float qbThreshold = 15.0f, int iterations = 10, float outlierRejectionPct = 0.20f);

} ;"""

if "void resampleAll(int targetPoints)" not in content:
    content = content.replace("} ;", methods_to_add)

with open("../src/bundlesData.h", "w") as f:
    f.write(content)

# 2. Fix default arguments in bundlesData.cxx
with open("../src/bundlesData.cxx", "r") as f:
    content = f.read()

content = content.replace(
    "void BundlesData::registerFast(BundlesData& target, float qbThreshold = 15.0f, int iterations = 10, float outlierRejectionPct = 0.20f)",
    "void BundlesData::registerFast(BundlesData& target, float qbThreshold, int iterations, float outlierRejectionPct)"
)

with open("../src/bundlesData.cxx", "w") as f:
    f.write(content)

