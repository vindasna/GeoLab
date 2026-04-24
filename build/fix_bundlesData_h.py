import re

with open("../src/bundlesData.h", "r") as f:
    content = f.read()

# 1. Add Cluster struct before class BundlesData
cluster_struct = """struct Cluster {
    std::vector<float> centroid; // Flattened array of 3 * numPoints
    int count;                   // Number of streamlines in this cluster
    std::vector<int> indices;    // Original indices of the streamlines
};

class BundlesData"""

if "struct Cluster" not in content:
    content = content.replace("class BundlesData", cluster_struct)

# 2. Add methods before the closing bracket of class BundlesData
# The class ends with "} ;" and then usually end of file or blank lines.
# I'll place them just before the final "} ;" at line 455 or similar.
methods_to_add = """
  void registerTo(const BundlesData& target, int iterations = 10);

  std::vector<Cluster> computeQuickBundles(float distance_threshold);

} ;"""

if "void registerTo(const BundlesData& target" not in content:
    content = content.replace("\n} ;", "\n" + methods_to_add)

with open("../src/bundlesData.h", "w") as f:
    f.write(content)

