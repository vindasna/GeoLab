import re

with open("../src/bundlesData.cxx", "r") as f:
    content = f.read()

# 1. Remove default argument from registerTo definition
content = content.replace(
    "void BundlesData::registerTo(const BundlesData& target, int iterations = 10)",
    "void BundlesData::registerTo(const BundlesData& target, int iterations)"
)

# 2. Remove redefinition of struct Cluster
old_cluster_struct = """struct Cluster {
    std::vector<float> centroid; // Flattened array of 3 * numPoints
    int count;                   // Number of streamlines in this cluster
    std::vector<int> indices;    // Original indices of the streamlines
};"""

if old_cluster_struct in content:
    content = content.replace(old_cluster_struct, "")
else:
    # Try with slightly different spacing if it failed
    content = re.sub(
        r'struct Cluster\s*\{[^}]+\};',
        "",
        content
    )

with open("../src/bundlesData.cxx", "w") as f:
    f.write(content)

