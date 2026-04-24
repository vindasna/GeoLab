
with open("../apps/getFeaturesFibersAtlas.cxx", "r") as f:
    content = f.read()

content = content.replace("std::vector<float> xVector = {1, 0, 0 } ;", "std::array<float, 3> xVector = {1, 0, 0 } ;")

with open("../apps/getFeaturesFibersAtlas.cxx", "w") as f:
    f.write(content)

