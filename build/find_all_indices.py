with open("../apps/ProjectAtlasGeoLab.cxx", "r") as f:
    content = f.read()

indices = [16494, 67212, 69498, 96960, 97849]
for idx in indices:
    sub = content[:idx]
    line = sub.count('\n') + 1
    print(f"--- Index {idx} around line {line} ---")
    print(content[idx-100:idx+100])
    print("\n")

