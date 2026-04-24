with open("../apps/ProjectAtlasGeoLab.cxx", "r") as f:
    content = f.read()

index = 16494
sub = content[:index]
line = sub.count('\n') + 1
print(f"Index 16494 is around line {line}")

# Let's print the local context around it
print("--- Context ---")
print(content[index-200:index+200])

