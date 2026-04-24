with open("../apps/ProjectAtlasGeoLab.cxx", "r") as f:
    content = f.read()

indices = [16494, 67212, 69498, 96960, 97849]

# Verify they are all '}'
for idx in indices:
    if content[idx] != '}':
        print(f"Index {idx} is NOT a closing brace! It is '{content[idx]}'")
        exit(1)

# Remove in reverse order
content_new = content
for idx in sorted(indices, reverse=True):
    content_new = content_new[:idx] + content_new[idx+1:]

with open("../apps/ProjectAtlasGeoLab.cxx", "w") as f:
    f.write(content_new)

print("Removed 5 extra closing braces in reverse order.")

# Verify count again
open_count = content_new.count('{')
close_count = content_new.count('}')
print(f"New Open: {open_count}, New Close: {close_count}")

