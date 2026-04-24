with open("../apps/ProjectAtlasGeoLab.cxx", "r") as f:
    content = f.read()

open_count = content.count('{')
close_count = content.count('}')

print(f"Open braces: {open_count}")
print(f"Close braces: {close_count}")
print(f"Difference (Open - Close): {open_count - close_count}")

# Let's also do a quick stack parse to find where it breaks!
stack = []
for i, char in enumerate(content):
    if char == '{':
        stack.append(i)
    elif char == '}':
        if stack:
            stack.pop()
        else:
            print(f"Extra closing brace at index {i}")

if stack:
    print(f"Open braces left unclosed: {len(stack)}")
    # Print the last unclosed brace line
    sub = content[:stack[-1]]
    line = sub.count('\n') + 1
    print(f"Last unclosed brace around line {line}")

