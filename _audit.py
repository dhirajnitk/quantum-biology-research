import re

with open('papers/p0_sub_tubulin/sub_tubular_quantum.tex', encoding='utf-8') as f:
    text = f.read()

# 1. Section structure
print('=== SECTIONS ===')
for m in re.finditer(r'\\(?:sub)*section\{([^}]+)\}', text):
    print(f'  {m.group(1)}')

# 2. Duplicate headings
headings = [m.group() for m in re.finditer(r'\\(?:sub)*section\{([^}]+)\}', text)]
dupes = [h for h in headings if headings.count(h) > 1]
print(f'\nDUPLICATE HEADINGS:')
print(f'  {set(dupes) if dupes else "NONE"}')

# 3. Citation integrity
cites = set()
for m in re.finditer(r'\\cite\{([^}]+)\}', text):
    for key in m.group(1).split(','):
        cites.add(key.strip())
bibs = set()
for m in re.finditer(r'\\bibitem\{(\w+)\}', text):
    bibs.add(m.group(1))
print(f'\nCITATION INTEGRITY:')
print(f'  Cites: {len(cites)}, Bibs: {len(bibs)}')
if bibs - cites:
    print(f'  ORPHANED: {bibs - cites}')
if cites - bibs:
    print(f'  MISSING: {cites - bibs}')

# 4. Stale values
stale = ['5,000', '5\,000', '7.84', 'per-core', 'E_bit', '16-fold', '16-fold', 'Landauer limit']
for pat in stale:
    if pat in text:
        print(f'\nFOUND STALE: "{pat}"')

# 5. Figures
for m in re.finditer(r'\\includegraphics\[([^\]]*)\]\{([^}]+)\}', text):
    print(f'\nFIGURE: {m.group(2)}')

# 6. Label/ref
labels = {m.group(1) for m in re.finditer(r'\\label\{([^}]+)\}', text)}
refs = {m.group(1) for m in re.finditer(r'\\ref\{([^}]+)\}', text)}
if refs - labels:
    print(f'\nBROKEN REFS: {refs - labels}')
if labels - refs:
    print(f'\nUNREFERENCED LABELS: {labels - refs}')

print('\n=== DONE ===')
