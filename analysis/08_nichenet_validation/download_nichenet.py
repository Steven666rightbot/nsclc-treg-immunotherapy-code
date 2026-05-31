import urllib.request
import os
import sys

url = 'https://zenodo.org/record/3260758/files/ligand_target_matrix.rds'
out = 'D:/Research/tomato/data/nichenet_ligand_target_matrix.rds'
os.makedirs(os.path.dirname(out), exist_ok=True)

print(f'Downloading from {url}')
print(f'Saving to {out}')

req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})

with urllib.request.urlopen(req, timeout=120) as response:
    total = int(response.headers.get('content-length', 0))
    print(f'Total size: {total / (1024*1024):.1f} MB')
    
    downloaded = 0
    chunk_size = 256 * 1024
    with open(out, 'wb') as f:
        while True:
            chunk = response.read(chunk_size)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            if downloaded % (5 * 1024 * 1024) == 0:
                pct = downloaded / total * 100 if total else 0
                print(f'  Progress: {downloaded/(1024*1024):.1f} MB ({pct:.1f}%)')

size = os.path.getsize(out) / (1024*1024)
print(f'Download complete: {size:.1f} MB')
