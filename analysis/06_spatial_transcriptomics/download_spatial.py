"""Download spatial transcriptomics datasets with progress tracking."""
import os
import sys
from urllib.request import urlopen


def download_file(url, out_path, chunk_size=1024*1024):
    """Download file with progress bar."""
    print(f"Downloading: {url}")
    print(f"To: {out_path}")
    
    req = urlopen(url, timeout=120)
    total_size = req.headers.get('Content-Length')
    if total_size:
        total_size = int(total_size)
        print(f"Size: {total_size / (1024*1024):.1f} MB")
    else:
        print("Size: unknown")
    
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    downloaded = 0
    
    with open(out_path, 'wb') as f:
        while True:
            chunk = req.read(chunk_size)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            if total_size:
                pct = downloaded / total_size * 100
                mb = downloaded / (1024*1024)
                print(f"\r  Progress: {mb:.1f} MB / {total_size/(1024*1024):.1f} MB ({pct:.1f}%)", end='', flush=True)
            else:
                print(f"\r  Downloaded: {downloaded / (1024*1024):.1f} MB", end='', flush=True)
    
    print(f"\nDone: {out_path}")


if __name__ == '__main__':
    # GSE233293 PDAC Visium (216MB)
    url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233293/suppl/GSE233293_RAW.tar'
    out = 'data/spatial/pdac_visium/GSE233293_RAW.tar'
    download_file(url, out)
