import requests
import os
from tqdm import tqdm

# NCBI RefSeq: Yarrowia lipolytica CLIB122 (W29)
# 这是一个非常稳定的永久链接
URL = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/525/GCF_000002525.2_ASM252v1/GCF_000002525.2_ASM252v1_rna.fna.gz"
REF_DIR = "reference"
FILENAME = "Yarrowia_lipolytica.cdna.all.fa.gz" # 我们统一重命名为这个
FILE_PATH = os.path.join(REF_DIR, FILENAME)

def download_ref():
    if not os.path.exists(REF_DIR):
        os.makedirs(REF_DIR)
        
    print(f"[*] Downloading Reference from NCBI...")
    print(f"    Target: {FILENAME}")
    
    try:
        response = requests.get(URL, stream=True, timeout=30)
        response.raise_for_status() # 检查 404
        
        total_size = int(response.headers.get('content-length', 0))
        
        with open(FILE_PATH, "wb") as f:
            with tqdm(total=total_size, unit='B', unit_scale=True, unit_divisor=1024) as pbar:
                for data in response.iter_content(1024):
                    f.write(data)
                    pbar.update(len(data))
        
        print(f"[*] Download complete. File saved to {FILE_PATH}")
        
    except Exception as e:
        print(f"[!] Error: {e}")
        if os.path.exists(FILE_PATH):
            os.remove(FILE_PATH)

if __name__ == "__main__":
    download_ref()
