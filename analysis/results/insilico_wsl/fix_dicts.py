import urllib.request
import os

dest = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/")
os.makedirs(dest, exist_ok=True)

base = "https://huggingface.co/ctheodoris/Geneformer/resolve/main/geneformer/gene_dictionaries_30m"
files = [
    "gene_median_dictionary_gc30M.pkl",
    "token_dictionary_gc30M.pkl",
    "gene_name_id_dict_gc30M.pkl",
    "ensembl_mapping_dict_gc30M.pkl",
]

for f in files:
    url = f"{base}/{f}"
    out = os.path.join(dest, f)
    print(f"Downloading {f}...")
    urllib.request.urlretrieve(url, out)
    size = os.path.getsize(out)
    # Check if it's an LFS pointer
    with open(out, 'rb') as fh:
        header = fh.read(10)
    if header.startswith(b'version'):
        print(f"  WARNING: Still LFS pointer! {size} bytes")
    else:
        print(f"  OK: {size} bytes")

print("Done!")
