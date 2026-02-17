from huggingface_hub import hf_hub_download
import shutil, os

dest = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/")
repo = "ctheodoris/Geneformer"
files = [
    "geneformer/gene_median_dictionary.pkl",
    "geneformer/token_dictionary.pkl", 
    "geneformer/gene_name_id_dict.pkl",
    "geneformer/ensembl_mapping_dict.pkl"
]
for f in files:
    name = f.split("/")[-1]
    print(f"Downloading {name}...")
    path = hf_hub_download(repo_id=repo, filename=f)
    shutil.copy(path, os.path.join(dest, name))
    print(f"  Saved to {dest}{name} ({os.path.getsize(os.path.join(dest, name))} bytes)")
print("Done!")
