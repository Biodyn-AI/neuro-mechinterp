from huggingface_hub import hf_hub_download
import shutil, os

dest = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/")
repo = "ctheodoris/Geneformer"

# V1 (30M) dictionaries
files = {
    "geneformer/gene_dictionaries_30m/gene_median_dictionary_gc30M.pkl": "gene_median_dictionary.pkl",
    "geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl": "token_dictionary.pkl",
    "geneformer/gene_dictionaries_30m/gene_name_id_dict_gc30M.pkl": "gene_name_id_dict.pkl",
    "geneformer/gene_dictionaries_30m/ensembl_mapping_dict_gc30M.pkl": "ensembl_mapping_dict.pkl",
}

for remote_path, local_name in files.items():
    print(f"Downloading {local_name} from {remote_path}...")
    try:
        path = hf_hub_download(repo_id=repo, filename=remote_path)
        dest_path = os.path.join(dest, local_name)
        shutil.copy(path, dest_path)
        print(f"  Saved: {os.path.getsize(dest_path)} bytes")
    except Exception as e:
        print(f"  ERROR: {e}")

print("Done!")
