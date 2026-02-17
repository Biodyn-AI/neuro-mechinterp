import pickle, os
from datasets import load_from_disk

# Load tokenized data
ds = load_from_disk("/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl/brain.dataset")
print(f"Dataset: {len(ds)} cells")
print(f"First cell tokens (first 20): {ds[0]['input_ids'][:20]}")

# Load token dict
td_path = os.path.expanduser("~/.local/lib/python3.12/site-packages/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl")
with open(td_path, "rb") as f:
    token_dict = pickle.load(f)
reverse_dict = {v: k for k, v in token_dict.items()}

# Intelligence gene tokens
genes = {
    "MEF2C": 1532, "BDNF": 669, "GRIN2B": 24049, "CADM2": 13697,
    "NRXN1": 14485, "CAMK2A": 1135, "GRIN2A": 15103, "SHANK3": 20642,
    "HOMER1": 9482, "APP": 8172, "NEGR1": 13147, "NLGN1": 12574,
    "TCF4": 16602, "MAPT": 15882, "FOXO3": 4828, "CREB1": 4785,
    "FMR1": 2621, "SYN1": 201, "SCN1A": 8504, "SLC6A4": 3545, "COMT": 1996
}

# Check which genes appear in any cell
all_tokens = set()
for i in range(len(ds)):
    all_tokens.update(ds[i]['input_ids'])

print(f"\nTotal unique tokens in dataset: {len(all_tokens)}")
for name, tok in genes.items():
    present = tok in all_tokens
    print(f"  {name} (token {tok}): {'PRESENT' if present else 'MISSING'}")
