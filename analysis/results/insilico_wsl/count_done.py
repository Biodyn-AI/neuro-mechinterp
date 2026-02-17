import os
base = "/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl"
done = []
empty = []
for d in sorted(os.listdir(base)):
    if d.startswith("perturb_"):
        gene = d.replace("perturb_", "")
        path = os.path.join(base, d)
        files = [f for f in os.listdir(path) if f.endswith('.pickle')]
        if files:
            done.append(gene)
        else:
            empty.append(gene)
print(f"Done ({len(done)}): {', '.join(done)}")
print(f"Empty ({len(empty)}): {', '.join(empty)}")
