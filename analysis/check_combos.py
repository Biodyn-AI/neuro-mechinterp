import os
base = '/mnt/d/openclaw/intelligence-augmentation/analysis/results/insilico_wsl'
for d in sorted(os.listdir(base)):
    if d.startswith('perturb_combo_'):
        path = os.path.join(base, d)
        files = [f for f in os.listdir(path) if f.endswith('.pickle')]
        print(f"{d}: {len(files)} pickle files")
        for f in files:
            print(f"  {f}")
