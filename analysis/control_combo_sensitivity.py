#!/usr/bin/env python3
"""Sensitivity analysis: compare intel vs control combos excluding low-n pairs."""
import json, numpy as np
from scipy import stats

with open('/mnt/d/openclaw/intelligence-augmentation/analysis/CONTROL_COMBO_RESULTS.json') as f:
    data = json.load(f)

combos = data['control_combos']
# Filter by n >= 20
filtered = {k: v for k, v in combos.items() if v['n'] >= 20}
print(f"Control pairs with n>=20: {len(filtered)}/{len(combos)}")
for k, v in sorted(filtered.items(), key=lambda x: -x[1]['ratio']):
    print(f"  {k}: ratio={v['ratio']:.3f}, n={v['n']}")

filtered_ratios = [v['ratio'] for v in filtered.values()]
intel_ratios = data['intel_ratios']

print(f"\nFiltered control: mean={np.mean(filtered_ratios):.3f} ± {np.std(filtered_ratios):.3f}, n={len(filtered_ratios)}")
print(f"Intelligence: mean={np.mean(intel_ratios):.3f}, n={len(intel_ratios)}")

u, p = stats.mannwhitneyu(intel_ratios, filtered_ratios, alternative='greater')
print(f"Mann-Whitney U (intel > filtered control): U={u}, p={p:.6f}")

# Also n>=10
filtered10 = {k: v for k, v in combos.items() if v['n'] >= 10}
r10 = [v['ratio'] for v in filtered10.values()]
u10, p10 = stats.mannwhitneyu(intel_ratios, r10, alternative='greater')
print(f"\nn>=10: {len(filtered10)} pairs, mean={np.mean(r10):.3f}, U={u10}, p={p10:.6f}")
