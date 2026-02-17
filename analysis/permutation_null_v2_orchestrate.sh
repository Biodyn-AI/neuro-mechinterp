#!/bin/bash
# Orchestrate permutation null experiment: run one gene at a time in fresh subprocesses
set -e

ANALYSIS_DIR="/mnt/d/openclaw/intelligence-augmentation/analysis"
SELECTION_JSON="$ANALYSIS_DIR/permutation_null_v2_selection.json"
RESULTS_DIR="$ANALYSIS_DIR/results/insilico_wsl/permutation_null_v2"
INCREMENTAL_JSON="$ANALYSIS_DIR/permutation_null_v2_incremental.json"
RUNNER="$ANALYSIS_DIR/permutation_null_v2_run_one.py"

mkdir -p "$RESULTS_DIR"

# Initialize incremental results file if it doesn't exist
if [ ! -f "$INCREMENTAL_JSON" ]; then
    echo '{}' > "$INCREMENTAL_JSON"
fi

# Read selected controls from JSON and run each one
echo "=== PERMUTATION NULL V2: Running 50 control genes one at a time ==="
echo ""

# Use python to parse the selection JSON and iterate
python3 -c "
import json, subprocess, sys, os, time

selection_path = '$SELECTION_JSON'
incremental_path = '$INCREMENTAL_JSON'
runner_path = '$RUNNER'
results_dir = '$RESULTS_DIR'

with open(selection_path) as f:
    selection = json.load(f)

# Load existing incremental results
with open(incremental_path) as f:
    incremental = json.load(f)

controls = selection['selected_controls']
total = len(controls)
n_done = 0
n_skip = 0
n_fail = 0

for i, ctrl in enumerate(controls):
    ensembl = ctrl['ensembl']
    symbol = ctrl['symbol']

    # Skip if already in incremental results with success
    if symbol in incremental and incremental[symbol].get('status') in ('success', 'reused'):
        print(f'[{i+1}/{total}] {symbol}: ALREADY DONE (shift={incremental[symbol].get(\"mean_shift\", \"?\"):.6f})')
        n_skip += 1
        continue

    print(f'')
    print(f'[{i+1}/{total}] Running {symbol} ({ensembl})...')

    # Output file for this gene
    gene_json = os.path.join(results_dir, f'result_{symbol}.json')

    # Run in fresh subprocess
    t0 = time.time()
    try:
        result = subprocess.run(
            ['python3', runner_path, ensembl, symbol, gene_json],
            capture_output=False,
            timeout=600,  # 10 min max per gene
        )

        # Load result
        if os.path.exists(gene_json):
            with open(gene_json) as f:
                gene_result = json.load(f)
            gene_result['bin'] = ctrl['bin']
            gene_result['frequency'] = ctrl['frequency']
            gene_result['matched_to'] = ctrl['matched_to']
            incremental[symbol] = gene_result

            if gene_result.get('status') in ('success', 'reused'):
                n_done += 1
                print(f'  -> OK: shift={gene_result[\"mean_shift\"]:.6f}')
            else:
                n_fail += 1
                print(f'  -> FAILED: {gene_result.get(\"error\", \"unknown\")}')
        else:
            n_fail += 1
            incremental[symbol] = {
                'ensembl': ensembl,
                'symbol': symbol,
                'status': 'failed',
                'error': f'No output file produced (exit code {result.returncode})',
            }
            print(f'  -> No output file produced')
    except subprocess.TimeoutExpired:
        n_fail += 1
        incremental[symbol] = {
            'ensembl': ensembl,
            'symbol': symbol,
            'status': 'failed',
            'error': 'Timeout (600s)',
        }
        print(f'  -> TIMEOUT')
    except Exception as e:
        n_fail += 1
        incremental[symbol] = {
            'ensembl': ensembl,
            'symbol': symbol,
            'status': 'failed',
            'error': str(e),
        }
        print(f'  -> ERROR: {e}')

    # Save incremental results after every gene
    with open(incremental_path, 'w') as f:
        json.dump(incremental, f, indent=2)

    elapsed = time.time() - t0
    print(f'  Time: {elapsed:.1f}s | Done: {n_done+n_skip} | Failed: {n_fail} | Remaining: {total-i-1}')

print(f'')
print(f'=== COMPLETE ===')
print(f'Successful: {n_done} new + {n_skip} previously done = {n_done+n_skip}')
print(f'Failed: {n_fail}')
print(f'Results saved to: {incremental_path}')
"
