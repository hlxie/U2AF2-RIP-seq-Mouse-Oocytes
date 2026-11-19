#!/usr/bin/env python3
"""Reconstructed featureCounts invocation; paired-end, unstranded fragment mode."""
import argparse
import json
import subprocess
from pathlib import Path
import pandas as pd
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--samples', type=Path, required=True)
p.add_argument('--root', type=Path, default=Path('results'))
p.add_argument('--threads', type=int, default=8)
a = p.parse_args()
s = pd.read_csv(a.samples, sep='\t')
if len(s) != 11 or s.count_column.duplicated().any():
    raise ValueError('Require eleven independent BAMs, no merged controls')
for name in s.count_column:
    if not Path(name).is_file():
        raise FileNotFoundError(name)
out = a.root / 'RIP_peak_counts_ge2'
out.mkdir(parents=True, exist_ok=True)
cmd = ['featureCounts', '-T', str(a.threads), '-F', 'SAF', '-a', str(out/'peaks.saf'),
       '-o', str(out/'featureCounts.tsv'), '-p', '--countReadPairs', '-B', '-C',
       '--primary', '-s', '0', '-Q', '0', '-O', *s.count_column]
(out/'featureCounts_command.json').write_text(json.dumps(cmd, indent=2))
with (out/'featureCounts.log').open('w') as log:
    subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
