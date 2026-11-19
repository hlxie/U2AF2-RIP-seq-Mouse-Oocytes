#!/usr/bin/env python3
"""Reconstructed fragment-count normalization with both Input and IgG controls.

Denominator = sum of ALL featureCounts summary categories per library,
not sum of counts assigned to selected peaks. No statistical enrichment test.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts', type=Path, required=True, help='Raw featureCounts output with six annotation columns')
    parser.add_argument('--summary', type=Path, required=True, help='Matching featureCounts .summary')
    parser.add_argument('--samples', type=Path, required=True, help='TSV: sample_id, group, count_column (exact BAM header)')
    parser.add_argument('--root', type=Path, default=Path('results'))
    parser.add_argument('--pseudocount', type=float, required=True, help='Explicit CPM pseudocount; historical final value needs verification')
    args = parser.parse_args()
    if not np.isfinite(args.pseudocount) or args.pseudocount <= 0:
        parser.error('Pseudocount must be finite and positive')
    samples = pd.read_csv(args.samples, sep='\t')
    if samples.groupby('group').size().to_dict() != {'IP':5, 'Input':3, 'IgG':3}:
        raise ValueError('Require 5 IP, 3 Input and 3 IgG independent libraries')
    if samples.sample_id.duplicated().any() or samples.count_column.duplicated().any():
        raise ValueError('Duplicate sample or count column')
    raw = pd.read_csv(args.counts, sep='\t', comment='#').set_index('Geneid')
    summary = pd.read_csv(args.summary, sep='\t', index_col=0)
    count_dir = args.root/'RIP_peak_counts_ge2'
    coords = pd.read_csv(count_dir/'peak_coordinates.tsv', sep='\t', dtype={'chr':str})
    if raw.index.duplicated().any() or set(coords.peak_id) != set(raw.index):
        raise ValueError('Peak IDs do not match coordinates')
    raw = raw.loc[coords.peak_id]
    # Confirm identical intervals, not just matching generated peak IDs.
    if not (raw.Chr.astype(str).to_numpy() == coords.chr.to_numpy()).all() or not np.array_equal(raw.Start.to_numpy()-1, coords.start_0based.to_numpy()) or not np.array_equal(raw.End.to_numpy(), coords.end.to_numpy()):
        raise ValueError('featureCounts coordinates differ from consensus coordinates')
    values = raw[samples.count_column].apply(pd.to_numeric).to_numpy(float)
    totals = summary[samples.count_column].apply(pd.to_numeric).sum(axis=0).to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any() or not np.isfinite(totals).all() or (totals <= 0).any():
        raise ValueError('Invalid counts or library totals')
    counts = pd.DataFrame(values, columns=samples.sample_id)
    counts.insert(0, 'peak_id', coords.peak_id.to_numpy())
    counts.to_csv(count_dir/'peak_counts.tsv', sep='\t', index=False)
    samples.to_csv(count_dir/'samples.tsv', sep='\t', index=False)
    audit = samples.copy()
    audit['processed_fragments'] = totals
    audit['sum_peak_assignments'] = values.sum(axis=0)
    audit.to_csv(count_dir/'library_sizes.tsv', sep='\t', index=False)
    cpm = values / totals * 1e6
    out = args.root/'RIP_enrichment_ge2'
    out.mkdir(parents=True, exist_ok=True)
    table = coords.copy()
    for i, sample in enumerate(samples.sample_id):
        table[sample+'_CPM'] = cpm[:, i]
    table.to_csv(out/'peak_CPM.tsv', sep='\t', index=False)
    for group in ['Input', 'IgG']:
        result = table.copy()
        control = cpm[:, (samples.group == group).to_numpy()].mean(axis=1)
        result['control_mean_CPM'] = control
        result['pseudocount_CPM'] = args.pseudocount
        for i in np.where(samples.group.to_numpy() == 'IP')[0]:
            ratio = (cpm[:, i]+args.pseudocount)/(control+args.pseudocount)
            result[samples.sample_id.iloc[i]+'_log2_enrichment'] = np.log2(ratio)
        result.to_csv(out/f'enrichment_vs_{group}.tsv', sep='\t', index=False)
    (out/'settings.json').write_text(json.dumps({k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}, indent=2))
    print('Saved per-library CPM and descriptive IP/control log2 enrichment')


if __name__ == '__main__':
    main()
