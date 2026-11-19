#!/usr/bin/env python3
"""Reconstructed consensus builder, starting from FIVE final per-IP BED files.

Input BEDs must already represent the chosen control/background strategy.
This does not call peaks, intersect Input/IgG peak calls, or assign strand.
"""
import argparse
from collections import defaultdict
from pathlib import Path
import pandas as pd


def merge_intervals(intervals):
    result = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError(f"Invalid BED interval: {start}, {end}")
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def read_bed(path):
    regions = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            chrom, start, end = line.split()[:3]
            regions[chrom].append((int(start), int(end)))
    return {chrom: merge_intervals(v) for chrom, v in regions.items()}


def build_segments(beds):
    events = defaultdict(lambda: defaultdict(list))
    for i, bed in enumerate(beds):
        for chrom, intervals in bed.items():
            for start, end in intervals:
                events[chrom][start].append((i, 1))
                events[chrom][end].append((i, -1))
    segments = []
    for chrom in sorted(events):
        active = [0] * len(beds)
        previous = None
        for position in sorted(events[chrom]):
            if previous is not None and position > previous and any(active):
                flags = [int(n > 0) for n in active]
                members = ','.join(str(i+1) for i, flag in enumerate(flags) if flag)
                segments.append([chrom, previous, position, sum(flags), members, *flags])
            for i, change in events[chrom][position]:
                active[i] += change
            previous = position
    return segments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--peaks', type=Path, required=True, help='TSV: sample_id, bed (path relative to current directory)')
    parser.add_argument('--root', type=Path, default=Path('results'))
    args = parser.parse_args()
    manifest = pd.read_csv(args.peaks, sep='\t')
    if len(manifest) != 5 or manifest.sample_id.duplicated().any():
        raise ValueError('Exactly five unique IP replicates are required')
    segments = build_segments([read_bed(p) for p in manifest.bed])
    selected = [r for r in segments if r[3] >= 2]
    if not selected:
        raise ValueError('No intervals supported by at least two IPs')
    qc = args.root / 'consensus_qc'
    counts = args.root / 'RIP_peak_counts_ge2'
    qc.mkdir(parents=True, exist_ok=True)
    counts.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(segments).to_csv(qc/'all.segments.tsv', sep='\t', header=False, index=False)
    pd.DataFrame(selected).to_csv(qc/'support_ge2.segments.tsv', sep='\t', header=False, index=False)
    (qc/'ip_order.txt').write_text('\n'.join(manifest.sample_id) + '\n')
    manifest.to_csv(qc/'input_peaks.tsv', sep='\t', index=False)
    grouped = defaultdict(list)
    for chrom, start, end, *_ in selected:
        grouped[chrom].append((start, end))
    rows = []
    for chrom in sorted(grouped):
        for start, end in merge_intervals(grouped[chrom]):
            rows.append([f'peak_{len(rows)+1:06d}', chrom, start, end, end-start])
    coords = pd.DataFrame(rows, columns=['peak_id', 'chr', 'start_0based', 'end', 'length'])
    coords.to_csv(counts/'peak_coordinates.tsv', sep='\t', index=False)
    coords[['chr','start_0based','end','peak_id']].to_csv(qc/'consensus_ge2.bed', sep='\t', index=False, header=False)
    saf = pd.DataFrame({'GeneID':coords.peak_id, 'Chr':coords.chr, 'Start':coords.start_0based+1, 'End':coords.end, 'Strand':'.'})
    saf.to_csv(counts/'peaks.saf', sep='\t', index=False)
    print(f'{len(coords)} merged consensus peaks; support segments retain original IP identities')


if __name__ == '__main__':
    main()
