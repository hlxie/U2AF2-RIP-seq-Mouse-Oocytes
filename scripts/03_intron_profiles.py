#!/usr/bin/env python3
"""Reconstructed, strand-oriented intronic RIP profiles from existing CPM bigWigs.

Analyze consensus-peak-overlapping introns; subtract ALL annotated exon bases.
These are conditional, peak-masked profiles, not an unbiased binding-site assay.
"""
import argparse
from bisect import bisect_right
from collections import defaultdict
from contextlib import ExitStack
import json
from pathlib import Path
import re
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from importlib.util import spec_from_file_location, module_from_spec

spec = spec_from_file_location('consensus', Path(__file__).with_name('01_consensus.py'))
intervals = module_from_spec(spec)
spec.loader.exec_module(intervals)


def annotation(path):
    transcripts = defaultdict(list)
    exons = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if len(fields) != 9 or fields[2] != 'exon':
                continue
            attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"', fields[8]))
            chrom, strand = fields[0], fields[6]
            if strand not in ['+', '-'] or 'transcript_id' not in attrs:
                raise ValueError('Exons need transcript_id and known strand')
            start, end = int(fields[3])-1, int(fields[4])
            transcripts[chrom, strand, attrs['transcript_id']].append((start,end))
            exons[chrom].append((start,end))
    introns = set()
    for (chrom, strand, tx), pieces in transcripts.items():
        pieces = intervals.merge_intervals(pieces)
        for left, right in zip(pieces, pieces[1:]):
            if left[1] < right[0]:
                introns.add((chrom,left[1],right[0],strand))
    return sorted(introns), {ch:intervals.merge_intervals(v) for ch,v in exons.items()}


def index_bed(regions):
    return {ch:(v,[b for a,b in v]) for ch,v in regions.items()}


def mask_for(index, chrom, start, end):
    mask = np.zeros(end-start, dtype=bool)
    regions, ends = index.get(chrom, ([],[]))
    i = bisect_right(ends,start)
    for left,right in regions[i:]:
        if left >= end:
            break
        mask[max(left,start)-start:min(right,end)-start] = True
    return mask


def orient(values, strand):
    return values if strand == '+' else values[::-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gtf',type=Path,required=True)
    parser.add_argument('--peaks',type=Path,required=True,help='Consensus BED, exact same genome and chromosome names')
    parser.add_argument('--samples',type=Path,required=True,help='TSV with sample_id, group, bigwig (unmasked CPM)')
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--min-intron-length',type=int,default=100)
    parser.add_argument('--bins',type=int,default=100)
    parser.add_argument('--upstream-bp',type=int,default=100)
    args = parser.parse_args()
    if min(args.min_intron_length,args.bins,args.upstream_bp) <= 0:
        parser.error('Lengths and bins must be positive')
    samples = pd.read_csv(args.samples,sep='\t')
    if samples.groupby('group').size().to_dict() != {'IP':5,'Input':3,'IgG':3} or samples.sample_id.duplicated().any():
        raise ValueError('Expected 5 IP, 3 Input and 3 IgG unique libraries')
    introns, exons = annotation(args.gtf)
    exon_index = index_bed(exons)
    peak_regions = intervals.read_bed(args.peaks)
    if not set(peak_regions) & set(exons):
        raise ValueError('No common chromosomes between GTF and peaks')
    peak_index = index_bed(peak_regions)
    # The SAME selected introns are used for every sample and both profiles.
    cutoff = max(args.min_intron_length,args.bins,args.upstream_bp)
    selected = []
    for chrom,start,end,strand in introns:
        if end-start < cutoff:
            continue
        mask = mask_for(peak_index,chrom,start,end) & ~mask_for(exon_index,chrom,start,end)
        if mask.any():
            selected.append((chrom,start,end,strand))
    if not selected:
        raise ValueError('No eligible peak-overlapping introns after exon exclusion')
    args.out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(selected,columns=['chr','start_0based','end','strand']).to_csv(args.out/'selected_introns.tsv',sep='\t',index=False)
    acceptors = sorted({(ch, end if st=='+' else start, st) for ch,start,end,st in selected})
    pd.DataFrame(acceptors,columns=['chr','acceptor_boundary_0based','strand']).to_csv(args.out/'unique_acceptors.tsv',sep='\t',index=False)
    sample_rows = []
    with ExitStack() as stack:
        handles = []
        for path in samples.bigwig:
            handle = pyBigWig.open(str(path))
            stack.callback(handle.close)
            handles.append(handle)
        scaled_sum = np.zeros((len(samples),args.bins))
        for chrom,start,end,strand in selected:
            mask = mask_for(peak_index,chrom,start,end) & ~mask_for(exon_index,chrom,start,end)
            for i,handle in enumerate(handles):
                if chrom not in handle.chroms() or end > handle.chroms(chrom):
                    raise ValueError(f'bigWig/GTF coordinate mismatch: {chrom}:{start}-{end}')
                values = np.asarray(handle.values(chrom,start,end),dtype=float)
                if len(values) != end-start:
                    raise ValueError('Empty or invalid bigWig')
                values = np.nan_to_num(values,nan=0.0) * mask
                values = orient(values,strand)
                scaled_sum[i] += np.array([part.mean() for part in np.array_split(values,args.bins)])
        upstream_sum = np.zeros((len(samples),args.upstream_bp))
        # Deduplicate acceptors: transcripts sharing an acceptor get one vote here.
        for chrom,boundary,strand in acceptors:
            start,end = (boundary-args.upstream_bp,boundary) if strand=='+' else (boundary,boundary+args.upstream_bp)
            mask = mask_for(peak_index,chrom,start,end) & ~mask_for(exon_index,chrom,start,end)
            for i,handle in enumerate(handles):
                values = np.nan_to_num(np.asarray(handle.values(chrom,start,end),dtype=float),nan=0.0)*mask
                upstream_sum[i] += orient(values,strand)
        for kind,matrix,x,n in [
            ('scaled_intron',scaled_sum/len(selected),(np.arange(args.bins)+.5)/args.bins,len(selected)),
            ('acceptor_upstream',upstream_sum/len(acceptors),np.arange(-args.upstream_bp,0),len(acceptors))]:
            for i,sample in samples.iterrows():
                for position,value in zip(x,matrix[i]):
                    sample_rows.append([kind,sample.sample_id,sample.group,position,value,n])
    profile = pd.DataFrame(sample_rows,columns=['profile','sample_id','group','position','mean_CPM_coverage','n_features'])
    profile.to_csv(args.out/'sample_profiles.tsv',sep='\t',index=False)
    stats = profile.groupby(['profile','group','position']).mean_CPM_coverage.agg(['mean','std','count']).reset_index()
    stats.to_csv(args.out/'group_profiles.tsv',sep='\t',index=False)
    colors = {'IP':'#376DB5','Input':'#E99B28','IgG':'#777777'}
    for kind in profile.profile.unique():
        fig,ax = plt.subplots(figsize=(6,4))
        for group in ['IP','Input','IgG']:
            d = profile[(profile.profile==kind)&(profile.group==group)]
            for _,r in d.groupby('sample_id'):
                ax.plot(r.position,r.mean_CPM_coverage,color=colors[group],lw=.6,alpha=.3)
            g = stats[(stats.profile==kind)&(stats.group==group)]
            ax.plot(g.position,g['mean'],color=colors[group],label=group)
        ax.set_ylabel('Peak-masked mean CPM coverage')
        ax.set_xlabel('Intron position (5′ donor → 3′ acceptor)' if kind=='scaled_intron' else 'Position upstream of 3′ acceptor (nt; last intronic base = −1)')
        ax.legend(frameon=False)
        fig.tight_layout()
        for extension in ['pdf','png']:
            fig.savefig(args.out/f'{kind}.{extension}',dpi=220)
        plt.close(fig)
    settings = {k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}
    settings.update(effective_min_intron_length=cutoff,n_introns=len(selected),n_unique_acceptors=len(acceptors),missing_bigwig_bases='zero',outside_peak_or_inside_any_exon='zero')
    (args.out/'settings.json').write_text(json.dumps(settings,indent=2))
    print(f'{len(selected)} introns; {len(acceptors)} unique acceptors')


if __name__ == '__main__':
    main()
