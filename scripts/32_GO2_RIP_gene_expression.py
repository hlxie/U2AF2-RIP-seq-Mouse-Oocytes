#!/usr/bin/env python3
# Python >=3.7. Run on the HPC, from the directory containing the inputs.
import argparse, csv, json, re, subprocess, sys
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

REF='reference/Mus_musculus.GRCm38.101.gtf'
def sid(s):
    s=str(s).strip()
    return re.sub(r'\.\d+$','',s) if s.startswith(('ENSMUSG','ENSMUST')) else s

def table(path):
    with open(str(path),'rb') as f: magic=f.read(8)
    if magic.startswith(b'PK\x03\x04') or magic.startswith(b'\xd0\xcf\x11\xe0'):
        return pd.read_excel(path,dtype=str).fillna('')
    for enc in ['utf-8-sig','gb18030']:
        try:
            d=pd.read_csv(path,sep='\t',dtype=str,encoding=enc).fillna('')
            if d.shape[1]<2: raise ValueError('Expected tab-delimited table: '+str(path))
            return d
        except UnicodeDecodeError: pass
    raise ValueError('Cannot decode '+str(path))

def attrs(s):
    return {m.group(1):m.group(2).strip().strip('"') for m in re.finditer(r'(\w+)\s+([^;]+)',s)}

def gtf(path):
    with open(path) as f:
        for line in f:
            if line.startswith('#'): continue
            p=line.rstrip('\n').split('\t')
            if len(p)!=9: raise ValueError('GTF must have 9 tab-separated fields: '+str(path))
            a=attrs(p[8]); yield p,a

def write(d,path,index=False): d.to_csv(path,sep='\t',index=index)

def prepare(a,o):
    genes={}; rt=defaultdict(set)
    for p,z in gtf(a.reference):
        g=sid(z.get('gene_id','')); t=sid(z.get('transcript_id',''))
        if not g: continue
        if p[2]=='gene':
            genes[g]=dict(gene_id=g,gene_name=z.get('gene_name',g),chrom=p[0].removeprefix('chr') if sys.version_info>=(3,9) else re.sub('^chr','',p[0]),start=int(p[3])-1,end=int(p[4]),strand=p[6])
        if t: rt[t].add(g)
    if not genes: raise ValueError('Reference must include gene features.')
    mapping=defaultdict(set)
    # The full known annotation, NOT SUPPA-filtered annotation.
    for p,z in gtf(a.known):
        t=sid(z.get('transcript_id','')); g=sid(z.get('gene_id',''))
        if t:
            if g in genes: mapping[t].add(g)
            elif len(rt.get(t,set()))==1: mapping[t].update(rt[t])
    novel=set()
    for p,z in gtf(a.novel):
        t=sid(z.get('transcript_id',''))
        if t: novel.add(t)
    cl=table(a.classification); cl.columns=cl.columns.str.strip()
    required={'isoform','associated_gene','structural_category'}
    if not required.issubset(cl.columns): raise ValueError('SQANTI needs columns '+str(required))
    if cl.isoform.duplicated().any(): raise ValueError('Duplicate isoform IDs in classification.')
    for _,r in cl.iterrows():
        t=sid(r.isoform)
        gs=[sid(v) for v in re.split(r'[,;|\s]+',r.associated_gene.strip()) if v]
        if t in novel and r.structural_category in ['full-splice_match','incomplete-splice_match','novel_in_catalog','novel_not_in_catalog'] and len(gs)==1 and gs[0] in genes:
            mapping[t].add(gs[0])
    txmap={t:next(iter(gs)) for t,gs in mapping.items() if len(gs)==1}
    write(pd.DataFrame([(t,','.join(sorted(gs)),len(gs)==1) for t,gs in mapping.items()],columns=['transcript_id','gene_ids','unique']),o/'transcript_gene_mapping.tsv')
    ann=pd.DataFrame(genes.values()).set_index('gene_id'); ann['gene_length']=ann.end-ann.start
    write(ann,o/'gene_annotation.tsv',True)
    metas=[]
    for group,path in [('WT',a.wt),('CKO',a.cko)]:
        d=table(path); d.columns=d.columns.str.strip()
        for c in ['batch','Name','Cell_Line','mouse_id','Pass_QC']:
            if c not in d: raise ValueError('Missing metadata column '+c)
        for c in d.columns: d[c]=d[c].astype(str).str.strip()
        d=d[(d.Cell_Line==a.stage)&(pd.to_numeric(d.Pass_QC,errors='coerce')==1)].copy()
        if (d.mouse_id=='').any(): raise ValueError('Empty mouse_id')
        d['cell_id']=d.batch+'_'+d.Name; d['genotype']=group; d['sample_id']=group+'__'+d.mouse_id
        metas.append(d)
    meta=pd.concat(metas,ignore_index=True)
    if meta.cell_id.duplicated().any(): raise ValueError('Duplicate cell_id across metadata; cannot safely assign.')
    with open(a.matrix,encoding='utf-8-sig',newline='') as f:
        reader=csv.reader(f,delimiter='\t'); h=next(reader)
        # Header may contain two annotation names, or only the cell names.
        first=next(reader); n=len(first)
        if len(h)==n: cells=h[2:]
        elif len(h)==n-2: cells=h
        else: raise ValueError('Matrix header/data width mismatch: %s/%s'%(len(h),n))
        if len(cells)!=len(set(cells)): raise ValueError('Duplicate matrix cell columns')
        write(meta.assign(matched_matrix=meta.cell_id.isin(cells)),o/'metadata_match_audit.tsv')
        meta=meta[meta.cell_id.isin(cells)].copy()
        nc=meta.groupby('sample_id').size(); keep=nc[nc>=a.min_cells_per_mouse].index
        meta=meta[meta.sample_id.isin(keep)].copy()
        write(meta,o/'selected_cells.tsv')
        samples=meta.groupby(['sample_id','genotype','mouse_id']).size().reset_index(name='n_cells').set_index('sample_id')
        if any((samples.genotype==g).sum()<2 for g in ['WT','CKO']): raise ValueError('Need >=2 mice in each genotype')
        write(samples,o/'mouse_samples.tsv',True)
        sample_names=list(samples.index); si={s:i for i,s in enumerate(sample_names)}
        ci={c:i+2 for i,c in enumerate(cells)}
        loc=[ci[c] for c in meta.cell_id]; targets=np.array([si[s] for s in meta.sample_id])
        sums=defaultdict(lambda:np.zeros(len(samples))); total=np.zeros(len(samples)); mapped=np.zeros(len(samples)); audit=[]; seen=set()
        def consume(row):
            if len(row)!=n: raise ValueError('Inconsistent matrix row width')
            t=sid(row[0])
            if t in seen: raise ValueError('Duplicate normalized transcript ID: '+t)
            seen.add(t)
            v=np.array([float(row[i]) for i in loc])
            if not np.isfinite(v).all() or (v<0).any(): raise ValueError('Invalid counts: '+t)
            b=np.bincount(targets,weights=v,minlength=len(samples)); total[:]+=b
            g=txmap.get(t)
            if g: sums[g]+=b; mapped[:]+=b
            else: audit.append((t,row[1],float(b.sum())))
        consume(first)
        for row in reader: consume(row)
    mat=pd.DataFrame.from_dict(sums,orient='index',columns=sample_names); mat.index.name='gene_id'
    if mat.empty or (mat.sum()==0).any(): raise ValueError('Empty mapped count matrix/sample')
    write(mat,o/'gene_mouse_raw_counts.tsv',True)
    write(pd.DataFrame(audit,columns=['transcript_id','matrix_annotation','GO2_count_sum']),o/'unmapped_transcripts.tsv')
    write(pd.DataFrame({'sample_id':sample_names,'all_transcript_count':total,'mapped_gene_count':mapped,'mapped_fraction':mapped/np.maximum(total,1)}),o/'count_mapping_coverage.tsv')
    print('GO2 cells/mice:',meta.groupby('genotype').agg({'cell_id':'size','sample_id':'nunique'}).to_string(),flush=True)
    print('Mapped genes:',len(mat),'Unmapped transcript rows:',len(audit),flush=True)
    peaks(a,o,ann)

def peaks(a,o,ann):
    p=table(a.peaks); p=p.rename(columns={'chr':'chrom'}); p['chrom']=p.chrom.str.replace(r'^chr','',regex=True)
    for c in ['start_0based','end']: p[c]=pd.to_numeric(p[c],errors='raise').astype(int)
    if p.peak_id.duplicated().any() or (p.end<=p.start_0based).any(): raise ValueError('Invalid/duplicate peak intervals')
    cols=[c for c in p if c.endswith('_CPM') and c.startswith('RIP-U2-')]
    if len(cols)!=5: raise ValueError('Expected 5 RIP-U2 IP CPM columns, got '+str(cols))
    vals=p[cols].apply(pd.to_numeric,errors='raise').values
    if not np.isfinite(vals).all() or (vals<0).any(): raise ValueError('Invalid IP CPM')
    # Fail rather than double count overlapping consensus intervals.
    for ch,d in p.groupby('chrom'):
        d=d.sort_values('start_0based')
        if len(d)>1 and (d.start_0based.values[1:]<np.maximum.accumulate(d.end.values[:-1])).any():
            raise ValueError('Overlapping consensus peaks on '+ch+': merge and recount first; cannot safely sum CPM.')
    hits=[]; assigned=[]
    for ch,d in p.groupby('chrom'):
        g=ann[ann.chrom==ch].sort_values('start'); starts=g.start.values; ends=g.end.values; ids=g.index.values
        for idx,r in d.iterrows():
            k=np.searchsorted(starts,r.end,side='left'); gs=ids[:k][ends[:k]>r.start_0based]
            hits.append((r.peak_id,','.join(gs),len(gs)))
            if len(gs)==1: assigned.append((idx,gs[0]))
    write(pd.DataFrame(hits,columns=['peak_id','overlap_gene_ids','n_overlapping_genes']),o/'peak_gene_assignment.tsv')
    metrics=[]
    for g,pairs in pd.DataFrame(assigned,columns=['idx','gene_id']).groupby('gene_id'):
        inds=pairs.idx.values.astype(int); c=vals[inds]; lengths=(p.loc[inds,'end']-p.loc[inds,'start_0based']).values/1000.; den=c/lengths[:,None]
        r={'gene_id':g,'n_peaks':len(inds),'peak_bp':lengths.sum()*1000}
        for agg,fn in [('median',np.median),('mean',np.mean)]:
            v=fn(den,axis=1)
            for name,fun in [('median_density',np.median),('max_density',np.max),('mean_density',np.mean),('sum_density',np.sum)]:
                r[agg+'IP__'+name]=float(fun(v))
            r[agg+'IP__total_CPM']=float(fn(c.sum(axis=0)))
            r[agg+'IP__weighted_density']=float(fn(c.sum(axis=0)/lengths.sum()))
        metrics.append(r)
    if not metrics: raise ValueError('No unambiguous peak-associated genes; check chromosome naming.')
    write(pd.DataFrame(metrics).set_index('gene_id'),o/'gene_RIP_metrics.tsv',True)
    write(pd.DataFrame({'IP_column':cols}),o/'IP_columns.tsv')

def python_expression(a,o):
    x=pd.read_csv(o/'gene_mouse_raw_counts.tsv',sep='\t',index_col=0)
    samples=pd.read_csv(o/'mouse_samples.tsv',sep='\t',index_col=0).loc[x.columns]
    cov=pd.read_csv(o/'count_mapping_coverage.tsv',sep='\t').set_index('sample_id').loc[x.columns]
    # Full transcript UMI totals, including unassigned transcripts, avoid renormalizing
    # to a potentially different mapped subset in each mouse.
    lib=cov.all_transcript_count
    if (lib<=0).any(): raise ValueError('Zero mouse library total')
    cpm=x.div(lib,axis=1)*1e6
    masks={g:(samples.genotype==g).values for g in ['WT','CKO']}
    if any(m.sum()<a.min_mice for m in masks.values()): raise ValueError('Insufficient mice')
    supports={g:(cpm.loc[:,m]>=a.min_cpm).sum(axis=1) for g,m in masks.items()}
    keep=(supports['WT']>=a.min_mice)|(supports['CKO']>=a.min_mice)
    write(pd.DataFrame({'eligible':keep,'n_WT_CPM_pass':supports['WT'],'n_CKO_CPM_pass':supports['CKO']}),o/'gene_expression_filter.tsv',True)
    means={g:cpm.loc[:,m].mean(axis=1) for g,m in masks.items()}
    d=pd.DataFrame({'WT_mean_CPM':means['WT'],'CKO_mean_CPM':means['CKO']})
    for pc in sorted(set([.1,.5,1.,a.pseudocount])):
        d['log2FC_pc'+str(pc)]=np.log2((means['WT']+pc)/(means['CKO']+pc))
    d['log2_WT_over_CKO']=np.log2((means['WT']+a.pseudocount)/(means['CKO']+a.pseudocount))
    d['method']='equal_mouse_mean_full_library_CPM_descriptive'
    d.index.name='gene_id'
    # Filename retained for downstream interface; no DE P values are invented.
    write(d.loc[keep],o/'gene_DE_WT_vs_CKO.tsv',True)
    write(d,o/'gene_expression_changes_all.tsv',True)
    write(cpm,o/'gene_mouse_CPM.tsv',True)
    write(pd.DataFrame({'full_library_total':lib}),o/'normalization_factors.tsv',True)
    (o/'EXPRESSION_METHOD.txt').write_text('Python descriptive expression fold change; NOT edgeR differential testing.\nEach mouse is normalized by full transcript UMI total then mice receive equal weight.\nPseudocount CPM = '+str(a.pseudocount)+'\nNo gene DE P values or FDR. See correlation_summary for cross-gene association statistics.\n')
    print('Expression eligible:',int(keep.sum()),'; Python descriptive fold change, not formal DE.',flush=True)

def main():
    ap=argparse.ArgumentParser(description='GO2 gene pseudobulk + WT RIP association. Python3.7+')
    ap.add_argument('--matrix',default='rc_trans_merge.txt'); ap.add_argument('--known',default='gtf_known.gtf'); ap.add_argument('--novel',default='gtf_novel.gtf')
    ap.add_argument('--classification',default='SQANTI3_novel_reannotation_fixed/novel_reannot_classification.txt')
    ap.add_argument('--reference',default=REF); ap.add_argument('--wt',default='WT.sample.info.mouse_id.xls'); ap.add_argument('--cko',default='CKO.sample.info.mouse_id.xls')
    ap.add_argument('--peaks',default='RIP_enrichment_ge2/enrichment_vs_IgG.tsv'); ap.add_argument('--stage',default='GO2')
    ap.add_argument('--min-cells-per-mouse',type=int,default=1); ap.add_argument('--min-cpm',type=float,default=1.); ap.add_argument('--min-mice',type=int,default=2)
    ap.add_argument('--method',choices=['python','edgeR'],default='python'); ap.add_argument('--pseudocount',type=float,default=0.5); ap.add_argument('--out',default='GO2_RIP_gene_expression_v2'); ap.add_argument('--rscript',default='Rscript'); ap.add_argument('--prepare-only',action='store_true')
    a=ap.parse_args()
    if a.min_mice<2 or a.min_cpm<=0 or a.pseudocount<=0 or a.min_cells_per_mouse<1: ap.error('Invalid threshold')
    if a.method=='edgeR' and not a.prepare_only:
        try:
            subprocess.run([a.rscript,'-e',"if (!requireNamespace('edgeR', quietly=TRUE)) stop('Missing edgeR; install it in this R environment or use --method python')"],check=True)
        except (OSError,subprocess.CalledProcessError):
            sys.exit('edgeR preflight failed BEFORE reading the matrix. Use --method python, or set --rscript to an R environment containing edgeR.')
    o=Path(a.out).resolve(); o.mkdir(parents=True,exist_ok=True)
    (o/'settings.json').write_text(json.dumps(vars(a),indent=2))
    prepare(a,o)
    if a.prepare_only: return
    here=Path(__file__).resolve().parent
    if a.method=='edgeR':
        subprocess.run([a.rscript,str(here/'33_GO2_gene_DE_edgeR.R'),str(o),str(a.min_cpm),str(a.min_mice)],check=True)
    else:
        python_expression(a,o)
    subprocess.run([sys.executable,str(here/'34_RIP_gene_correlations.py'),'--out',str(o)],check=True)
    print('Completed:',o)
if __name__=='__main__': main()
