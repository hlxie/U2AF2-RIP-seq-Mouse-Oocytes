#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def bh(p):
    p=np.asarray(p,float); q=np.full(len(p),np.nan); ix=np.where(np.isfinite(p))[0]
    order=ix[np.argsort(p[ix])]; n=len(order)
    if n: q[order]=np.minimum(1,np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1])
    return q

def corr(x,y):
    if len(x)<5 or np.ptp(x)==0 or np.ptp(y)==0: return np.nan,np.nan
    r,p=stats.spearmanr(x,y); return float(r),float(p)

def partial(x,y,z):
    # Rank residualization: descriptive partial Spearman, approximate t-test.
    x=stats.rankdata(x); y=stats.rankdata(y)
    z=np.column_stack([np.ones(len(x))]+[stats.rankdata(z[:,j]) for j in range(z.shape[1])])
    rank=np.linalg.matrix_rank(z); df=len(x)-rank-1
    rx=x-z.dot(np.linalg.lstsq(z,x,rcond=None)[0]); ry=y-z.dot(np.linalg.lstsq(z,y,rcond=None)[0])
    if df<=0 or np.std(rx)<1e-9 or np.std(ry)<1e-9: return np.nan,np.nan,rx,ry
    r=float(np.corrcoef(rx,ry)[0,1]); p=float(2*stats.t.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-15)),df))
    return r,p,rx,ry

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',default='GO2_RIP_gene_expression'); a=ap.parse_args(); o=Path(a.out)
    def read(n): return pd.read_csv(o/n,sep='\t').set_index('gene_id')
    rip=read('gene_RIP_metrics.tsv'); de=read('gene_DE_WT_vs_CKO.tsv'); ann=read('gene_annotation.tsv')
    d=rip.join(de,how='inner').join(ann,how='inner')
    metrics=[c for c in rip if '__' in c]
    cols=metrics+['log2_WT_over_CKO','WT_mean_CPM','gene_length','n_peaks']
    d=d[np.isfinite(d[cols].values).all(axis=1)].copy()
    if len(d)<10: raise ValueError('Fewer than 10 shared finite expressed peak-associated genes')
    d.to_csv(o/'gene_RIP_expression_merged.tsv',sep='\t')
    rows=[]; qc=[]; groups=[]
    fig,axes=plt.subplots(4,3,figsize=(13,14)); y=d.log2_WT_over_CKO.values
    z=np.column_stack([np.log1p(d.WT_mean_CPM),np.log(d.gene_length),d.n_peaks])
    for ax,m in zip(axes.flat,metrics):
        x=np.log2(1+d[m].values); r,p=corr(x,y); pr,pp,rx,ry=partial(x,y,z)
        rows.append(dict(metric=m,n_genes=len(d),spearman_rho=r,p_value=p,partial_rank_r=pr,partial_p_approx=pp))
        for c in ['WT_mean_CPM','gene_length','n_peaks']:
            qr,qp=corr(x,d[c]); qc.append(dict(metric=m,covariate=c,rho=qr,p=qp))
        ax.scatter(x,y,s=9,alpha=.35,color='#376AA5',rasterized=True); ax.axhline(0,color='grey',lw=.7)
        ax.set(title=m+'\nrho=%.3f; n=%d'%(r,len(d)),xlabel='log2(1 + RIP metric)',ylabel='log2(WT / CKO)')
        # Same genes for all plots; no significance/downregulation selection.
        f,b=plt.subplots(figsize=(5,4.5)); b.scatter(x,y,s=10,alpha=.4,color='#376AA5',rasterized=True)
        b.axhline(0,color='grey',lw=.8); b.set(xlabel='log2(1 + '+m+')',ylabel='log2(WT / CKO)',title='Spearman rho=%.3f; P=%.3g; n=%d'%(r,p,len(d)))
        f.tight_layout()
        for ext in ['pdf','svg','png']: f.savefig(o/(m+'_scatter.'+ext),dpi=250)
        plt.close(f)
        f,b=plt.subplots(figsize=(5,4.5)); b.scatter(rx,ry,s=10,alpha=.35,color='#376AA5',rasterized=True)
        b.set(xlabel='RIP rank residual',ylabel='Expression-change rank residual',title='Adjusted rank r=%.3f (descriptive)'%pr); f.tight_layout(); f.savefig(o/(m+'_adjusted.pdf')); plt.close(f)
        # Quantile ties stay together; group sizes can differ.
        bins=pd.qcut(d[m],3,labels=False,duplicates='drop')
        for k in sorted(bins.dropna().unique()):
            v=y[bins.values==k]
            groups.append(dict(metric=m,signal_quantile=int(k)+1,n_genes=len(v),median_log2_WT_over_CKO=np.median(v),q25=np.quantile(v,.25),q75=np.quantile(v,.75)))
    fig.tight_layout(); fig.savefig(o/'all_12_RIP_correlations.pdf'); fig.savefig(o/'all_12_RIP_correlations.png',dpi=220); plt.close(fig)
    res=pd.DataFrame(rows); res['BH_FDR_12_metrics']=bh(res.p_value); res['partial_BH_FDR_12_metrics']=bh(res.partial_p_approx)
    res.to_csv(o/'correlation_summary.tsv',sep='\t',index=False)
    pd.DataFrame(qc).to_csv(o/'RIP_covariate_correlations.tsv',sep='\t',index=False)
    pd.DataFrame(groups).to_csv(o/'RIP_signal_quantile_summary.tsv',sep='\t',index=False)
    sensitivity=[]
    for fc in [c for c in d if c.startswith('log2FC_pc')]:
        for m in metrics:
            r,p=corr(d[m].values,d[fc].values)
            sensitivity.append(dict(metric=m,fold_change=fc,n_genes=len(d),rho=r,p_value=p))
    if sensitivity:
        sen=pd.DataFrame(sensitivity); sen['BH_FDR']=bh(sen.p_value)
        sen.to_csv(o/'pseudocount_sensitivity.tsv',sep='\t',index=False)
    print(res.to_string(index=False))
if __name__=='__main__': main()
