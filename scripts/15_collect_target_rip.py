#!/usr/bin/env python3
"""Collect target-gene RIP evidence; run from the existing RIP analysis directory."""
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
import re

import pandas as pd

GTF = Path("reference/Mus_musculus.GRCm38.101.gtf")
TARGETS = ["Zp2", "Bmp15", "Rspo2", "Ezh2", "Atrx", "Mapk1", "Limk2"]

# This order MUST match the five original multiinter inputs.
REP_SAMPLES = ["RIP-U2-5", "RIP-U2-7", "RIP-U2-8", "RIP-U2-9", "RIP-U2-10"]


def read_peak_table(path):
    table = pd.read_csv(path, sep="\t")
    if "peak_id" not in table.columns:
        raise ValueError(f"{path} 缺少 peak_id 列；实际列名：{table.columns.tolist()}")
    if table["peak_id"].isna().any() or table["peak_id"].duplicated().any():
        raise ValueError(f"{path} 存在空白或重复的 peak_id")
    return table


def main(gtf_path=GTF, root=Path(".")):
    root = Path(root)
    count_dir = root / "RIP_peak_counts_ge2"
    enrich_dir = root / "RIP_enrichment_ge2"
    segment_path = root / "consensus_qc/support_ge2.segments.tsv"
    out = root / "RIP_target_genes"
    out.mkdir(parents=True, exist_ok=True)

    # 1. Read only per-peak tables. enrichment_summary.tsv is NOT a per-peak table.
    tables = {
        "raw_count": read_peak_table(count_dir / "peak_counts.tsv"),
        "peak_CPM": read_peak_table(enrich_dir / "peak_CPM.tsv"),
        "enrichment_vs_Input": read_peak_table(enrich_dir / "enrichment_vs_Input.tsv"),
        "enrichment_vs_IgG": read_peak_table(enrich_dir / "enrichment_vs_IgG.tsv"),
    }
    coords = pd.read_csv(
        count_dir / "peak_coordinates.tsv", sep="\t", dtype={"chr": str}
    )
    needed = {"peak_id", "chr", "start_0based", "end", "length"}
    if not needed.issubset(coords.columns) or coords["peak_id"].duplicated().any():
        raise ValueError("peak_coordinates.tsv 列名异常或 peak_id 重复")

    # Use the same main-analysis peak set as the existing CPM table.
    main_ids = set(tables["peak_CPM"]["peak_id"])
    if not main_ids.issubset(set(coords["peak_id"])):
        raise ValueError("部分 CPM peak 在坐标表中不存在")
    coords = coords.loc[coords["peak_id"].isin(main_ids)].copy()

    # 2. Obtain target gene-body coordinates from the SAME mm10 annotation.
    genes = []
    with Path(gtf_path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = dict(re.findall(r'(\S+)\s+"([^"]*)"', fields[8]))
            name = attributes.get("gene_name", "")
            if name in TARGETS:
                genes.append({
                    "gene_name": name,
                    "gene_id": attributes.get("gene_id", ""),
                    "chr": fields[0],
                    "gene_start_0based": int(fields[3]) - 1,
                    "gene_end": int(fields[4]),
                    "strand": fields[6],
                })

    missing = set(TARGETS) - {gene["gene_name"] for gene in genes}
    if missing:
        raise ValueError(f"GTF 中未找到目标基因：{sorted(missing)}")
    pd.DataFrame(genes).to_csv(out / "target_gene_coordinates.tsv", sep="\t", index=False)

    records, overview = [], []
    for gene in genes:
        subset = coords.loc[
            (coords["chr"] == gene["chr"])
            & (coords["start_0based"] < gene["gene_end"])
            & (coords["end"] > gene["gene_start_0based"])
        ]
        overview.append({**gene, "n_overlapping_consensus_peaks": len(subset)})
        for _, peak in subset.iterrows():
            left, right = int(peak["start_0based"]), int(peak["end"])
            records.append({
                **gene,
                "peak_id": peak["peak_id"],
                "peak_start_0based": left,
                "peak_end": right,
                "peak_length": right - left,
                "gene_overlap_bp": min(right, gene["gene_end"])
                    - max(left, gene["gene_start_0based"]),
            })

    overview = pd.DataFrame(overview)
    overview.to_csv(out / "target_gene_overview.tsv", sep="\t", index=False)
    if not records:
        pd.DataFrame(columns=["gene_name", "peak_id"]).to_csv(
            out / "target_peak_evidence.tsv", sep="\t", index=False
        )
        print("目标基因体内未检出主分析集合中的 consensus peak。")
        print(overview.to_string(index=False))
        return
    result = pd.DataFrame(records)

    # 3. Read the original >=2/5 multiinter segments, preserving per-IP support.
    segments = defaultdict(list)
    with segment_path.open() as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 10:
                raise ValueError(f"{segment_path}:{lineno} 预期为10列 multiinter 结果")
            left, right, support = map(int, fields[1:4])
            flags = tuple(map(int, fields[5:10]))
            if not (left < right and support >= 2 and all(x in (0, 1) for x in flags)
                    and sum(flags) == support):
                raise ValueError(f"multiinter 第 {lineno} 行坐标或支持信息异常")
            segments[fields[0]].append((left, right, support, flags))

    indices = {}
    for chrom, entries in segments.items():
        entries.sort()
        if any(entries[i][0] < entries[i-1][1] for i in range(1, len(entries))):
            raise ValueError(f"{chrom} 的 multiinter 分段存在重叠")
        indices[chrom] = (entries, [e[0] for e in entries], [e[1] for e in entries])

    support_rows = []
    for _, peak in result.drop_duplicates("peak_id").iterrows():
        chrom = peak["chr"]
        left, right = int(peak["peak_start_0based"]), int(peak["peak_end"])
        width = right - left
        per_ip = [0] * 5
        supported_bp = 0
        max_support = 0
        if chrom in indices:
            entries, starts, ends = indices[chrom]
            first, last = bisect_right(ends, left), bisect_left(starts, right)
            for a, b, support, flags in entries[first:last]:
                overlap = min(right, b) - max(left, a)
                if overlap > 0:
                    supported_bp += overlap
                    max_support = max(max_support, support)
                    for i, flag in enumerate(flags):
                        per_ip[i] += overlap * flag
        row = {
            "peak_id": peak["peak_id"],
            "ge2_supported_bp": supported_bp,
            "ge2_supported_fraction": supported_bp / width,
            "max_simultaneous_IP_support": max_support,
            "n_IP_with_any_supported_overlap": sum(bp > 0 for bp in per_ip),
            "IP_replicates_with_any_supported_overlap": ",".join(
                sample for sample, bp in zip(REP_SAMPLES, per_ip) if bp > 0
            ),
        }
        for sample, bp in zip(REP_SAMPLES, per_ip):
            row[f"{sample}__supported_bp"] = bp
            row[f"{sample}__supported_fraction"] = bp / width
        support_rows.append(row)

    result = result.merge(pd.DataFrame(support_rows), on="peak_id", validate="many_to_one")

    # 4. Keep every original per-peak column; prefixes distinguish the source table.
    target_ids = set(result["peak_id"])
    for label, table in tables.items():
        if not target_ids.issubset(set(table["peak_id"])):
            raise ValueError(f"{label} 中缺少目标 peak，请核对分析文件是否来自同一次运行")
        renamed = table.rename(columns={
            name: f"{label}__{name}" for name in table.columns if name != "peak_id"
        })
        result = result.merge(renamed, on="peak_id", how="left", validate="many_to_one")

    result.sort_values(["gene_name", "chr", "peak_start_0based"]).to_csv(
        out / "target_peak_evidence.tsv", sep="\t", index=False
    )
    print(overview[["gene_name", "chr", "strand", "n_overlapping_consensus_peaks"]].to_string(index=False))
    print(f"\n完成，结果目录：{out}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--targets", nargs="+", default=TARGETS)
    parser.add_argument("--ip-order", type=Path, help="One sample_id per line; exact multiinter column order")
    args = parser.parse_args()
    TARGETS = args.targets
    if args.ip_order:
        REP_SAMPLES = [x.strip() for x in args.ip_order.read_text().splitlines() if x.strip()]
        if len(REP_SAMPLES) != 5 or len(set(REP_SAMPLES)) != 5:
            parser.error("--ip-order must list exactly five unique IP IDs")
    main(args.gtf, args.root)
