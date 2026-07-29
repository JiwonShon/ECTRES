#!/usr/bin/env python3
"""
erbb2_sv_plot.py

erbb2_sv_compare.py가 만든 long-format TSV(bp_id, sample_id, status, AD, DP, VAF ...)를
가지고 두 가지 그림을 그린다:

1) {out_prefix}_clustermap.png
   - presence/absence (detected_filtered/detected_unfiltered = 1, not_detected = 0) 기반
   - 클론(parental 제외)들을 Jaccard distance로 클러스터링한 clustermap
   - "어떤 클론들이 비슷한 SV set을 갖고 있는지" (= 비슷한 시기에 일어난 이벤트인지)를 직접 보여줌

2) {out_prefix}_vaf_heatmap.png
   - parental + 클론(1번 클러스터링 순서대로 정렬) 전체를 annotated heatmap으로
   - 셀 안에 VAF 값 표시, AD>DP인 경우(=신뢰 불가) 빨간 테두리로 표시
   - not_detected는 회색으로 표시

사용:
python3 erbb2_sv_plot.py --long-tsv H2170_erbb2_sv_long.tsv --out-prefix H2170_erbb2
"""

import argparse
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist


def bp_label(row):
    return f"{row['chrom1']}:{row['pos1']:,}-{row['chrom2']}:{row['pos2']:,}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--long-tsv", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.long_tsv, sep="\t")
    df["bp_label"] = df.apply(bp_label, axis=1)

    # numeric 변환 (NA는 그대로 NaN)
    for c in ["AD", "DP"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["VAF_num"] = pd.to_numeric(df["VAF"], errors="coerce")
    df["unreliable"] = (df["AD"] > df["DP"]).fillna(False)

    samples = list(df["sample_id"].unique())
    if "parental" not in samples:
        sys.exit("[error] sample_id == 'parental' 행이 없습니다. erbb2_sv_compare.py 결과인지 확인하세요.")
    clones = [s for s in samples if s != "parental"]

    bp_labels = list(df["bp_label"].unique())
    if len(bp_labels) == 0:
        sys.exit("[error] breakpoint가 0개입니다. 매칭 단계부터 다시 확인하세요.")

    # ---------- 1) presence/absence matrix (clone만) ----------
    detected = df["status"].isin(["detected_filtered", "detected_unfiltered"]).astype(int)
    df["detected"] = detected

    pres_abs = (
        df[df["sample_id"] != "parental"]
        .pivot_table(index="bp_label", columns="sample_id", values="detected", fill_value=0)
        .reindex(index=bp_labels)
    )

    if len(clones) >= 2 and len(bp_labels) >= 2:
        try:
            g = sns.clustermap(
                pres_abs, cmap=["#e8e8e8", "#2b6cb0"], cbar_pos=None,
                figsize=(max(8, 0.35 * len(clones)), max(6, 0.4 * len(bp_labels))),
                linewidths=0.4, linecolor="white",
                xticklabels=True, yticklabels=True,
                metric="jaccard", method="average",
            )
            g.ax_heatmap.set_xlabel("clone")
            g.ax_heatmap.set_ylabel("ERBB2-linked breakpoint")
            g.fig.suptitle("Presence/absence of ERBB2-linked SV breakpoints across clones\n"
                            "(clustered by Jaccard distance -> similar columns = similar event timing)",
                            y=1.02, fontsize=10)
            clone_order = [t.get_text() for t in g.ax_heatmap.get_xticklabels()]
            g.savefig(f"{args.out_prefix}_clustermap.png", dpi=200, bbox_inches="tight")
            plt.close(g.fig)
            print(f"[done] {args.out_prefix}_clustermap.png", file=sys.stderr)
        except Exception as e:
            print(f"[warn] clustermap 생성 실패({e}), 클러스터링 없이 원래 순서 사용", file=sys.stderr)
            clone_order = clones
    else:
        print("[warn] 클론 또는 breakpoint가 2개 미만이라 클러스터링을 스킵합니다.", file=sys.stderr)
        clone_order = clones

    # ---------- 2) parental + clone(clustered order) VAF annotated heatmap ----------
    sample_order = ["parental"] + clone_order
    status_code = {"source": 3, "detected_filtered": 2, "detected_unfiltered": 1, "not_detected": 0}
    df["status_code"] = df["status"].map(status_code).fillna(0)

    code_mat = df.pivot_table(index="bp_label", columns="sample_id", values="status_code", fill_value=0)
    code_mat = code_mat.reindex(index=bp_labels, columns=sample_order)

    vaf_mat = df.pivot_table(index="bp_label", columns="sample_id", values="VAF_num")
    vaf_mat = vaf_mat.reindex(index=bp_labels, columns=sample_order)

    unrel_mat = df.pivot_table(index="bp_label", columns="sample_id", values="unreliable", aggfunc="max", fill_value=False)
    unrel_mat = unrel_mat.reindex(index=bp_labels, columns=sample_order).fillna(False)

    annot = vaf_mat.copy().astype(object)
    for r in annot.index:
        for c in annot.columns:
            v = vaf_mat.loc[r, c]
            annot.loc[r, c] = "NA" if pd.isna(v) else f"{v:.2f}"

    fig_w = max(8, 0.6 * len(sample_order))
    fig_h = max(4, 0.5 * len(bp_labels))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = matplotlib.colors.ListedColormap(["#e8e8e8", "#bee3f8", "#63b3ed", "#2c5282"])
    sns.heatmap(code_mat, cmap=cmap, vmin=0, vmax=3, cbar=False, annot=annot, fmt="",
                linewidths=0.5, linecolor="white", ax=ax)

    # AD>DP (불신뢰) 셀에 빨간 테두리
    for i, r in enumerate(unrel_mat.index):
        for j, c in enumerate(unrel_mat.columns):
            if unrel_mat.loc[r, c]:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor="red", lw=2))

    ax.set_title("ERBB2-linked SV breakpoints: parental vs clones\n"
                  "(color: dark=parental-source, blue shades=detected, grey=not_detected; "
                  "red box = AD>DP, VAF unreliable)", fontsize=10)
    ax.set_xlabel("sample (parental first, then clustered clone order)")
    ax.set_ylabel("breakpoint")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    fig.savefig(f"{args.out_prefix}_vaf_heatmap.png", dpi=200, bbox_inches="tight")
    print(f"[done] {args.out_prefix}_vaf_heatmap.png", file=sys.stderr)


if __name__ == "__main__":
    main()
