#!/usr/bin/env python3
"""
erbb2_sv_arcs.py

erbb2_sv_compare.py의 long-format TSV를 가지고
AmpliconArchitect amplicon.png 느낌의 "염색체 위치 트랙 + breakpoint를 곡선으로 연결"
다이어그램을 그린다.

두 가지 그림:
1) {out_prefix}_arcs_summary.png
   - 염색체별 horizontal track + 모든 ERBB2-linked breakpoint를 곡선(arc)으로 연결
   - 곡선 색=진할수록 더 많은 클론에서 검출(truncal에 가까움), 연할수록 일부 클론만(subclonal)
2) {out_prefix}_arcs_facet.png
   - 샘플(parental + 지정한 클론들)별로 작은 패널 나눠서, 그 샘플에 있는 junction만 표시
   - 클론 간 구조 차이를 한눈에 비교

사용:
python3 erbb2_sv_arcs.py --long-tsv H2170_erbb2_sv_long.tsv --out-prefix H2170_erbb2
# 클론이 너무 많으면 facet에 보여줄 샘플을 직접 골라도 됨 (parental은 항상 포함됨)
python3 erbb2_sv_arcs.py --long-tsv H2170_erbb2_sv_long.tsv --out-prefix H2170_erbb2 \
  --facet-samples NCI_9,NCI_3,NCI_22,NCI_31
"""

import argparse
import re
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.cm as cm


def chrom_sort_key(c):
    c2 = c[3:] if c.lower().startswith("chr") else c
    if c2.isdigit():
        return (0, int(c2))
    order = {"X": 100, "Y": 101, "M": 102, "MT": 102}
    return (1, order.get(c2.upper(), 999), c2)


def build_tracks(df, pad_frac=0.1, min_window_frac=0.04):
    """chrom -> (y_pos, min_pos, max_pos) 결정.
    breakpoint가 1개뿐인 염색체는 span=0이 되어 트랙이 안 보이므로,
    전체 x축 범위의 min_window_frac 만큼 최소 길이를 보장한다."""
    spans = {}
    for _, r in df.iterrows():
        for c, p in [(r["chrom1"], r["pos1"]), (r["chrom2"], r["pos2"])]:
            lo, hi = spans.get(c, (p, p))
            spans[c] = (min(lo, p), max(hi, p))

    chroms = sorted(spans.keys(), key=chrom_sort_key)

    # 1차: pad만 적용
    raw = {}
    for c in chroms:
        lo, hi = spans[c]
        span = max(hi - lo, 1)
        pad = span * pad_frac
        raw[c] = (lo - pad, hi + pad)

    global_lo = min(v[0] for v in raw.values())
    global_hi = max(v[1] for v in raw.values())
    min_window = max((global_hi - global_lo) * min_window_frac, 1.0)

    tracks = {}
    for i, c in enumerate(chroms):
        lo, hi = raw[c]
        if (hi - lo) < min_window:
            center = (lo + hi) / 2
            lo, hi = center - min_window / 2, center + min_window / 2
        tracks[c] = {"y": i, "lo": lo, "hi": hi}
    return tracks


def draw_arc(ax, x1, y1, x2, y2, color="#2c5282", lw=1.5, alpha=0.8, zorder=2):
    """y1==y2(같은 염색체)면 위로 솟는 반원형 arc, 다르면 S자 곡선으로 연결"""
    if y1 == y2:
        mid_y = y1 + 0.35
        verts = [(x1, y1), (x1, mid_y), (x2, mid_y), (x2, y2)]
    else:
        mid_y = (y1 + y2) / 2
        verts = [(x1, y1), (x1, mid_y), (x2, mid_y), (x2, y2)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    patch = PathPatch(Path(verts, codes), facecolor="none", edgecolor=color,
                       lw=lw, alpha=alpha, zorder=zorder)
    ax.add_patch(patch)


def draw_tracks(ax, tracks):
    for c, t in tracks.items():
        ax.plot([t["lo"], t["hi"]], [t["y"], t["y"]], color="#444444", lw=4,
                solid_capstyle="round", zorder=1)
        ax.text(t["lo"], t["y"] - 0.25, c, fontsize=9, ha="left", va="top")


def mb(x):
    return x / 1_000_000


def plot_summary(df, out_path):
    bp_summary = (
        df[df["sample_type"] == "clone"]
        .groupby("bp_id")
        .apply(lambda g: (g["status"].isin(["detected_filtered", "detected_unfiltered"])).mean())
        .rename("frac_clones_detected")
        .reset_index()
    )
    meta = df.drop_duplicates("bp_id")[["bp_id", "chrom1", "pos1", "chrom2", "pos2"]]
    bp_summary = bp_summary.merge(meta, on="bp_id")

    tracks = build_tracks(bp_summary)
    fig_h = max(3, 1.0 * len(tracks) + 1)
    fig, ax = plt.subplots(figsize=(10, fig_h))
    draw_tracks(ax, tracks)

    cmap = matplotlib.colormaps["Blues"]
    for _, r in bp_summary.iterrows():
        t1, t2 = tracks[r["chrom1"]], tracks[r["chrom2"]]
        frac = r["frac_clones_detected"]
        color = cmap(0.3 + 0.7 * frac)
        draw_arc(ax, r["pos1"], t1["y"], r["pos2"], t2["y"],
                  color=color, lw=1 + 3 * frac, alpha=0.85)

    ax.set_ylim(-1, len(tracks))
    ax.set_xlabel("genomic position")
    ax.set_yticks([])
    ax.set_title("ERBB2-linked SV breakpoints (parental-defined)\n"
                  "darker/thicker arc = detected in more clones (closer to truncal)",
                  fontsize=10)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(0, 1))
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("fraction of clones detected")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[done] {out_path}", file=sys.stderr)


def plot_facet(df, out_path, samples):
    meta = df.drop_duplicates("bp_id")[["bp_id", "chrom1", "pos1", "chrom2", "pos2"]]
    tracks = build_tracks(meta)

    n = len(samples)
    ncols = min(4, n)
    nrows = -(-n // ncols)
    fig_h = max(3, 1.0 * len(tracks) + 0.5)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, fig_h * nrows), squeeze=False)

    for idx, s in enumerate(samples):
        ax = axes[idx // ncols][idx % ncols]
        draw_tracks(ax, tracks)
        sub = df[df["sample_id"] == s]
        for _, r in sub.iterrows():
            detected = r["status"] in ("source", "detected_filtered", "detected_unfiltered")
            if not detected:
                continue
            t1, t2 = tracks[r["chrom1"]], tracks[r["chrom2"]]
            color = "#742a2a" if r["status"] == "source" else "#2c5282"
            lw = 2.2 if r["status"] == "detected_filtered" or r["status"] == "source" else 1.2
            alpha = 0.9 if r["status"] != "detected_unfiltered" else 0.5
            draw_arc(ax, r["pos1"], t1["y"], r["pos2"], t2["y"], color=color, lw=lw, alpha=alpha)
        ax.set_ylim(-1, len(tracks))
        ax.set_yticks([])
        ax.set_title(s, fontsize=10)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle("ERBB2-linked SV breakpoints per sample\n"
                  "(red=parental, dark blue=detected_filtered, light blue=detected_unfiltered only)",
                  fontsize=11, y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[done] {out_path}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--long-tsv", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--facet-samples", default=None,
                     help="콤마로 구분된 클론 sample_id 목록 (없으면 전체 클론, 너무 많으면 4개씩 자동 wrap). "
                          "parental은 자동으로 포함됨")
    args = ap.parse_args()

    df = pd.read_csv(args.long_tsv, sep="\t")
    if df.empty:
        sys.exit("[error] 입력 TSV가 비어있습니다.")

    plot_summary(df, f"{args.out_prefix}_arcs_summary.png")

    if args.facet_samples:
        clones = [s.strip() for s in args.facet_samples.split(",") if s.strip()]
    else:
        clones = [s for s in df["sample_id"].unique() if s != "parental"]
    samples = ["parental"] + clones
    plot_facet(df, f"{args.out_prefix}_arcs_facet.png", samples)


if __name__ == "__main__":
    main()
