#!/usr/bin/env python3
"""Build each poster panel as a separate high-resolution figure."""
from __future__ import annotations

import ast
import csv
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

BASE = Path("/home/theoa/yup/panglofi_results")
HYB = BASE / "i3_ST258_PRJNA1133668_rustaroo_iqtree" / "results"
SHORT = BASE / "i3_ST258_PRJNA1133668_rustaroo_iqtree_short" / "results"
OUT = BASE / "poster_figures"
OUT.mkdir(exist_ok=True)
FONT_DIR = OUT / "fonts"

# Poster-friendly palette (not purple-on-white / cream-serif defaults)
HYBRID = "#0B6E4F"  # deep green
SHORT_C = "#C45C26"  # terracotta/rust
SHARED = "#2F4858"  # slate
HYB_ONLY = "#0B6E4F"
SHORT_ONLY = "#C45C26"
GRID = "#4A4A4A"
TEXT = "#1C1C1C"
MUTED = "#5C5C5C"

# Shared canvas + type scale for every panel
FIGSIZE = (10.0, 4.82)  # 5% taller than 4.59
DPI = 300
PAD_INCHES = 0.25
FS_SUP = 18
FS_TITLE = 17
FS_LABEL = 16
FS_TICK = 15
FS_LEGEND = 15
FS_ANNOT = 14


def apply_rc():
    # Prefer Helvetica; fall back to bundled Nimbus Sans (Helvetica clone), then Arial.
    for otf in sorted(FONT_DIR.glob("*.otf")) if FONT_DIR.exists() else []:
        font_manager.fontManager.addfont(str(otf))
    plt.rcParams.update(
        {
            "font.size": FS_TICK,
            "axes.titlesize": FS_TITLE,
            "axes.labelsize": FS_LABEL,
            "xtick.labelsize": FS_TICK,
            "ytick.labelsize": FS_TICK,
            "legend.fontsize": FS_LEGEND,
            "figure.titlesize": FS_SUP,
            "axes.titleweight": "medium",
            "figure.titleweight": "medium",
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Helvetica",
                "Arial",
                "Nimbus Sans",
                "NimbusSans",
                "TeX Gyre Heros",
                "Liberation Sans",
                "DejaVu Sans",
            ],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def new_figure(nrows=1, ncols=2):
    fig, axes = plt.subplots(nrows, ncols, figsize=FIGSIZE, constrained_layout=True)
    fig.patch.set_facecolor("none")
    fig.patch.set_alpha(0.0)
    return fig, axes


def style_axes(ax):
    ax.set_facecolor("none")
    ax.patch.set_alpha(0.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#B0B0B0")
    ax.spines["bottom"].set_color("#B0B0B0")
    ax.tick_params(colors=TEXT, labelsize=FS_TICK)
    ax.yaxis.grid(True, color=GRID, linewidth=0.9, zorder=0)
    ax.set_axisbelow(True)


def panel_title(ax, title: str):
    ax.set_title(title, loc="left", fontsize=FS_TITLE, fontweight=500, color=TEXT, pad=8)


def axis_label(ax, which: str, text: str):
    if which == "x":
        ax.set_xlabel(text, fontsize=FS_LABEL, color=TEXT)
    else:
        ax.set_ylabel(text, fontsize=FS_LABEL, color=TEXT)


def fig_title(fig, title: str):
    fig.suptitle(title, fontsize=FS_SUP, fontweight=500, color=TEXT)


def classify_edit(delete_len: int, inserted: str) -> str:
    ins = len(inserted or "")
    delete_len = int(delete_len)
    if delete_len == 1 and ins == 1:
        return "snp"
    if delete_len == 0 and ins > 0:
        return "insertion"
    if delete_len > 0 and ins == 0:
        return "deletion"
    if delete_len == ins and delete_len > 1:
        return "mnp"
    if delete_len == 0 and ins == 0:
        return "noop"
    return "complex_indel"


def indel_frame(delete_len: int, inserted: str) -> str | None:
    ins = len(inserted or "")
    delete_len = int(delete_len)
    if delete_len == ins:
        return None
    return "inframe" if abs(delete_len - ins) % 3 == 0 else "frameshift"


def parse_codes(cell: str):
    cell = (cell or "").strip()
    if cell in {"", "[]"}:
        return []
    try:
        parsed = ast.literal_eval(cell)
    except Exception:
        return None
    return parsed if isinstance(parsed, list) else None


def load_lof(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def analyze(rows, replicon_mode: str):
    """replicon_mode: 'chromosome' or 'plasmid'."""
    if replicon_mode == "chromosome":
        rows = [r for r in rows if r.get("Replicon") == "chromosome"]
    else:
        rows = [r for r in rows if r.get("Replicon") and r.get("Replicon") != "chromosome"]

    edit_counts = Counter()
    frame_counts = Counter()
    total_ref_bp = 0
    events = defaultdict(set)  # class -> {(sample, gene)}
    edits_per_sample = Counter()

    for row in rows:
        cell = row.get("Mutation_Codes", "")
        if cell.strip() == "":
            continue
        try:
            ref_len = int(row.get("Reference_Length") or 0)
        except ValueError:
            ref_len = 0
        total_ref_bp += ref_len
        codes = parse_codes(cell)
        if not codes:
            continue
        sample = row["Sample"].split("__", 1)[0]
        gene = row["Gene"]
        for edit in codes:
            if not (isinstance(edit, (list, tuple)) and len(edit) >= 3):
                continue
            _, delete_len, inserted = edit[0], edit[1], edit[2]
            et = classify_edit(delete_len, inserted or "")
            if et == "noop":
                continue
            edit_counts[et] += 1
            edits_per_sample[sample] += 1
            # poster classes
            if et == "snp":
                events["snp"].add((sample, gene))
            elif et == "mnp":
                events["mnp"].add((sample, gene))
            elif et in {"insertion", "deletion", "complex_indel"}:
                events["indel"].add((sample, gene))
            fr = indel_frame(delete_len, inserted or "")
            if fr == "frameshift":
                frame_counts["frameshift"] += 1
                events["frameshift"].add((sample, gene))
            elif fr == "inframe":
                frame_counts["inframe"] += 1

    # rates for poster classes
    indel_n = (
        edit_counts["insertion"] + edit_counts["deletion"] + edit_counts["complex_indel"]
    )
    rates = {
        "SNP": 1000 * edit_counts["snp"] / total_ref_bp if total_ref_bp else 0,
        "MNP": 1000 * edit_counts["mnp"] / total_ref_bp if total_ref_bp else 0,
        "Indel": 1000 * indel_n / total_ref_bp if total_ref_bp else 0,
        "Frameshift": 1000 * frame_counts["frameshift"] / total_ref_bp if total_ref_bp else 0,
    }
    counts = {
        "SNP": edit_counts["snp"],
        "MNP": edit_counts["mnp"],
        "Indel": indel_n,
        "Frameshift": frame_counts["frameshift"],
    }
    return {
        "rates": rates,
        "counts": counts,
        "events": events,
        "edits_per_sample": edits_per_sample,
        "total_ref_bp": total_ref_bp,
        "samples": sorted(edits_per_sample.keys() | set()),
    }


def concordance(hyb_events, short_events, key):
    a = hyb_events.get(key, set())
    b = short_events.get(key, set())
    return len(a & b), len(a - b), len(b - a)


def parse_bandage(path: Path) -> dict[str, float]:
    data = {}
    for line in path.read_text().splitlines():
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        data[k.strip()] = float(v.strip().replace("%", "").split()[0])
    return {
        "nodes": data["Node count"],
        "n50": data["N50 (bp)"],
    }


def load_bandage(root: Path):
    out = {}
    for p in (root / "bandage").glob("*_bandage.info.txt"):
        sid = p.name.replace("_bandage.info.txt", "")
        out[sid] = parse_bandage(p)
    return out


def save(fig, name: str):
    # Fixed canvas size so every panel exports at the same pixel dimensions.
    fig.set_size_inches(*FIGSIZE)
    png = OUT / f"{name}.png"
    pdf = OUT / f"{name}.pdf"
    for path in (png, pdf):
        fig.savefig(
            path,
            dpi=DPI,
            facecolor="none",
            edgecolor="none",
            transparent=True,
            bbox_inches="tight",
            pad_inches=PAD_INCHES,
        )
    plt.close(fig)
    print(f"wrote {png}")
    print(f"wrote {pdf}")


def fig1a(hyb_chrom, short_chrom, hyb_plas, short_plas):
    classes = ["SNP", "MNP", "Indel", "Frameshift"]
    x = np.arange(len(classes))
    width = 0.35

    fig, axes = new_figure(1, 2)
    panels = [
        (axes[0], hyb_chrom, short_chrom, "Chromosome"),
        (axes[1], hyb_plas, short_plas, "Plasmid"),
    ]
    for ax, hyb, short, title in panels:
        style_axes(ax)
        hvals = [hyb["rates"][c] for c in classes]
        svals = [short["rates"][c] for c in classes]
        ax.bar(
            x - width / 2,
            hvals,
            width,
            label="Hybrid",
            color=HYBRID,
            edgecolor="none",
            zorder=3,
        )
        ax.bar(
            x + width / 2,
            svals,
            width,
            label="Short-only",
            color=SHORT_C,
            edgecolor="none",
            zorder=3,
        )
        ax.set_xticks(x)
        ax.set_xticklabels(classes)
        axis_label(ax, "x", title)
        axis_label(ax, "y", "Edits / kb CDS")

    axes[0].legend(frameon=False, fontsize=FS_LEGEND, loc="upper right")
    save(fig, "fig1a_mutation_rates")


def fig1b(hyb_chrom, short_chrom, hyb_plas, short_plas):
    classes = ["SNP", "MNP", "Indel", "Frameshift"]
    keys = ["snp", "mnp", "indel", "frameshift"]

    fig, axes = new_figure(1, 2)
    for ax, hyb, short, title in [
        (axes[0], hyb_chrom, short_chrom, "Chromosome"),
        (axes[1], hyb_plas, short_plas, "Plasmid"),
    ]:
        style_axes(ax)
        shared, hyb_only, short_only = [], [], []
        jaccards = []
        for key in keys:
            s, h, o = concordance(hyb["events"], short["events"], key)
            shared.append(s)
            hyb_only.append(h)
            short_only.append(o)
            union = s + h + o
            jaccards.append(s / union if union else 0)

        x = np.arange(len(classes))
        ax.bar(x, shared, color=SHARED, label="Shared", zorder=3)
        ax.bar(x, hyb_only, bottom=shared, color=HYB_ONLY, label="Hybrid-only", zorder=3)
        bottoms = [a + b for a, b in zip(shared, hyb_only)]
        ax.bar(x, short_only, bottom=bottoms, color=SHORT_ONLY, label="Short-only", zorder=3)

        for i, j in enumerate(jaccards):
            top = shared[i] + hyb_only[i] + short_only[i]
            ax.text(
                i,
                top * 1.02 if top else 0.5,
                f"J={j:.2f}",
                ha="center",
                va="bottom",
                fontsize=FS_ANNOT,
                color=MUTED,
            )

        ax.set_xticks(x)
        ax.set_xticklabels(classes)
        axis_label(ax, "x", title)
        axis_label(ax, "y", "Sample–gene events")
        ymax = max(shared[i] + hyb_only[i] + short_only[i] for i in range(len(classes)))
        ax.set_ylim(0, ymax * 1.18 if ymax else 1)

    axes[0].legend(frameon=False, fontsize=FS_LEGEND, loc="upper right")
    save(fig, "fig1b_concordance")


def fig1c(hyb_b, short_b):
    common = sorted(set(hyb_b) & set(short_b))
    hyb_n50 = [hyb_b[s]["n50"] / 1e6 for s in common]
    short_n50 = [short_b[s]["n50"] / 1e6 for s in common]
    hyb_nodes = [hyb_b[s]["nodes"] for s in common]
    short_nodes = [short_b[s]["nodes"] for s in common]

    fig, axes = new_figure(1, 2)
    rng = np.random.default_rng(0)

    ax = axes[0]
    style_axes(ax)
    vp = ax.violinplot(
        [hyb_n50, short_n50],
        positions=[0, 1],
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )
    for body, color in zip(vp["bodies"], [HYBRID, SHORT_C]):
        body.set_facecolor(color)
        body.set_alpha(0.85)
        body.set_edgecolor("none")
    vp["cmedians"].set_color(TEXT)
    vp["cmedians"].set_linewidth(1.4)
    for i, vals, color in [(0, hyb_n50, HYBRID), (1, short_n50, SHORT_C)]:
        jitter = rng.uniform(-0.08, 0.08, size=len(vals))
        ax.scatter(
            np.full(len(vals), i) + jitter,
            vals,
            s=10,
            color="white",
            edgecolors=color,
            linewidths=0.6,
            alpha=0.7,
            zorder=4,
        )
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Hybrid", "Short-only"])
    axis_label(ax, "y", "N50 (Mb)")

    ax = axes[1]
    style_axes(ax)
    vp = ax.violinplot(
        [hyb_nodes, short_nodes],
        positions=[0, 1],
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )
    for body, color in zip(vp["bodies"], [HYBRID, SHORT_C]):
        body.set_facecolor(color)
        body.set_alpha(0.85)
        body.set_edgecolor("none")
    vp["cmedians"].set_color(TEXT)
    vp["cmedians"].set_linewidth(1.4)
    for i, vals, color in [(0, hyb_nodes, HYBRID), (1, short_nodes, SHORT_C)]:
        jitter = rng.uniform(-0.08, 0.08, size=len(vals))
        ax.scatter(
            np.full(len(vals), i) + jitter,
            vals,
            s=10,
            color="white",
            edgecolors=color,
            linewidths=0.6,
            alpha=0.7,
            zorder=4,
        )
    ax.set_yscale("log")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Hybrid", "Short-only"])
    axis_label(ax, "y", "Graph nodes (log)")

    save(fig, "fig1c_assembly_contiguity")


def fig2(hyb_chrom, short_chrom):
    presence = HYB / "rustaroo" / "chromosome" / "rustroo" / "presence.tsv"
    all_samples = presence.read_text().splitlines()[0].split("\t")[1:]
    xs = [hyb_chrom["edits_per_sample"].get(s, 0) for s in all_samples]
    ys = [short_chrom["edits_per_sample"].get(s, 0) for s in all_samples]

    # Same canvas as the dual-panel figures; scatter fills the full width.
    fig, ax = new_figure(1, 1)
    style_axes(ax)

    lim = max(max(xs), max(ys)) * 1.05
    ax.plot([0, lim], [0, lim], color="#B0B0B0", linewidth=1.0, linestyle="--", zorder=1)
    ax.scatter(
        xs,
        ys,
        s=36,
        color=HYBRID,
        edgecolors="white",
        linewidths=0.6,
        alpha=0.85,
        zorder=3,
    )
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    axis_label(ax, "x", "Hybrid chromosome edits / sample")
    axis_label(ax, "y", "Short-only chromosome edits / sample")

    diffs = [(abs(y - x), s, x, y) for s, x, y in zip(all_samples, xs, ys)]
    diffs.sort(reverse=True)
    for _, s, x, y in diffs[:3]:
        ax.annotate(
            s.replace("SAMN", ""),
            (x, y),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=FS_ANNOT,
            color=MUTED,
        )

    ax.text(
        0.04,
        0.96,
        "n = 70 isolates",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=FS_ANNOT,
        color=MUTED,
    )
    save(fig, "fig2_per_sample_scatter")


def main():
    apply_rc()
    print("Loading LoF matrices…")
    hyb_lof = load_lof(
        HYB / "lof_matrix" / "ST258_PRJNA1133668_rustaroo_iqtree_lof_matrix.tsv"
    )
    short_lof = load_lof(SHORT / "lof_matrix" / "lof_matrix.tsv")

    print("Computing mutation summaries…")
    hyb_chrom = analyze(hyb_lof, "chromosome")
    short_chrom = analyze(short_lof, "chromosome")
    hyb_plas = analyze(hyb_lof, "plasmid")
    short_plas = analyze(short_lof, "plasmid")

    print("Rates chromosome hybrid", hyb_chrom["rates"])
    print("Rates chromosome short ", short_chrom["rates"])
    print("Rates plasmid hybrid", hyb_plas["rates"])
    print("Rates plasmid short ", short_plas["rates"])

    fig1a(hyb_chrom, short_chrom, hyb_plas, short_plas)
    fig1b(hyb_chrom, short_chrom, hyb_plas, short_plas)

    hyb_b = load_bandage(HYB)
    short_b = load_bandage(SHORT)
    fig1c(hyb_b, short_b)
    fig2(hyb_chrom, short_chrom)

    # one-line caption file for convenience
    caption = OUT / "figure_captions.txt"
    caption.write_text(
        "Fig 1A. Mutation rates by edit type (edits per kb CDS) for chromosome and plasmid, "
        "hybrid vs short-only assemblies.\n"
        "Fig 1B. Concordance of sample–gene mutation events (shared / hybrid-only / short-only) "
        "with Jaccard index; rates can match while call identity diverges.\n"
        "Fig 1C. Assembly contiguity (N50 and graph nodes) explaining why mutation calls change.\n"
        "Fig 2. Per-sample chromosome edit counts; points off the diagonal are isolates whose "
        "mutation burden shifts most between hybrid and short-only assemblies.\n"
        "One-line result: Hybrid assembly barely changes chromosome mutation rates, but changes "
        "a large fraction of which calls are made — and reshapes plasmid mutation calls much more.\n"
    )
    print(f"wrote {caption}")


if __name__ == "__main__":
    main()
