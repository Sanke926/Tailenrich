from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "sample_input" / "MSBB-BM36"
RESULT_FILE = ROOT / "results" / "selected_genes.tsv"
OUT_FILE = ROOT / "figures" / "MSBB-BM36_selected_pe_curves.png"

COLOR_LEFT = "#00798c"
COLOR_RIGHT = "#d1495b"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_csv(path: Path) -> list[list[str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.reader(fh))


def load_expression() -> tuple[dict[str, list[float]], list[str]]:
    rows = read_csv(INPUT_DIR / "gene_TPM_by_salmon_covAdjusted.csv")
    header = rows[0]
    sample_ids = header[1:]
    expr = {}
    for row in rows[1:]:
        expr[row[0]] = [float(x) for x in row[1:]]
    return expr, sample_ids


def load_groups(sample_ids: list[str]) -> list[int]:
    rows = read_tsv(INPUT_DIR / "used_samples_group.tsv")
    group_map = {row["sampleID"]: int(float(row["group"])) for row in rows}
    return [group_map[sample_id] for sample_id in sample_ids]


def compute_curve(values: list[float], groups: list[int], scan_mode: str) -> tuple[list[float], list[float], int]:
    pairs = list(zip(values, groups))
    reverse = scan_mode == "H2L"
    pairs.sort(key=lambda x: x[0], reverse=reverse)

    y_sorted = [float(g) for _, g in pairs]
    mean_y = sum(y_sorted) / len(y_sorted)
    pos_sum = sum(v - mean_y for v in y_sorted if v - mean_y > 0)
    running = []
    cumulative = 0.0
    for g in y_sorted:
        cumulative += g - mean_y
        running.append(0.0 if pos_sum == 0 else cumulative / pos_sum)

    xvals = [idx / (len(running) - 1) for idx in range(len(running))]
    peak_idx = max(range(len(running)), key=lambda i: running[i])

    # Mirror the display for right-tail enrichment so the enriched tail sits
    # on the right side of the panel while preserving the original curve shape.
    if scan_mode == "H2L":
        xvals = [1.0 - x for x in xvals][::-1]
        running = running[::-1]
        peak_idx = len(running) - 1 - peak_idx

    return xvals, running, peak_idx


def add_dimension_marks(ax, xvals: list[float], yvals: list[float], peak_idx: int, color: str, direction: str) -> None:
    x_peak = xvals[peak_idx]
    y_peak = yvals[peak_idx]
    ymin = min(yvals)
    y_abs = max(abs(y) for y in yvals) or 1.0
    x_arrow_y = min(-0.06 * y_abs, ymin - 0.08 * y_abs)

    ax.annotate("", xy=(x_peak, y_peak), xytext=(x_peak, 0.0), arrowprops=dict(arrowstyle="<->", color=color, lw=1.1))
    if direction == "left":
        x_start, x_end = 0.0, x_peak
    else:
        x_start, x_end = x_peak, 1.0
    ax.annotate("", xy=(x_end, x_arrow_y), xytext=(x_start, x_arrow_y), arrowprops=dict(arrowstyle="<->", color=color, lw=1.1))
    ax.text(x_peak + 0.02, y_peak / 2 if y_peak != 0 else 0.02, "h", color=color, fontsize=9, va="center")
    ax.text((x_start + x_end) / 2, x_arrow_y - 0.04 * y_abs, "x", color=color, fontsize=9, ha="center", va="top")


def panel_label(scan_mode: str) -> str:
    return "left scan" if scan_mode == "L2H" else "right scan"


def enriched_label(direction: str) -> str:
    return "left-tail enriched gene" if direction == "left" else "right-tail enriched gene"


def main() -> None:
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    expr_map, sample_ids = load_expression()
    groups = load_groups(sample_ids)
    selected_rows = read_tsv(RESULT_FILE)

    fig, axes = plt.subplots(len(selected_rows), 2, figsize=(11.8, 13.8))

    for row_idx, row in enumerate(selected_rows):
        gene = row["gene_id"]
        direction = row["direction"]
        row_label = enriched_label(direction)

        for col_idx, scan_mode in enumerate(("L2H", "H2L")):
            ax = axes[row_idx, col_idx]
            scan_direction = "left" if scan_mode == "L2H" else "right"
            color = COLOR_LEFT if scan_direction == "left" else COLOR_RIGHT
            xvals, yvals, peak_idx = compute_curve(expr_map[gene], groups, scan_mode)

            ax.plot(xvals, yvals, color=color, lw=2.2)
            ax.axhline(0.0, color="#7c8a9f", linestyle="--", lw=1)
            ax.scatter([xvals[peak_idx]], [yvals[peak_idx]], color=color, s=28, zorder=3)
            add_dimension_marks(ax, xvals, yvals, peak_idx, color, scan_direction)

            y_abs = max(abs(y) for y in yvals) or 1.0
            ax.set_xlim(0, 1.03)
            ax.set_ylim(min(min(yvals) - 0.18 * y_abs, -0.18 * y_abs), max(yvals) + 0.18 * y_abs)
            ax.set_title(f"{gene}\n{row_label} | {panel_label(scan_mode)}", fontsize=10)
            ax.set_xlabel("Ordered sample fraction", fontsize=9)
            ax.set_ylabel("Running score", fontsize=9)
            ax.text(
                0.02,
                0.96,
                f"log2FC = {float(row['tailenrich_log2FC']):.3f}\nTailEnrich FDR = {float(row['tailenrich_fdr']):.4g}",
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=8.3,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#d9d9d9", alpha=0.9),
            )
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

    fig.suptitle("Selected PE curves in both scan directions", fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
