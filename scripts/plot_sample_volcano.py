from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
TAILENRICH_FILE = ROOT / "sample_output" / "MSBB-BM36" / "tailEnrich.tsv"
OUT_FILE = ROOT / "figures" / "MSBB-BM36_volcano.png"
LOG2FC_CUTOFF = math.log2(1.5)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    te_rows = read_tsv(TAILENRICH_FILE)

    xs_bg = []
    ys_bg = []
    xs_hit = []
    ys_hit = []

    for row in te_rows:
        log2fc = float(row["log2FC"]) if row["log2FC"] not in ("", "NA", "NaN") else 0.0
        fdr = float(row["FDR"])
        yval = -math.log10(max(fdr, 1e-300))
        if fdr < 0.05 and abs(log2fc) >= LOG2FC_CUTOFF:
            xs_hit.append(log2fc)
            ys_hit.append(yval)
        else:
            xs_bg.append(log2fc)
            ys_bg.append(yval)

    plt.figure(figsize=(8.6, 6.0))
    plt.scatter(xs_bg, ys_bg, s=10, c="#d9d9d9", alpha=0.55, linewidths=0, label="Other genes")
    plt.scatter(xs_hit, ys_hit, s=16, c="#d1495b", alpha=0.8, linewidths=0, label="TailEnrich significant")
    plt.axhline(-math.log10(0.05), color="#465c69", linestyle=":", linewidth=1)
    plt.axvline(LOG2FC_CUTOFF, color="#465c69", linestyle="--", linewidth=1)
    plt.axvline(-LOG2FC_CUTOFF, color="#465c69", linestyle="--", linewidth=1)

    plt.xlabel("log2FC")
    plt.ylabel("-log10(FDR)")
    plt.title("MSBB-BM36 volcano plot")
    plt.text(
        0.98,
        0.02,
        "|log2FC| >= log2(1.5)",
        transform=plt.gca().transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="#465c69",
    )
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(OUT_FILE, dpi=180)
    plt.close()


if __name__ == "__main__":
    main()
