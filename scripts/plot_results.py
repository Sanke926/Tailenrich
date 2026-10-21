"""Plot each gene's selected enrichment direction with an independent color strip."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from real_gene_colorstrip import (
    COLOR_AD, COLOR_CONTROL, add_external_group_color_strip,
    normalize_curve_x, normalized_draw_curve,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("display_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Times New Roman", "svg.fonttype": "path"})
    genes = pd.read_csv(args.display_dir / "display_genes.tsv", sep="\t")
    curves = pd.read_csv(args.display_dir / "display_curves.tsv", sep="\t")
    audit = []
    for gene in genes.itertuples(index=False):
        direction = gene.direction
        fig, ax = plt.subplots(figsize=(6.8, 5.4))
        fig.subplots_adjust(left=0.14, right=0.97, bottom=0.19, top=0.90)
        data = curves[(curves.gene_id == gene.gene_id) &
                      (curves.direction == direction)].sort_values("position").copy()
        if data.empty or not np.isfinite(data.running_score).all():
            raise ValueError(f"Missing or invalid curve: {gene.gene_id} {direction}")
        peak_record = data.iloc[int(np.argmax(data.running_score.to_numpy()))]
        original_x = float(peak_record.scan_x)
        height = float(peak_record.running_score)
        pe = float(gene.enrichment_score)
        if not np.isclose(height * (1 - original_x), pe, atol=1e-12):
            raise ValueError("Curve statistic differs from archived result")
        data["group_code"] = data["group"]
        # Mirror right scans before assigning normalized sample-cell widths.
        data["sample_fraction_plot"] = data.scan_x if direction == "left" else 1 - data.scan_x
        normalized = normalize_curve_x(data)
        # Preserve the FIRST maximum in scan order, including tied maxima.
        peak = normalized.loc[normalized.position == peak_record.position].iloc[0]
        widths = normalized.cell_right - normalized.cell_left
        for group in (-1, 1):
            if not np.isclose(widths[normalized.group_code == group].sum(), 0.5):
                raise ValueError("Each group must span half the normalized axis")
        if not np.all(np.diff(normalized.position) == (1 if direction == "left" else -1)):
            raise ValueError("Curve or color strip has incorrect direction")
        row = pd.Series({"direction_best": "L2H" if direction == "left" else "H2L",
                         "gene_symbol": gene.gene_symbol, "tail_fdr": gene.fdr,
                         "edger_fdr": gene.edger_fdr})
        normalized_draw_curve(ax, normalized, peak, row, show_hd_text=True)
        ax.set_title(f"{gene.gene_symbol} ({direction.capitalize()}-tail enrichment)", fontsize=14)
        add_external_group_color_strip(ax, normalized, row.direction_best)
        distance = float(peak.sample_fraction_plot) if direction == "left" else 1 - float(peak.sample_fraction_plot)
        audit.append({"gene_id": gene.gene_id, "gene_symbol": gene.gene_symbol, "direction": direction,
                      "peak_position": int(peak_record.position), "peak_height": height,
                      "peak_height_percent": 100 * height,
                      "original_scan_x": original_x, "enrichment_score": pe,
                      "display_peak_x": float(peak.sample_fraction_plot),
                      "display_tail_distance": distance,
                      "display_tail_distance_percent": 100 * distance,
                      "AD_total_width": float(widths[normalized.group_code == 1].sum()),
                      "control_total_width": float(widths[normalized.group_code == -1].sum())})
        fig.legend(handles=[Patch(color=COLOR_AD, label="AD"),
                            Patch(color=COLOR_CONTROL, label="Control")],
                   loc="lower center", ncol=2, frameon=False, fontsize=12)
        for suffix in ("png", "svg", "pdf"):
            fig.savefig(args.output_dir / f"{gene.gene_symbol}.{suffix}",
                        dpi=300, bbox_inches="tight", pad_inches=0.04)
        plt.close(fig)
    pd.DataFrame(audit).to_csv(args.output_dir / "curve_plot_validation.tsv", sep="\t", index=False)
    print(f"Verified {len(audit)} scans for {len(genes)} genes; curves only")


if __name__ == "__main__":
    main()
