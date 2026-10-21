"""Curve-only adaptation of source_final_real_gene_colorstrip.py.

Retains its group-normalized coordinates, independent color strip, colors,
peak-height and tail-distance arrows. Display distance is not statistic x.
"""
import matplotlib
matplotlib.use("Agg")
from matplotlib.colors import ListedColormap
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

STRIP_HEIGHT_SCALE = 2.0
COLOR_AD = "#d62828"
COLOR_CONTROL = "#1f9e4d"
STRIP_CMAP = ListedColormap([COLOR_CONTROL, COLOR_AD])
ARROW_CAP_HALF_HEIGHT = 0.035
H_CAP_HALF_WIDTH = 0.018
D_ARROW_Y_RATIO = 0.05
LEFT_LOWER_BLOCK_FRACTION = 0.22
STRIP_HEIGHT_FRACTION = 0.04


def normalize_curve_x(curve: pd.DataFrame) -> pd.DataFrame:
    c = curve.copy()
    c = c.sort_values("sample_fraction_plot").reset_index(drop=True)
    group = c["group_code"].to_numpy(dtype=int)
    n_ad = int(np.sum(group == 1))
    n_ctrl = int(np.sum(group == -1))
    width = np.where(group == 1, 0.5 / n_ad, 0.5 / n_ctrl)
    cum = np.cumsum(width)
    left = cum - width
    right = cum
    c["sample_fraction_plot"] = np.clip(left + width / 2.0, 0.0, 1.0)
    c["cell_left"] = np.clip(left, 0.0, 1.0)
    c["cell_right"] = np.clip(right, 0.0, 1.0)
    c = c.sort_values("sample_fraction_plot").reset_index(drop=True)
    return c


def _strip_cell_values(group_code) -> np.ndarray:
    arr = np.asarray(group_code, dtype=int)
    return np.where(arr == 1, 1, 0).astype(float).reshape(1, -1)


def format_p(value: float) -> str:
    value = float(value)
    if not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 1e-3 else f"{value:.3f}"


def fdr_only_label(row: pd.Series) -> str:
    return "\n".join([
        f"TailEnrich: FDR={format_p(row['tail_fdr'])}",
        f"edgeR: FDR={format_p(row['edger_fdr'])}",
    ])


def add_external_group_color_strip(ax, curve: pd.DataFrame, direction: str):
    curve = curve.sort_values("sample_fraction_plot")
    group = np.asarray(curve["group_code"], dtype=int)
    if "cell_left" in curve.columns and "cell_right" in curve.columns:
        edges = np.concatenate([
            np.asarray(curve["cell_left"], dtype=float),
            [float(curve["cell_right"].iloc[-1])],
        ])
        edges = np.clip(edges, 0.0, 1.0)
    else:
        n = len(curve)
        edges = np.linspace(0.0, 1.0, n + 1)
    full_box = ax.get_position()
    strip_height = full_box.height * STRIP_HEIGHT_FRACTION * STRIP_HEIGHT_SCALE
    lower_block_fraction = (
        LEFT_LOWER_BLOCK_FRACTION
        + STRIP_HEIGHT_FRACTION * (STRIP_HEIGHT_SCALE - 1.0)
    )
    lower_block = full_box.height * lower_block_fraction
    ax.set_position([
        full_box.x0, full_box.y0 + lower_block,
        full_box.width, full_box.height - lower_block,
    ])
    strip_ax = ax.figure.add_axes([
        full_box.x0, full_box.y0, full_box.width, strip_height,
    ])
    strip_ax.pcolormesh(
        edges, np.array([0.0, 1.0]), _strip_cell_values(group),
        cmap=STRIP_CMAP, shading="flat", edgecolors="none", linewidth=0,
    )
    strip_ax.set_xlim(0.0, 1.0)
    strip_ax.set_ylim(0.0, 1.0)
    strip_ax.set_xticks([])
    strip_ax.set_yticks([])
    for spine in strip_ax.spines.values():
        spine.set_visible(False)
    return strip_ax


def normalized_draw_curve(ax, curve, peak, row, show_hd_text=False):
    direction = str(row["direction_best"])
    line_color = COLOR_CONTROL if direction == "L2H" else COLOR_AD
    peak_x = float(peak["sample_fraction_plot"])
    peak_y = float(peak["running_score"])
    if abs(peak_y) < 1e-12:
        peak_y = 0.0
    score_min = float(curve["running_score"].min())
    score_max = float(curve["running_score"].max())
    score_span = max(score_max - score_min, abs(score_max) * 0.25, 1e-6)
    x_mark_start, x_mark_end = (peak_x, 1.0) if direction == "H2L" else (0.0, peak_x)
    d_arrow_y = -D_ARROW_Y_RATIO * score_span
    x_cap_halfheight = ARROW_CAP_HALF_HEIGHT * score_span
    y_min = min(score_min - 0.06 * score_span, d_arrow_y - 4.0 * x_cap_halfheight)
    # Reserve text space for opposite-direction panels with no positive peak.
    y_top = max(score_max + 0.10 * score_span, 0.30 * score_span)
    ax.axhline(0, color="#707070", linewidth=0.9, zorder=1)
    cx = np.concatenate([[0.0], curve["sample_fraction_plot"].to_numpy(), [1.0]])
    cy = np.concatenate([[0.0], curve["running_score"].to_numpy(), [0.0]])
    ax.plot(cx, cy, color=line_color, linewidth=1.6, zorder=2)
    ax.plot([peak_x, peak_x], [0, peak_y], color=line_color,
            linewidth=0.9, linestyle="--", zorder=2)
    ax.scatter([peak_x], [peak_y], color=line_color, s=20, zorder=4)
    ax.annotate("", xy=(peak_x, peak_y), xytext=(peak_x, 0),
                arrowprops={"arrowstyle": "<->", "color": line_color,
                            "linewidth": 0.9, "mutation_scale": 9}, zorder=3)
    for height in (0, peak_y):
        ax.plot([peak_x - H_CAP_HALF_WIDTH, peak_x + H_CAP_HALF_WIDTH], [height, height],
                color=line_color, linewidth=1.0)
    ax.annotate("", xy=(x_mark_end, d_arrow_y), xytext=(x_mark_start, d_arrow_y),
                arrowprops={"arrowstyle": "<->", "color": line_color,
                            "linewidth": 0.9, "mutation_scale": 9}, zorder=3)
    for xv in (x_mark_start, x_mark_end):
        ax.plot([xv, xv], [d_arrow_y - x_cap_halfheight, d_arrow_y + x_cap_halfheight],
                color=line_color, linewidth=1.0)
    if show_hd_text:
        label_x = max(peak_x - 0.022, 0.05) if direction == "H2L" else min(peak_x + 0.022, 0.95)
        if peak_y > 0:
            ax.text(label_x, 0.22 * peak_y, r"$h_{\max}$", color=line_color,
                    fontsize=11, ha="right" if direction == "H2L" else "left", va="center")
        else:
            ax.text(0.97 if direction == "L2H" else 0.03, 0.70,
                    r"$h_{\max}=0$", transform=ax.transAxes,
                    ha="right" if direction == "L2H" else "left", va="top",
                    color=line_color, fontsize=11)
        ax.text((x_mark_start + x_mark_end) / 2, d_arrow_y - 1.25 * x_cap_halfheight,
                "d", color=line_color, fontsize=11, ha="center", va="top")
    p_on_left = direction == "H2L"
    ax.text(0.02 if p_on_left else 0.98, 0.94, fdr_only_label(row), color="#333333",
            fontsize=11, linespacing=1.0,
            ha="left" if p_on_left else "right", va="top", transform=ax.transAxes)
    ax.set_xlim(0, 1)
    ax.set_ylim(y_min, y_top)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["100", "75", "50", "25", "0"] if direction == "H2L"
                       else ["0", "25", "50", "75", "100"])
    ax.set_xlabel("All samples (%)", fontsize=12, color="#333333")
    ax.set_ylabel("Disease samples (%)", fontsize=12, color="#333333", labelpad=8)
    # Running scores remain fractions; display the axis in percentage units.
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0, decimals=0, symbol=""))
    ax.tick_params(axis="both", labelsize=12, colors="#333333")
    ax.grid(True, color="#e5e5e5", linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("#333333")
        spine.set_linewidth(0.7)
