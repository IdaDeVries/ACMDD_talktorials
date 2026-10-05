"""Generate the schematic figures used in Talktorial T033.

The figures are deliberately drawn with matplotlib (rather than in a drawing
program) so that they can be regenerated and tweaked together with the notebook.

Usage
-----
    python scripts/make_figures.py

Writes the PNG files into ``../images``.
"""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

matplotlib.use("Agg")

HERE = Path(__file__).resolve().parent
IMAGES = HERE.parent / "images"
IMAGES.mkdir(exist_ok=True)

# Colour roles. The three categorical hues are the first three slots of a
# palette that was validated for colour-vision deficiency on all pairs.
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
TASK_COLORS = [BLUE, ORANGE, AQUA]

INK = "#0b0b0b"
INK_SOFT = "#52514e"
INK_MUTED = "#8a8880"
SURFACE = "#ffffff"
FILL_SOFT = "#f0efec"
GRID = "#d9d7d1"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "savefig.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
    }
)


# --------------------------------------------------------------------------- #
# Small drawing helpers                                                        #
# --------------------------------------------------------------------------- #
def blank_axes(fig, rect=(0, 0, 1, 1)):
    ax = fig.add_axes(rect)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_axis_off()
    return ax


def rounded_box(ax, x, y, w, h, *, label="", color=INK, face=SURFACE, lw=2.0,
                fontsize=11, fontweight="normal", text_color=None, zorder=3,
                label_dy=0.0, ls="-"):
    """Draw a rounded rectangle centred text box; (x, y) is the lower-left corner."""
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0,rounding_size=1.4",
        linewidth=lw,
        edgecolor=color,
        facecolor=face,
        linestyle=ls,
        zorder=zorder,
    )
    ax.add_patch(box)
    if label:
        ax.text(
            x + w / 2,
            y + h / 2 + label_dy,
            label,
            ha="center",
            va="center",
            fontsize=fontsize,
            fontweight=fontweight,
            color=text_color or color,
            zorder=zorder + 1,
        )
    return box


def arrow(ax, xy_from, xy_to, *, color=INK_SOFT, lw=1.8, ls="-", zorder=2,
          shrink=2.0, alpha=1.0):
    ax.add_patch(
        FancyArrowPatch(
            xy_from,
            xy_to,
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=lw,
            linestyle=ls,
            color=color,
            shrinkA=shrink,
            shrinkB=shrink,
            zorder=zorder,
            alpha=alpha,
        )
    )


def descriptor_vector(ax, x, y, n_cells, cell_w, cell_h, *, color=INK_SOFT,
                      lw=1.3, zorder=3, filled=None):
    """A row of small squares, echoing the descriptor vectors of Talktorial T032."""
    for i in range(n_cells):
        face = SURFACE if filled is None else (color if filled[i] else SURFACE)
        ax.add_patch(
            Rectangle(
                (x + i * cell_w, y),
                cell_w,
                cell_h,
                linewidth=lw,
                edgecolor=color,
                facecolor=face,
                zorder=zorder,
            )
        )


def molecule_glyph(ax, cx, cy, r, *, color=INK_SOFT, lw=1.6, zorder=3):
    """A tiny hexagon-with-substituents glyph standing in for 'a molecule'."""
    angles = np.linspace(0, 2 * np.pi, 7)[:-1] + np.pi / 6
    xs, ys = cx + r * np.cos(angles), cy + r * np.sin(angles)
    ax.plot(
        np.append(xs, xs[0]),
        np.append(ys, ys[0]),
        color=color,
        lw=lw,
        zorder=zorder,
        solid_joinstyle="round",
    )
    for idx, dxy in ((0, (1.7, 0.0)), (3, (-1.7, 0.0))):
        ax.plot(
            [xs[idx], xs[idx] + dxy[0] * r * 0.6],
            [ys[idx], ys[idx] + dxy[1] * r * 0.6],
            color=color,
            lw=lw,
            zorder=zorder,
        )
    ax.add_patch(Circle((xs[1], ys[1]), r * 0.22, color=color, zorder=zorder + 1))


def panel_title(ax, x, y, text, *, color=INK, fontsize=13, ha="center"):
    ax.text(x, y, text, ha=ha, va="center", fontsize=fontsize,
            fontweight="bold", color=color)


def save(fig, name):
    path = IMAGES / name
    fig.savefig(path, dpi=160, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(f"wrote {path}")


# --------------------------------------------------------------------------- #
# Figure 1 - from single-task QSPR to multitask modelling                      #
# --------------------------------------------------------------------------- #
def figure_concept():
    fig = plt.figure(figsize=(13.5, 6.2))
    ax = blank_axes(fig)

    task_names = ["Property A", "Property B", "Property C"]
    rows_y = [72, 50, 28]

    # ---------------- left: single-task ----------------
    panel_title(ax, 24, 94, "Single-task models: one model per property")
    for y, name, color in zip(rows_y, task_names, TASK_COLORS):
        molecule_glyph(ax, 4.5, y + 3.0, 2.4, color=color)
        descriptor_vector(ax, 8.5, y + 1.4, 7, 1.5, 3.2, color=color)
        arrow(ax, (20.0, y + 3.0), (25.0, y + 3.0), color=color)
        rounded_box(ax, 25.0, y - 1.0, 12.0, 8.0, label="model", color=color,
                    fontsize=10.5)
        arrow(ax, (37.0, y + 3.0), (41.5, y + 3.0), color=color)
        ax.text(42.5, y + 3.0, name, ha="left", va="center", fontsize=11,
                color=color, fontweight="bold")

    ax.text(
        24,
        12,
        "Each model sees only the compounds measured\n"
        "for its own property, and learns its own representation.",
        ha="center", va="center", fontsize=10.5, color=INK_SOFT,
    )

    # ---------------- divider ----------------
    ax.plot([54, 54], [6, 92], color=GRID, lw=1.6, ls=(0, (4, 4)), zorder=1)

    # ---------------- right: multitask ----------------
    panel_title(ax, 79, 94, "Multitask model: one model, several properties")

    for y, color in zip(rows_y, TASK_COLORS):
        molecule_glyph(ax, 59.0, y + 3.0, 2.4, color=color)
        descriptor_vector(ax, 63.0, y + 1.4, 7, 1.5, 3.2, color=color)

    hub_in, hub_out = 78.0, 92.5
    shared_x, shared_y, shared_w, shared_h = 80.5, 27.0, 10.0, 51.0
    rounded_box(ax, shared_x, shared_y, shared_w, shared_h, face=FILL_SOFT,
                color=INK_SOFT, lw=2.0)
    ax.text(shared_x + shared_w / 2, shared_y + shared_h / 2,
            "shared\nmodel", ha="center", va="center", fontsize=11.5,
            color=INK, fontweight="bold")

    # feature vectors funnel into a junction, then a single trunk enters the model
    for y, color in zip(rows_y, TASK_COLORS):
        ax.plot([74.0, hub_in], [y + 3.0, 52.5], color=color, lw=1.8,
                alpha=0.85, zorder=2, solid_capstyle="round")
    arrow(ax, (hub_in, 52.5), (shared_x, 52.5), color=INK_SOFT, lw=2.2, shrink=0)

    # one trunk leaves the model and fans out into the task-specific outputs
    ax.plot([shared_x + shared_w, hub_out], [52.5, 52.5], color=INK_SOFT,
            lw=2.2, zorder=2, solid_capstyle="round")
    for y, name, color in zip(rows_y, task_names, TASK_COLORS):
        arrow(ax, (hub_out, 52.5), (95.5, y + 3.0), color=color, alpha=0.9,
              shrink=0)
        ax.text(96.5, y + 3.0, name, ha="left", va="center", fontsize=11,
                color=color, fontweight="bold")

    ax.text(
        79,
        12,
        "One shared representation is fitted on all compounds at once;\n"
        "properties that are physically related can borrow strength from each other.",
        ha="center", va="center", fontsize=10.5, color=INK_SOFT,
    )

    save(fig, "multitask_concept.png")


# --------------------------------------------------------------------------- #
# Figure 2 - three ways of getting several outputs out of scikit-learn         #
# --------------------------------------------------------------------------- #
def figure_strategies():
    fig = plt.figure(figsize=(14.5, 5.8))
    ax = blank_axes(fig)

    def node(x, y, *, color=INK_SOFT, size=11):
        ax.plot([x], [y], marker="o", markersize=size, markerfacecolor=SURFACE,
                markeredgecolor=color, markeredgewidth=1.8, zorder=5)

    def leaf_vector(cx, cy, *, cell_w=2.2, cell_h=2.6):
        for j, color in enumerate(TASK_COLORS):
            ax.add_patch(
                Rectangle((cx - 1.5 * cell_w + j * cell_w, cy), cell_w, cell_h,
                          facecolor=color, edgecolor=SURFACE, lw=0.8, zorder=4)
            )

    # ------------- panel A: wrapper -------------
    panel_title(ax, 16, 95, "A \u00b7 Wrapper")
    ax.text(16, 88, "MultiOutputRegressor / MultiOutputClassifier",
            ha="center", va="center", fontsize=9.3, color=INK_SOFT,
            family="monospace")

    # the label matrix Y: rows = compounds, columns = tasks
    mat_x, mat_top, cell = 2.5, 74.0, 2.9
    n_rows = 6
    for i in range(n_rows):
        for j, color in enumerate(TASK_COLORS):
            ax.add_patch(
                Rectangle((mat_x + j * cell, mat_top - (i + 1) * cell), cell, cell,
                          facecolor=color, edgecolor=SURFACE, lw=1.0, zorder=3)
            )
    ax.text(mat_x + 1.5 * cell, mat_top + 2.0,
            "$Y$  (one column per task)", ha="center", va="bottom",
            fontsize=9.2, color=INK_SOFT)

    for j, color in enumerate(TASK_COLORS):
        y = 66 - j * 19
        col_x = mat_x + (j + 0.5) * cell
        ax.add_patch(
            FancyArrowPatch(
                (col_x, mat_top - n_rows * cell),
                (14.5, y + 5.5),
                connectionstyle="arc3,rad=-0.25",
                arrowstyle="-|>", mutation_scale=13, linewidth=1.6,
                color=color, shrinkA=3, shrinkB=2, zorder=2,
            )
        )
        rounded_box(ax, 14.5, y, 17.0, 11.0,
                    label=f"copy {j + 1}\nfit on column {j + 1}", color=color,
                    fontsize=9.0)

    ax.text(16, 10,
            "T independent clones of one estimator.\n"
            "Nothing is shared: identical to fitting T models by hand,\n"
            "but wrapped in a single scikit-learn object.",
            ha="center", va="center", fontsize=9.5, color=INK_SOFT)

    ax.plot([37, 37], [6, 92], color=GRID, lw=1.4, ls=(0, (4, 4)), zorder=1)

    # ------------- panel B: native multi-output tree -------------
    panel_title(ax, 56, 95, "B \u00b7 Native multi-output")
    ax.text(56, 88, "RandomForestRegressor, ExtraTrees, KNeighbors, Ridge, MLP",
            ha="center", va="center", fontsize=8.4, color=INK_SOFT,
            family="monospace")

    root = (56, 74)
    kids = [(48, 58), (64, 58)]
    leaves = [(42, 38), (52, 38), (60, 38), (70, 38)]
    for parent, child in ((root, kids[0]), (root, kids[1]),
                          (kids[0], leaves[0]), (kids[0], leaves[1]),
                          (kids[1], leaves[2]), (kids[1], leaves[3])):
        ax.plot(*zip(parent, child), color=INK_SOFT, lw=1.6, zorder=2)
    for cx, cy in [root] + kids:
        node(cx, cy)
    ax.text(59.5, 76.0, "each split is scored on the\nimpurity summed over all tasks",
            ha="left", va="center", fontsize=8.8, color=INK_MUTED)
    for cx, cy in leaves:
        leaf_vector(cx, cy - 3.2)
    ax.text(56, 28, "every leaf stores a vector of T values",
            ha="center", va="center", fontsize=9.5, color=INK_SOFT)
    ax.text(56, 14,
            "One tree structure is shared by all tasks.\n"
            "Splits are compromises that serve the tasks jointly.",
            ha="center", va="center", fontsize=9.5, color=INK_SOFT)

    ax.plot([76, 76], [6, 92], color=GRID, lw=1.4, ls=(0, (4, 4)), zorder=1)

    # ------------- panel C: xgboost -------------
    panel_title(ax, 89, 95, "C \u00b7 Boosted")
    ax.text(89, 88, "XGBoost  multi_strategy=...", ha="center", va="center",
            fontsize=9.0, color=INK_SOFT, family="monospace")

    ax.text(89, 79, "'one_output_per_tree'", ha="center", va="center",
            fontsize=9.0, color=INK, family="monospace")
    for j, color in enumerate(TASK_COLORS):
        x = 82 + j * 7.0
        ax.plot([x, x - 2.2], [71.5, 66.5], color=color, lw=1.4, zorder=2)
        ax.plot([x, x + 2.2], [71.5, 66.5], color=color, lw=1.4, zorder=2)
        node(x, 72.3, color=color, size=8)
        for dx in (-3.2, 1.2):
            ax.add_patch(Rectangle((x + dx, 63.8), 2.0, 2.4, facecolor=color,
                                   edgecolor=SURFACE, lw=0.6, zorder=4))
    ax.text(89, 58, "T trees per boosting round\n(the default)", ha="center",
            va="center", fontsize=8.8, color=INK_MUTED)

    ax.text(89, 47, "'multi_output_tree'", ha="center", va="center",
            fontsize=9.0, color=INK, family="monospace")
    ax.plot([89, 84], [40.5, 34.5], color=INK_SOFT, lw=1.6, zorder=2)
    ax.plot([89, 94], [40.5, 34.5], color=INK_SOFT, lw=1.6, zorder=2)
    node(89, 41.3, size=10)
    for cx in (84, 94):
        leaf_vector(cx, 31.5, cell_w=2.0, cell_h=2.4)
    ax.text(89, 25, "one tree per round,\nvector-valued leaves", ha="center",
            va="center", fontsize=8.8, color=INK_MUTED)

    ax.text(89, 13,
            "The gradient statistics of all tasks\n"
            "are summed when a split is scored.",
            ha="center", va="center", fontsize=9.5, color=INK_SOFT)

    save(fig, "multitask_strategies.png")


# --------------------------------------------------------------------------- #
# Figure 3 - the sparse label matrix                                           #
# --------------------------------------------------------------------------- #
def figure_label_matrix():
    rng = np.random.default_rng(7)
    n_rows, n_tasks = 14, 3
    measured = rng.random((n_rows, n_tasks)) < 0.62
    measured[0] = True
    measured[5] = True
    measured[9] = True
    # guarantee at least one measurement per compound
    for i in range(n_rows):
        if not measured[i].any():
            measured[i, rng.integers(n_tasks)] = True

    complete = measured.all(axis=1)

    fig = plt.figure(figsize=(13.5, 5.6))
    ax = blank_axes(fig)

    cell_w, cell_h = 6.5, 4.4
    task_names = ["Property A", "Property B", "Property C"]

    def draw_matrix(x0, y0, keep_rows, *, dim_dropped=False, mask_missing=False):
        for i in range(n_rows):
            kept = keep_rows[i]
            for j in range(n_tasks):
                x = x0 + j * cell_w
                y = y0 - i * cell_h
                if measured[i, j]:
                    face = TASK_COLORS[j] if kept else FILL_SOFT
                    edge = TASK_COLORS[j] if kept else GRID
                    hatch = None
                else:
                    face = SURFACE
                    edge = GRID
                    hatch = "///"
                ax.add_patch(
                    Rectangle(
                        (x, y),
                        cell_w * 0.92,
                        cell_h * 0.86,
                        facecolor=face,
                        edgecolor=edge,
                        hatch=hatch,
                        linewidth=1.0,
                        alpha=0.35 if (dim_dropped and not kept) else 1.0,
                        zorder=3,
                    )
                )
                if mask_missing and not measured[i, j]:
                    ax.plot(
                        [x + 0.9, x + cell_w * 0.92 - 0.9],
                        [y + 0.7, y + cell_h * 0.86 - 0.7],
                        color=INK_MUTED, lw=1.0, zorder=4,
                    )
                    ax.plot(
                        [x + 0.9, x + cell_w * 0.92 - 0.9],
                        [y + cell_h * 0.86 - 0.7, y + 0.7],
                        color=INK_MUTED, lw=1.0, zorder=4,
                    )
        for j, name in enumerate(task_names):
            ax.text(x0 + j * cell_w + cell_w * 0.46, y0 + cell_h * 1.0, name,
                    ha="center", va="bottom", fontsize=8.6, rotation=30,
                    color=INK_SOFT)

    top = 72
    left_x = 6
    mid_x = 42
    right_x = 78

    panel_title(ax, left_x + 9, 92, "The label matrix is sparse")
    draw_matrix(left_x, top, np.ones(n_rows, dtype=bool))
    ax.text(left_x + 9, 5,
            "Hatched cells are compounds that were\nnever measured for that property.",
            ha="center", va="center", fontsize=9.5, color=INK_SOFT)

    panel_title(ax, mid_x + 9, 92, "Option 1 · complete cases")
    draw_matrix(mid_x, top, complete, dim_dropped=True)
    ax.text(mid_x + 9, 5,
            f"Keep only rows measured for every task\n"
            f"({complete.sum()} of {n_rows} rows here). Simple, but\n"
            "most of the data is thrown away.",
            ha="center", va="center", fontsize=9.5, color=INK_SOFT)

    panel_title(ax, right_x + 9, 92, "Option 2 · per-task masking")
    draw_matrix(right_x, top, np.ones(n_rows, dtype=bool), mask_missing=True)
    ax.text(right_x + 9, 5,
            "Keep every row; when scoring a task,\n"
            "use only the rows with a label for it.\n"
            "Requires a model that tolerates gaps.",
            ha="center", va="center", fontsize=9.5, color=INK_SOFT)

    save(fig, "sparse_label_matrix.png")


# --------------------------------------------------------------------------- #
# Figure 4 - splitting strategies                                              #
# --------------------------------------------------------------------------- #
def figure_splits():
    rng = np.random.default_rng(3)

    # a toy 2-D "chemical space" made of a handful of scaffold clusters
    centres = np.array([[0.22, 0.75], [0.52, 0.86], [0.78, 0.62],
                        [0.30, 0.30], [0.66, 0.22], [0.88, 0.88]])
    sizes = [26, 18, 24, 22, 20, 12]
    points, cluster_id = [], []
    for k, (c, n) in enumerate(zip(centres, sizes)):
        points.append(rng.normal(c, 0.055, size=(n, 2)))
        cluster_id += [k] * n
    points = np.vstack(points)
    cluster_id = np.asarray(cluster_id)
    n = len(points)

    # split definitions -----------------------------------------------------
    random_test = np.zeros(n, dtype=bool)
    random_test[rng.choice(n, size=int(0.2 * n), replace=False)] = True

    cluster_test = np.isin(cluster_id, [2, 5])

    # "scaffold": whole small clusters to test, mimicking Bemis-Murcko groups
    scaffold_test = np.isin(cluster_id, [1, 4])

    # stratified on a synthetic label: keep the label balance in both sets
    label = (points[:, 0] + points[:, 1] > 1.05)
    strat_test = np.zeros(n, dtype=bool)
    for value in (False, True):
        idx = np.flatnonzero(label == value)
        strat_test[rng.choice(idx, size=max(1, int(0.2 * len(idx))),
                              replace=False)] = True

    panels = [
        ("Random split", random_test,
         "Rows are assigned at random.\nOptimistic: near-duplicates land on both sides."),
        ("Scaffold split", scaffold_test,
         "Groups = Bemis–Murcko scaffolds.\nTests generalisation to unseen core structures."),
        ("Cluster split", cluster_test,
         "Groups = fingerprint similarity clusters.\nTests generalisation to unseen chemical series."),
        ("Stratified random split", strat_test,
         "Random, but the balance of value bins\n(or classes) is preserved on both sides."),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(15.0, 4.3))
    for ax, (title, test_mask, caption) in zip(axes, panels):
        ax.scatter(points[~test_mask, 0], points[~test_mask, 1], s=34,
                   facecolor=BLUE, edgecolor=SURFACE, linewidth=0.6,
                   label="training set", zorder=3)
        ax.scatter(points[test_mask, 0], points[test_mask, 1], s=40,
                   facecolor=ORANGE, edgecolor=SURFACE, linewidth=0.6,
                   marker="D", label="test set", zorder=4)
        ax.set_title(title, fontsize=12, fontweight="bold", color=INK, pad=9)
        ax.set_xlim(0.03, 1.02)
        ax.set_ylim(0.03, 1.02)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color(GRID)
        ax.set_xlabel(caption, fontsize=9.0, color=INK_SOFT, labelpad=8)

    axes[0].text(-0.06, 0.5, "chemical space\n(2 descriptor axes)", rotation=90,
                 transform=axes[0].transAxes, ha="center", va="center",
                 fontsize=9, color=INK_MUTED)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
               fontsize=10.5, bbox_to_anchor=(0.5, 1.07))
    fig.subplots_adjust(top=0.80, bottom=0.18, wspace=0.12)
    save(fig, "splitting_methods.png")


if __name__ == "__main__":
    figure_concept()
    figure_strategies()
    figure_label_matrix()
    figure_splits()
