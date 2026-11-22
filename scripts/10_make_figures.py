"""Stage 10. Manuscript figures (1-7, S1), Table 1 and the supplementary tables.

Reads results/numbers.json, results/event_studies.json, results/sensitivity.json, the ingredient table and a Natural Earth country layer;
writes results/figures/Figure_N.{png,pdf}, results/tables/Table_1.{csv,md}, results/tables/Table_S7.csv, results/tables/ingredient_sample.csv
and results/supplement/ (Figure_S1 and Supplementary_Tables.xlsx).
"""

import json

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Patch

from pesticide_evidence import geography
from pesticide_evidence.config import path
from pesticide_evidence.stats import ABSENT, AUTHORIZED, GROUPS as NORTHERN_GROUPS, analysis_sample, design, pct_change, ppml

TEXT, TEXT2, MUTED = "#0b0b0b", "#52514e", "#7a7973"
GROUPS = [AUTHORIZED, *NORTHERN_GROUPS]
GROUP_LABEL = {"authorized": "Authorized", "withdrawn: hazard/treaty": "Withdrawn", "lapsed: other reason": "Lapsed", "absent from Northern registers": "Absent"}
GROUP_COLOR = {"authorized": "#a7a6a1", "withdrawn: hazard/treaty": "#1197b3", "lapsed: other reason": "#0a7a42", "absent from Northern registers": "#7a44a8"}
NORTH_COLOR, SOUTH_COLOR = "#184f95", "#3987e5"
NO_DATA, NORTH_FILL = "#efeeeb", "#cfcec9"

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.titleweight": "bold",
    "axes.labelsize": 7.5,
    "axes.labelweight": "bold",
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.frameon": False,
    "legend.fontsize": 7,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "axes.grid": True,
    "grid.linewidth": 0.4,
    "grid.alpha": 0.15,
    "axes.axisbelow": True,
    "axes.edgecolor": TEXT2,
    "axes.labelcolor": TEXT,
    "xtick.color": TEXT2,
    "ytick.color": TEXT2,
    "text.color": TEXT,
}
DOT_SIZE, DOT_LW, ERR_LW = 16, 0.4, 1.0
COL2 = 180 / 25.4

OUT_FIG = path("results", "figures")
OUT_TAB = path("results", "tables")
OUT_SUP = path("results", "supplement")
NUMBERS = json.load(open(path("results", "numbers.json")))
EVENTS = json.load(open(path("results", "event_studies.json")))


def save(fig, name):
    """Save a figure as PNG and PDF in the figures folder and close it."""
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT_FIG / f"{name}.{ext}", facecolor="white")
    plt.close(fig)


def ratio_axis(ax, lo=0.04, hi=4):
    """Log-scaled ratio axis with ticks from 0.05 to 2 and a reference line at 1."""
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.set_xticks([0.05, 0.1, 0.25, 0.5, 1, 2])
    ax.set_xticklabels(["0.05", "0.1", "0.25", "0.5", "1", "2"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.axvline(1, color=TEXT2, lw=0.6, zorder=1)
    ax.grid(axis="y", visible=False)


def dot_ci(ax, y, est, lo, hi, color):
    """Point estimate with a horizontal confidence interval at height y."""
    ax.hlines(y, lo, hi, color=color, lw=ERR_LW, zorder=2)
    ax.scatter([est], [y], s=DOT_SIZE, facecolor=color, edgecolor=TEXT, linewidth=DOT_LW, zorder=3)


RATIO_KEYS = {"absent": ABSENT, "lapsed": "lapsed: other reason", "withdrawn": "withdrawn: hazard/treaty"}


def panel_letters(fig, axes, y=0.985):
    """Bold panel letters at the left edge of each axis."""
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    for ax, letter in zip(axes, "abcdef"):
        x0 = inv.transform(ax.get_tightbbox(fig.canvas.get_renderer()).p0)[0]
        fig.text(max(x0, 0.004), y, letter, fontsize=8, fontweight="bold", ha="left", va="top", color=TEXT)


def forest_rows(ax, ys, rows, groups=("absent", "lapsed", "withdrawn"), spread=0.22):
    """Dots and intervals for the absent, lapsed and withdrawn ratios of each row, offset around the row height."""
    offsets = np.linspace(spread, -spread, len(groups))
    for y, row in zip(ys, rows):
        for g, off in zip(groups, offsets):
            dot_ci(ax, y + off, *row[g], GROUP_COLOR[RATIO_KEYS[g]])


def group_handles(groups):
    """Legend markers for the named groups."""
    return [marker(GROUP_COLOR[g], GROUP_LABEL[g]) for g in groups]


def marker(color, label, face=None, edge=TEXT, lw=DOT_LW):
    """Round legend marker with a label."""
    return Line2D([], [], marker="o", ls="", markersize=5, markerfacecolor=face or color, markeredgecolor=edge, markeredgewidth=lw, label=label)


def figure1():
    """Figure 1: map of the share of each country's reported mass in pesticides absent from Northern registers."""
    by = NUMBERS["country_absent_share"]["by_country"]
    share = {v["iso3"]: v["share"] * 100 for v in by.values()}
    world = gpd.read_file(f"zip://{path('raw', 'geography', 'ne_110m_admin_0_countries.zip')}")
    world = world[world.CONTINENT != "Antarctica"].copy()
    world["iso3"] = world.ADM0_A3
    world["share"] = world.iso3.map(share)
    world["north"] = world.iso3.isin(geography.NORTH_ISO)
    world = world.to_crs("EPSG:8857")
    cmap = LinearSegmentedColormap.from_list("absent", ["#ece4f3", "#b79bd1", "#7a44a8", "#3b1f62"])
    vmax = 30
    fig, ax = plt.subplots(figsize=(COL2, 3.4))
    world[world.north].plot(ax=ax, color=NORTH_FILL, edgecolor="white", linewidth=0.3)
    world[~world.north & world.share.isna()].plot(ax=ax, color=NO_DATA, edgecolor="white", linewidth=0.3)
    world[world.share.notna()].plot(ax=ax, column="share", cmap=cmap, vmin=0, vmax=vmax, edgecolor="white", linewidth=0.3)
    ax.set_axis_off()
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, vmax))
    cax = ax.inset_axes([0.035, 0.27, 0.2, 0.025])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.set_ticks([0, 10, 20, 30])
    cb.ax.tick_params(length=2, labelsize=6.5)
    cax.set_title("National mass in pesticides absent\nfrom Northern registers (%)", fontsize=7, fontweight="normal", loc="left", pad=3)
    for i, (color, label) in enumerate(((NORTH_FILL, "Global North"), (NO_DATA, "No reported use data"))):
        ax.add_patch(plt.Rectangle((0.035 + 0.0, 0.12 - i * 0.06), 0.016, 0.03, transform=ax.transAxes, facecolor=color, edgecolor="none"))
        ax.text(0.06, 0.135 - i * 0.06, label, transform=ax.transAxes, fontsize=7, va="center", ha="left")
    arrow = ax.inset_axes([0.935, 0.04, 0.04, 0.17])
    arrow.set_xlim(-1, 1)
    arrow.set_ylim(-1.4, 3.0)
    arrow.axis("off")
    arrow.fill([0, 0.55, 0], [1.6, -0.3, -0.05], facecolor=TEXT, edgecolor=TEXT, linewidth=0.7)
    arrow.fill([0, -0.55, 0], [1.6, -0.3, -0.05], facecolor="white", edgecolor=TEXT, linewidth=0.7)
    arrow.text(0, 2.0, "N", ha="center", va="bottom", fontsize=8, fontweight="bold")
    fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)
    save(fig, "Figure_1")


def figure2(s):
    """Figure 2: papers against Southern applied mass with the fitted model lines, and share of papers relative to share of mass by group."""
    X, names = design(s)
    b, _ = ppml(s.en_zh_tox.values.astype(float), X, np.arange(len(s)))
    base = X.mean(axis=0)
    for g in ("withdrawn", "lapsed", "absent"):
        base[names.index(g)] = 0
    grid = np.logspace(np.log10(50), np.log10(s.South.max()), 100)
    lines = {}
    for label, absent in (("authorized", 0), ("absent", 1)):
        x = np.tile(base, (len(grid), 1))
        x[:, names.index("ln_south")] = np.log(grid)
        x[:, names.index("absent")] = absent
        lines[label] = np.exp(x @ b)
    fig, (bx, ax) = plt.subplots(1, 2, figsize=(COL2, 3.5), gridspec_kw={"width_ratios": [1, 1.6], "wspace": 0.4})
    for g in GROUPS:
        d = s[s.grp4 == g]
        top = g == ABSENT
        ax.scatter(
            d.South,
            d.en_zh_tox,
            s=DOT_SIZE * (1.5 if top else 1),
            color=GROUP_COLOR[g],
            edgecolor=TEXT,
            linewidth=DOT_LW,
            alpha=0.95 if top else 0.8,
            zorder=4 if top else 2,
        )
    ax.plot(grid, lines["authorized"], color=TEXT2, lw=1.2, zorder=3)
    ax.plot(grid, lines["absent"], color=GROUP_COLOR[ABSENT], lw=1.6, zorder=5)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1, linscale=0.35)
    ax.set_ylim(0, 6000)
    ax.set_yticks([0, 1, 10, 100, 1000])
    ax.set_yticklabels(["0", "1", "10", "100", "1,000"])
    ax.set_xlabel("Southern applied mass (t, 2016–2020)")
    ax.set_ylabel("Toxicity papers, 2005–2025")
    gs = NUMBERS["group_shares"]
    order = ["authorized", "withdrawn: hazard/treaty", "lapsed: other reason", ABSENT]
    for i, g in enumerate(order):
        bx.scatter([gs[g]["papers_per_mass"]], [i], s=DOT_SIZE * 2.2, color=TEXT2, edgecolor=TEXT, linewidth=DOT_LW, zorder=3)
        bx.hlines(i, 1, gs[g]["papers_per_mass"], color=TEXT2, lw=1.0, zorder=2)
    bx.axvline(1, color=TEXT2, lw=0.6)
    bx.set_xscale("log")
    bx.set_xlim(0.3, 2.6)
    bx.set_xticks([0.5, 1, 2])
    bx.set_xticklabels(["0.5", "1", "2"])
    bx.xaxis.set_minor_formatter(NullFormatter())
    bx.set_yticks(range(len(order)))
    bx.set_yticklabels([GROUP_LABEL[g] for g in order])
    bx.set_ylim(len(order) - 0.5, -0.5)
    bx.grid(axis="y", visible=False)
    bx.tick_params(axis="y", length=0)
    bx.set_xlabel("Share of papers relative to\nshare of Southern mass")
    handles = group_handles(GROUPS) + [
        Line2D([], [], color=TEXT2, lw=1.2, label="Model fit, authorized"),
        Line2D([], [], color=GROUP_COLOR[ABSENT], lw=1.6, label="Model fit, absent"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=6, bbox_to_anchor=(0.5, 0.0), columnspacing=1.4, handletextpad=0.3)
    fig.subplots_adjust(left=0.1, right=0.985, top=0.97, bottom=0.27)
    panel_letters(fig, [bx, ax])
    save(fig, "Figure_2")


def figure3():
    """Figure 3: observed and expected papers for each absent pesticide."""
    df = pd.DataFrame(NUMBERS["missing_papers"]["by_ingredient"]).T.sort_values("south_t", ascending=True)
    fig, ax = plt.subplots(figsize=(COL2 * 0.62, 6.4))
    y = np.arange(len(df))
    ax.hlines(y, df.observed, df.expected, color=GROUP_COLOR[ABSENT], lw=1.0, alpha=0.45, zorder=1)
    ax.scatter(df.expected, y, s=DOT_SIZE, facecolor="white", edgecolor=GROUP_COLOR[ABSENT], linewidth=0.9, zorder=3)
    ax.scatter(df.observed, y, s=DOT_SIZE, color=GROUP_COLOR[ABSENT], edgecolor=TEXT, linewidth=DOT_LW, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(df.index)
    ax.set_ylim(-0.7, len(df) - 0.3)
    ax.set_xscale("symlog", linthresh=1, linscale=0.35)
    ax.set_xlim(0, 3000)
    ax.set_xticks([0, 1, 10, 100, 1000])
    ax.set_xticklabels(["0", "1", "10", "100", "1,000"])
    ax.set_xlabel("Toxicity papers, 2005–2025")
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", length=0)
    handles = [marker(GROUP_COLOR[ABSENT], "Observed"), marker(GROUP_COLOR[ABSENT], "Expected if authorized", face="white", edge=GROUP_COLOR[ABSENT], lw=0.9)]
    fig.legend(handles=handles, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.0), columnspacing=1.6, handletextpad=0.3)
    fig.subplots_adjust(left=0.2, right=0.98, top=0.99, bottom=0.1)
    save(fig, "Figure_3")


def section_rows(ax, sections):
    """Draw titled sections of labelled rows on an axis; returns the row positions."""
    y = 0
    ticks, labels = [], []
    for title, rows, plot in sections:
        ax.text(-0.03, y, title, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7, fontweight="bold", color=TEXT)
        y += 1
        ys = []
        for lab, row in rows:
            ys.append(y)
            ticks.append(y)
            labels.append(lab)
            y += 1
        plot(ys, [r for _, r in rows])
        y += 0.7
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_ylim(y - 0.5, -0.8)
    ax.tick_params(axis="y", length=0)


def figure4():
    """Figure 4: ratios by outcome definition and period, and by missing register."""
    h, w, r = NUMBERS["headline"], NUMBERS["windows"], NUMBERS["registers"]
    joint = r["joint_EU_US_CA"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(COL2, 4.2), gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.72})

    def three(ys, rows):
        """Draw the three-group rows on the left panel."""
        return forest_rows(ax, ys, rows)

    section_rows(
        ax,
        [
            (
                "Outcome definition",
                [
                    ("Primary", h["combined_tox_pubmed_plus_zh"]),
                    ("PubMed only", h["pubmed_tox_only"]),
                    ("All topics", h["all_topics_combined"]),
                    ("CAS-registry search", h["cas_registry_pubmed"]),
                ],
                three,
            ),
            ("Period", [("2005–2009", w["2005-09"]), ("2010–2014", w["2010-14"]), ("2015–2019", w["2015-19"]), ("2020–2025", w["2020-25"])], three),
        ],
    )
    ratio_axis(ax)
    ax.set_xlabel("Ratio of papers vs authorized (95% CI)")

    def single(ys, rows):
        """Draw single-estimate rows on the right panel."""
        return [dot_ci(bx, y, *row, TEXT2) for y, row in zip(ys, rows)]

    section_rows(
        bx,
        [
            (
                "Separate models",
                [
                    ("Not in EU", r["not_EU"]["not_EU"]),
                    ("Not in US", r["not_US"]["not_US"]),
                    ("Not in Canada", r["not_CA"]["not_CA"]),
                    ("In none of the three", r["none_of_three"]["none3"]),
                ],
                single,
            ),
            ("Joint model", [("Not in EU", joint["not_EU"]), ("Not in US", joint["not_US"]), ("Not in Canada", joint["not_CA"])], single),
        ],
    )
    ratio_axis(bx)
    bx.set_xlabel("Ratio of papers vs registered (95% CI)")
    fig.legend(
        handles=group_handles(["withdrawn: hazard/treaty", "lapsed: other reason", ABSENT]),
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.3, 0.0),
        columnspacing=1.4,
        handletextpad=0.3,
    )
    fig.subplots_adjust(left=0.165, right=0.985, top=0.97, bottom=0.14)
    panel_letters(fig, [ax, bx])
    save(fig, "Figure_4")


def figure5():
    """Figure 5: absent-group ratios by funder class, and the share of papers with a listed funder by group."""
    f = NUMBERS["funders_openalex"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(COL2, 3.0), gridspec_kw={"width_ratios": [1.05, 1], "wspace": 0.55})
    items = [
        ("All works", f["all"]),
        ("Northern funders", f["northern"]),
        ("Southern funders", f["southern"]),
        ("Industry funders", f["industry"]),
        ("Unclassified funders", f["unclassified"]),
    ]
    for i, (lab, (est, lo, hi)) in enumerate(items):
        dot_ci(ax, i, est, lo, hi, GROUP_COLOR[ABSENT])
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([i[0] for i in items])
    ax.set_ylim(len(items) - 0.5, -0.5)
    ax.tick_params(axis="y", length=0)
    ratio_axis(ax, 0.04, 2)
    ax.set_xlabel("Ratio of papers, absent vs authorized (95% CI)")
    north, south = f["share_works_northern_funder"], f["share_works_southern_funder"]
    wd = 0.38
    for i, g in enumerate(GROUPS):
        bx.bar(i - wd / 2, north[g] * 100, wd, color=NORTH_COLOR, edgecolor=TEXT, linewidth=DOT_LW)
        bx.bar(i + wd / 2, south[g] * 100, wd, color=SOUTH_COLOR, edgecolor=TEXT, linewidth=DOT_LW)
    bx.set_xticks(range(len(GROUPS)))
    bx.set_xticklabels([GROUP_LABEL[g] for g in GROUPS])
    bx.set_ylabel("Papers with a listed funder (%)")
    bx.set_ylim(0, 32)
    bx.grid(axis="x", visible=False)
    bx.legend(
        handles=[
            Patch(facecolor=NORTH_COLOR, edgecolor=TEXT, linewidth=0.5, label="Northern funder"),
            Patch(facecolor=SOUTH_COLOR, edgecolor=TEXT, linewidth=0.5, label="Southern funder"),
        ],
        loc="upper left",
    )
    fig.subplots_adjust(left=0.17, right=0.985, top=0.97, bottom=0.18)
    panel_letters(fig, [ax, bx])
    save(fig, "Figure_5")


def figure6():
    """Figure 6: absent-group ratios for human-relevant evidence types, and the share of pesticides with no papers by group."""
    h = NUMBERS["human"]
    prim = NUMBERS["headline"]["combined_tox_pubmed_plus_zh"]
    rows = [
        ("All toxicity (primary)", prim),
        ("Human studies", h["human"]),
        ("Epidemiology and\nbiomonitoring", h["epi"]),
        ("Clinical poisoning", h["poison"]),
        ("Rodent toxicity", h["mamm"]),
    ]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(COL2, 3.2), gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.75})
    ys = list(range(len(rows)))
    forest_rows(ax, ys, [r for _, r in rows])
    ax.set_yticks(ys)
    ax.set_yticklabels([label for label, _ in rows])
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.tick_params(axis="y", length=0)
    ratio_axis(ax, 0.04, 3.5)
    ax.set_xlabel("Ratio of papers vs authorized (95% CI)")
    zero = h["zero_paper_share_by_group"]
    keys = [("human", "Human studies"), ("epi", "Epidemiology"), ("poison", "Clinical poisoning"), ("mamm", "Rodent toxicity")]
    for i, (k, lab) in enumerate(keys):
        for g, off in zip(GROUPS, np.linspace(-0.27, 0.27, 4)):
            bx.scatter([zero[k][g] * 100], [i + off], s=DOT_SIZE, color=GROUP_COLOR[g], edgecolor=TEXT, linewidth=DOT_LW, zorder=3)
    bx.set_yticks(range(len(keys)))
    bx.set_yticklabels([k[1] for k in keys])
    bx.set_ylim(len(keys) - 0.5, -0.5)
    bx.set_xlim(-2, 70)
    bx.set_xlabel("Pesticides with zero papers (%)")
    bx.grid(axis="y", visible=False)
    bx.tick_params(axis="y", length=0)
    fig.legend(handles=group_handles(GROUPS), loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.0), columnspacing=1.4, handletextpad=0.3)
    fig.subplots_adjust(left=0.17, right=0.985, top=0.97, bottom=0.24)
    panel_letters(fig, [ax, bx])
    save(fig, "Figure_6")


def figure7():
    """Figure 7: change in Northern-authored and Southern-authored papers after Northern exits, by design and by years from the exit."""
    st, rv, et = EVENTS["stacked"], EVENTS["review_timing"], EVENTS["event_time"]

    rows = [
        ("Exit from any register,\nstill used ≥200 t in the South", st["any exit, still used >=200 t"]),
        ("Exit from any register,\nnot treaty-listed", st["not treaty-listed, any use"]),
        ("EU non-renewal at\nscheduled review", None),
    ]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(COL2, 3.3), gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.5})
    for i, (lab, r) in enumerate(rows):
        for who, off, col in (("North", -0.17, NORTH_COLOR), ("South", 0.17, SOUTH_COLOR)):
            y = i + off
            if r is not None:
                ax.hlines(y, *r[f"{who}_cal_ci"], color=col, lw=0.6, ls=(0, (3, 1.5)), zorder=1)
                dot_ci(ax, y, r[who], *r[f"{who}_ci"], col)
            else:
                b, se = rv["real"]["bN" if who == "North" else "bS"], rv["real"]["seN" if who == "North" else "seS"]
                dot_ci(ax, y, pct_change(b), pct_change(b - 1.96 * se), pct_change(b + 1.96 * se), col)
    ax.axvline(0, color=TEXT2, lw=0.6)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([x[0] for x in rows])
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlim(-60, 65)
    ax.set_xlabel("Change in papers after Northern exit (%)")
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", length=0)
    centers = {"-5 to -4": -4.5, "-3 to -2": -2.5, "0": 0, "1 to 2": 1.5, "3 to 5": 4}
    for who, off, col in (("North", -0.12, NORTH_COLOR), ("South", 0.12, SOUTH_COLOR)):
        xs = np.array([centers[k] + off for k in et[who]])
        mid = np.array([v["pct"] for v in et[who].values()])
        bx.plot(xs, mid, color=col, lw=1.0, zorder=2)
        bx.vlines(xs, [v["lo"] for v in et[who].values()], [v["hi"] for v in et[who].values()], color=col, lw=ERR_LW, zorder=2)
        bx.scatter(xs, mid, s=DOT_SIZE, color=col, edgecolor=TEXT, linewidth=DOT_LW, zorder=3)
    bx.axhline(0, color=TEXT2, lw=0.6)
    bx.axvline(-1, color=TEXT2, lw=0.6, ls=(0, (3, 2)))
    bx.set_xticks([-4.5, -2.5, 0, 1.5, 4])
    bx.set_xticklabels(["−5 to −4", "−3 to −2", "0", "1 to 2", "3 to 5"])
    bx.set_xlabel("Years from Northern exit")
    bx.set_ylabel("Change in papers (%)")
    bx.grid(axis="x", visible=False)
    fig.legend(
        handles=[marker(NORTH_COLOR, "Northern authors"), marker(SOUTH_COLOR, "Southern authors")],
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, 0.0),
        columnspacing=1.6,
        handletextpad=0.3,
    )
    fig.subplots_adjust(left=0.2, right=0.985, top=0.97, bottom=0.26)
    panel_letters(fig, [ax, bx])
    save(fig, "Figure_7")


def figure_s1():
    """Figure S1: absent pesticides by crop, in share of mass and papers per kilotonne."""
    crops = pd.DataFrame(NUMBERS["crops"]).T.astype(float).sort_values("absent_share", ascending=False)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(COL2, 3.6), gridspec_kw={"width_ratios": [1, 1], "wspace": 0.45})
    y = np.arange(len(crops))
    ax.barh(y, crops.absent_share * 100, color=GROUP_COLOR[ABSENT], height=0.65)
    ax.set_yticks(y)
    ax.set_yticklabels(crops.index)
    ax.set_ylim(len(crops) - 0.5, -0.5)
    ax.set_xlabel("Share of crop mass in absent pesticides (%)")
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", length=0)
    bx.hlines(y, crops.ppk_absent, crops.ppk_authorized, color=TEXT2, lw=1.0, zorder=1)
    bx.scatter(crops.ppk_authorized, y, s=DOT_SIZE * 1.5, color=GROUP_COLOR["authorized"], edgecolor=TEXT, linewidth=DOT_LW, zorder=3)
    bx.scatter(crops.ppk_absent, y, s=DOT_SIZE * 1.5, color=GROUP_COLOR[ABSENT], edgecolor=TEXT, linewidth=DOT_LW, zorder=4)
    bx.set_yticks(y)
    bx.set_yticklabels([])
    bx.set_ylim(len(crops) - 0.5, -0.5)
    bx.set_xlabel("Toxicity papers per kt applied")
    bx.grid(axis="y", visible=False)
    bx.tick_params(axis="y", length=0)
    fig.legend(
        handles=[marker(GROUP_COLOR["authorized"], "Authorized"), marker(GROUP_COLOR[ABSENT], "Absent")],
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, 0.0),
        columnspacing=1.6,
        handletextpad=0.3,
    )
    fig.subplots_adjust(left=0.12, right=0.985, top=0.97, bottom=0.22)
    panel_letters(fig, [ax, bx])
    OUT_SUP.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT_SUP / f"Figure_S1.{ext}", facecolor="white")
    plt.close(fig)


def supplementary_tables(s):
    """Write Supplementary_Tables.xlsx from the results files."""
    sens = json.load(open(path("results", "sensitivity.json")))["sensitivity"]

    def ci(v):
        """Ratio and interval as text."""
        return f"{v[0]:.2f} ({v[1]:.2f}–{v[2]:.2f})"

    t1 = s.assign(group=s.grp4.map(GROUP_LABEL))[["group", "eu_ok", "us_ok", "ca_ok", "South", "North", "age", "target", "pm_tox", "oa_zh_tox", "en_zh_tox"]].rename(
        columns={
            "eu_ok": "EU register",
            "us_ok": "US register",
            "ca_ok": "Canadian register",
            "South": "Southern use (t)",
            "North": "Northern use (t)",
            "age": "Market age (years)",
            "target": "Target class",
            "pm_tox": "PubMed toxicity papers",
            "oa_zh_tox": "Chinese-language toxicity papers",
            "en_zh_tox": "Toxicity papers (primary)",
        }
    )
    t1 = t1.sort_values(["group", "Southern use (t)"], ascending=[True, False]).rename_axis("Ingredient").reset_index()
    names = {
        "baseline": "Primary specification",
        ">=200 t Southern use": "At least 200 t of Southern use",
        "first PubMed record 1985 or later (age<=40)": "First PubMed record in 1985 or later",
        "first-filter corpus counts": "Counts from the first-filter corpus",
        "selected corpus counts": "Counts from the selected corpus",
        "market age of at least 20 years": "Market age of at least 20 years",
        "market age of at least 30 years": "Market age of at least 30 years",
        "+ China share of Southern use": "Adjusted for the Chinese share of Southern use",
    }
    rows = [
        {"Specification": names.get(k, k[0].upper() + k[1:]), "Pesticides": v["n"], "Absent": ci(v["absent"]), "Lapsed": ci(v["lapsed"]), "Withdrawn": ci(v["withdrawn"])}
        for k, v in sens.items()
        if "absent" in v
    ]
    loo, perm = sens["leave-one-out (absent), range of ratio"], sens["permutation within target class (400)"]
    rows.append(
        {
            "Specification": "Leave one absent pesticide out (range of the absent ratio)",
            "Pesticides": 294,
            "Absent": f"{loo['min']:.2f}–{loo['max']:.2f}",
            "Lapsed": "",
            "Withdrawn": "",
        }
    )
    rows.append(
        {
            "Specification": "Group labels permuted within target class, 400 draws (median ratio; draws as low as observed)",
            "Pesticides": 295,
            "Absent": f"{perm['median_permuted']:.2f}; {perm['reached_observed']}",
            "Lapsed": "",
            "Withdrawn": "",
        }
    )
    t2 = pd.DataFrame(rows)
    cov = NUMBERS["coverage"]["by_region"]
    t3a = pd.DataFrame(
        [
            {
                "Region": k,
                "Atlas (kt)": v["atlas_kt"],
                "FAOSTAT (kt)": v["fao_kt"],
                "GloPUT (kt)": v["glopu_kt"],
                "Atlas / FAOSTAT": v["atlas_over_fao"],
                "Atlas / GloPUT": v["atlas_over_glopu"],
            }
            for k, v in cov.items()
            if k != "None"
        ]
    )
    im = NUMBERS["atlas_imputation"]
    t3b = pd.DataFrame(
        [
            {"Item": "Countries with reported data", "Value": im["countries_reported"]},
            {"Item": "Countries with imputed data only", "Value": im["countries_imputed"]},
            {"Item": "Share of atlas mass in imputed countries (%)", "Value": round(im["imputed_share_of_atlas_mass"] * 100, 1)},
            {"Item": "Sub-Saharan African countries, reported", "Value": im["ssa_countries_reported"]},
            {"Item": "Sub-Saharan African countries, imputed only", "Value": im["ssa_countries_imputed"]},
            {"Item": "Absent ratio with imputed use added", "Value": ci(im["ratios_with_imputed_use"]["absent"])},
        ]
    )
    rc = NUMBERS["register_crosscheck"]
    t3c = pd.DataFrame([{"Register": k, "Agreement with independent table": round(v, 3)} for k, v in rc["agreement"].items()])
    ev = []
    for k, v in EVENTS["stacked"].items():
        ev.append(
            {
                "Design": {"any exit, still used >=200 t": "Exits, ≥200 t Southern use", "not treaty-listed, any use": "Exits, not treaty-listed"}.get(k, k),
                "Exit events": v["events"],
                "Pesticides": v["nT"],
                "Northern change (%)": f"{v['North']:.1f} ({v['North_ci'][0]:.1f} to {v['North_ci'][1]:.1f})",
                "Southern change (%)": f"{v['South']:.1f} ({v['South_ci'][0]:.1f} to {v['South_ci'][1]:.1f})",
                "Difference, z": round(v["diff_z"], 2),
                "Placebo or permutation runs": v["placebo_n"],
                "Placebo or permutation SD, Northern": round(v["placebo_North"]["sd"], 2),
                "One-sided p, Northern": round(v["p_North_le"], 3),
                "One-sided p, difference": round(v["p_diff_le"], 3),
                "Calibrated interval, Northern (%)": f"{v['North_cal_ci'][0]:.0f} to {v['North_cal_ci'][1]:.0f}",
                "Calibrated interval, Southern (%)": f"{v['South_cal_ci'][0]:.0f} to {v['South_cal_ci'][1]:.0f}",
            }
        )
    rt = EVENTS["review_timing"]
    ev.append(
        {
            "Design": "EU review timing",
            "Exit events": "",
            "Pesticides": f"{int(rt['real']['nT'])} treated, {int(rt['real']['nC'])} controls",
            "Northern change (%)": f"{rt['North_pct']:.1f}",
            "Southern change (%)": f"{rt['South_pct']:.1f}",
            "Difference, z": "",
            "Placebo or permutation runs": rt["perm_n"],
            "Placebo or permutation SD, Northern": round(rt["perm_sd_N"], 2),
            "One-sided p, Northern": "<0.004" if rt["p_N_le"] == 0 else round(rt["p_N_le"], 3),
            "One-sided p, difference": round(rt["p_diff_le"], 3),
            "Calibrated interval, Northern (%)": "",
            "Calibrated interval, Southern (%)": "",
        }
    )
    sf = NUMBERS["selection_flow"]
    labels = {
        "no_ingredient": "No named pesticide active ingredient in the title or abstract",
        "drug_only": "Only compounds whose literature is mostly drug use (warfarin, ivermectin, niclosamide, moxidectin, nicotine)",
        "single_mention": "Ingredient named once in the title and abstract",
        "publication_type": "Editorial, letter, comment, news item, biography, interview or erratum",
        "tool_use": "Pesticide used to induce a Parkinson model",
        "efficacy": "Efficacy or resistance wording with a target-pest word and no non-target wording",
        "pest_in_title": "Target-pest word in the title and no non-target wording",
        "residue": "Residue or analytical wording with no strong toxicity wording",
    }
    rows = [
        ("PubMed queries run, 2005–2025", f"{sf['queries']:,}"),
        ("Records retrieved", f"{sf['records']:,}"),
        (
            "Removed by the first rule filter: target-pest efficacy, synthesis, delivery or detection, disease-model tool, residue or analytical",
            f"{sf['excluded_by_rule']:,}",
        ),
        ("Removed by the first rule filter: no named pesticide or animal model, outside 2005–2025, or retracted", f"{sf['excluded_other']:,}"),
        ("Passing the first rule filter", f"{sf['rule_retained']:,}"),
    ]
    rows += [(f"   removed: {labels[k]}", f"{v:,}") for k, v in sf["removed_by_step"].items()]
    rows.append(("Selected (event-study corpus)", f"{sf['selected']:,}"))
    t5 = pd.DataFrame(rows, columns=["Step", "Records"])
    spec = json.load(open(path("results", "sensitivity.json")))["model_specification"]

    def est(v):
        """Ratio and interval from a record, as text."""
        return f"{v['ratio']:.2f} ({v['lo']:.2f}–{v['hi']:.2f})"

    cn = spec["coefficients"]
    rows5 = [
        ("Ratio per unit of ln Southern use (elasticity)", est(cn["ln_south"])),
        ("Ratio per unit of ln Northern use + 1", est(cn["ln_north"])),
        ("Ratio per unit of ln market age", est(cn["ln_age"])),
        ("Withdrawn vs authorized", est(cn["withdrawn"])),
        ("Lapsed vs authorized", est(cn["lapsed"])),
        ("Absent vs authorized", est(cn["absent"])),
        (f"Variance inflation factor, maximum ({spec['vif_max']['name'].replace('crop_', '').lower()} share)", f"{spec['vif_max']['value']:.1f}"),
        ("Variance inflation factor, median", f"{spec['vif_median']:.1f}"),
    ]
    rows5 += [(f"Absent vs authorized, {k}", f"{v['absent'][0]:.2f} ({v['absent'][1]:.2f}–{v['absent'][2]:.2f})") for k, v in spec["alternative_forms"].items()]
    t6 = pd.DataFrame(rows5, columns=["Item", "Estimate (95% CI)"])
    bal = EVENTS["balance"]
    names6 = {
        "papers_per_year_before": "Papers per year in the five years before the exit (matched)",
        "first_record_year": "Year of first PubMed record (matched)",
        "change_in_papers_before": "Change in ln papers, years −5 to −4 to years −2 to −1",
        "ln_south_use": "ln Southern use (t)",
        "ln_north_use": "ln Northern use + 1 (t)",
        "insecticide": "Share of insecticides",
        "herbicide": "Share of herbicides",
        "fungicide": "Share of fungicides",
    }
    rows6 = [(names6[k], f"{bal[k]['treated']:.2f}", f"{bal[k]['control']:.2f}", f"{bal[k]['std_diff']:.2f}") for k in names6]
    rows6.append(("Units (pesticide by exit year)", str(bal["n_treated_units"]), str(bal["n_control_units"]), ""))
    t7 = pd.DataFrame(rows6, columns=["Variable", "Treated", "Control", "Standardized difference"])
    t4 = pd.DataFrame(ev)
    OUT_TAB.mkdir(parents=True, exist_ok=True)
    t1.to_csv(OUT_TAB / "ingredient_sample.csv", index=False)
    t_s7 = table_s7(s)
    t_s7.to_csv(OUT_TAB / "Table_S7.csv", index=False)
    OUT_SUP.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(OUT_SUP / "Supplementary_Tables.xlsx") as xl:
        t2.to_excel(xl, sheet_name="Table S1", index=False)
        t3a.to_excel(xl, sheet_name="Table S2", index=False, startrow=0)
        t3b.to_excel(xl, sheet_name="Table S2", index=False, startrow=len(t3a) + 3)
        t3c.to_excel(xl, sheet_name="Table S2", index=False, startrow=len(t3a) + len(t3b) + 6)
        t5.to_excel(xl, sheet_name="Table S4", index=False)
        t6.to_excel(xl, sheet_name="Table S5", index=False)
        t7.to_excel(xl, sheet_name="Table S6", index=False)
        t4.set_index("Design").T.rename_axis("Measure").reset_index().to_excel(xl, sheet_name="Table S3", index=False)
        t_s7.to_excel(xl, sheet_name="Table S7", index=False)


def table_s7(s):
    """Table S7: the sample pesticides with CAS number, target class, group, register status and basis of removal."""

    def basis(r):
        """Reasons a withdrawn pesticide left the registers, as text."""
        if r.grp4 != NORTHERN_GROUPS[0]:
            return ""
        return "; ".join(label for label, flag in (("EU hazard", r.eu_hazard), ("US regulatory", r.us_reg), ("Treaty", r.treaty)) if flag)

    t = pd.DataFrame(
        {
            "Pesticide": s.index,
            "CAS RN": s.cas.values,
            "Target class": s.target.values,
            "Group": s.grp4.map(GROUP_LABEL).values,
            "EU": s.eu.values,
            "US": s.us.fillna(0).astype(int).values,
            "Canada": s.ca.fillna(0).astype(int).values,
            "Basis of removal": [basis(r) for r in s.itertuples()],
        }
    )
    return t.sort_values("Pesticide", key=lambda c: c.str.lower()).reset_index(drop=True)


def table1(s):
    """Table 1: sample by authorization group."""
    rows = []
    for g in GROUPS:
        d = s[s.grp4 == g]
        rows.append(
            {
                "Group": GROUP_LABEL[g],
                "Pesticides": len(d),
                "Median Southern use (t)": round(d.South.median()),
                "Median toxicity papers": d.en_zh_tox.median(),
                "Pesticides with zero papers (%)": round((d.en_zh_tox == 0).mean() * 100, 1),
                "In EU register (%)": round(d.eu_ok.mean() * 100, 1),
                "In US register (%)": round(d.us_ok.mean() * 100, 1),
                "In Canadian register (%)": round(d.ca_ok.mean() * 100, 1),
            }
        )
    t = pd.DataFrame(rows)
    OUT_TAB.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT_TAB / "Table_1.csv", index=False)
    (OUT_TAB / "Table_1.md").write_text(t.to_markdown(index=False) + "\n")
    return t


def main():
    """Check the numbers file against the sample, then write every figure and table."""
    plt.rcParams.update(RC)
    s = analysis_sample(pd.read_parquet(path("processed", "ingredient_table.parquet")))
    assert len(s) == NUMBERS["sample"]["n"]
    medians = s.groupby("grp4").en_zh_tox.median()
    assert all(abs(medians[g] - NUMBERS["raw_median_papers"][g]) < 1e-9 for g in GROUPS)
    figure1()
    figure2(s)
    figure3()
    figure4()
    figure5()
    figure6()
    figure7()
    figure_s1()
    table1(s)
    supplementary_tables(s)


if __name__ == "__main__":
    main()
