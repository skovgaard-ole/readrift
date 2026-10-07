"""Summary figures.

These are new -- the PostScript version reported four numbers on its front
page and nothing else, so the shape of the data behind those numbers was
invisible.  Each figure answers one question:

1. classification -- what kinds of read are there, and how many?
2. coverage       -- is the reference evenly covered, and where are the holes?
3. read lengths   -- is ``--min-read-length`` throwing away useful data?
4. distances      -- is the short/long cut-off in a sensible place?
5. events         -- where are the structural events, all on one page?
6. identity       -- how well do the reads actually match? (``--identity``)
7. table          -- the same numbers in text, per contig.

Design notes (data-viz method):

* Form is chosen per question: magnitude -> bars, distribution -> histogram,
  position -> profile, location -> event track, exact values -> table.
* Colour is assigned last, by role: the three read classes take the reference
  palette's first three slots, inversions take the fixed ``critical`` status
  step because an inversion is a *state* of a divided read rather than a
  fourth class.  Status colour never appears without its label.
* Slot 3 (green) is below 3:1 on white, so the relief rule applies: every
  figure using it carries direct labels, and figure 7 is the table view.
* No dual axes anywhere; one measure per axis.
* The palette validator could not be run in this environment, so the colours
  are taken unchanged from the reference palette's documented pre-validated
  slots rather than mixed from ramps.  Re-run
  ``validate_palette.js "#1baf7a,#2a78d6,#eb6834" --mode light --surface "#ffffff"``
  if the palette is ever changed.
* There is no hover layer: this is a printed PDF, not an HTML chart.
"""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np
from matplotlib.backends.backend_pdf import FigureCanvasPdf
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from readrift.inputs.reference import Reference
from readrift.models import ReadClass
from readrift.params import Params
from readrift.render import theme
from readrift.stats import Stats

#: Most contigs to show in the per-contig figures before folding the rest into
#: a note.  Beyond this the panels stop being legible.
MAX_PANELS = 6
MAX_EVENT_ROWS = 10

_TITLE_SIZE = 21.0
_SUBTITLE_SIZE = 12.5
_AXIS_LABEL_SIZE = 13.0
_TICK_SIZE = 11.0
_ANNOTATION_SIZE = 11.5


# --------------------------------------------------------------------------
# Page scaffolding
# --------------------------------------------------------------------------


def _page() -> Figure:
    width, height = theme.A3_LANDSCAPE
    fig = Figure(figsize=(width / 72.0, height / 72.0), facecolor=theme.SURFACE)
    FigureCanvasPdf(fig)
    return fig


def _heading(fig: Figure, title: str, subtitle: str = "") -> None:
    fig.text(
        0.055,
        0.945,
        title,
        fontsize=_TITLE_SIZE,
        family=theme.SANS,
        color=theme.INK,
        ha="left",
        va="top",
    )
    if subtitle:
        fig.text(
            0.055,
            0.905,
            subtitle,
            fontsize=_SUBTITLE_SIZE,
            family=theme.SANS,
            color=theme.INK_SECONDARY,
            ha="left",
            va="top",
        )


def _style(ax, grid_axis: str = "y") -> None:
    """Recessive chrome: no box, hairline grid behind the data."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme.AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(
        colors=theme.INK_MUTED, labelsize=_TICK_SIZE, length=4, width=0.8
    )
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_family(theme.SANS)
    if grid_axis != "none":
        ax.grid(axis=grid_axis, color=theme.GRID, linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)


def _axis_label(ax, xlabel: str = "", ylabel: str = "") -> None:
    if xlabel:
        ax.set_xlabel(
            xlabel, fontsize=_AXIS_LABEL_SIZE, family=theme.SANS, color=theme.INK_SECONDARY
        )
    if ylabel:
        ax.set_ylabel(
            ylabel, fontsize=_AXIS_LABEL_SIZE, family=theme.SANS, color=theme.INK_SECONDARY
        )


def _empty(fig: Figure, message: str) -> Figure:
    fig.text(
        0.5,
        0.5,
        message,
        fontsize=16,
        family=theme.SANS,
        color=theme.INK_MUTED,
        ha="center",
        va="center",
    )
    return fig


def _scaled_axis(length: int) -> tuple[float, str]:
    """Divisor and unit label for a reference-coordinate axis."""
    if length >= 2_000_000:
        return 1_000_000.0, "Mb"
    if length >= 2_000:
        return 1_000.0, "kb"
    return 1.0, "bp"


#: Width of one character as a fraction of the font size, for the sans face
#: matplotlib falls back to.  Deliberately generous.
_CHAR_WIDTH = 0.58


def _label_room(labels: list[str], size: float, cap: float = 0.30) -> float:
    """Left margin, as a figure fraction, wide enough for *labels*.

    Every figure here places its axes with :meth:`Figure.add_axes` and an
    explicit rectangle, which is what keeps the pages consistent -- but it also
    means nothing widens the margin to fit a tick label.  A long category name
    simply runs off the left edge of the sheet, so the room it needs has to be
    reserved up front.
    """
    longest = max((len(text) for text in labels), default=0)
    room = (longest * size * _CHAR_WIDTH + 10.0) / theme.A3_LANDSCAPE[0]
    return min(room, cap)


# --------------------------------------------------------------------------
# 1. Classification
# --------------------------------------------------------------------------


def figure_classification(stats: Stats) -> Figure:
    fig = _page()
    _heading(
        fig,
        "How the reads map",
        "Every read that passed the length and match filters, by the shape of its "
        "alignment. Inversions are a property of divided reads, not a fourth class.",
    )
    if stats.total_reads == 0:
        return _empty(fig, "No reads were classified.")

    classes = list(theme.CLASS_ORDER)
    counts = [stats.counts[c] for c in classes]
    positions = [3.0, 2.0, 1.0]
    tick_labels = [theme.CLASS_LABELS[c] for c in classes] + [
        "Reads indicating an inversion"
    ]

    # The class names are the y tick labels and they are long -- "Divided reads,
    # short distance" is 29 characters, which at 13 pt runs a good 200 pt to the
    # left of the axes.  Reserve that instead of letting it run off the sheet.
    plot_right = 0.715          # the side panel starts at 0.78
    left = _label_room(tick_labels, _AXIS_LABEL_SIZE)
    ax = fig.add_axes((left, 0.18, plot_right - left, 0.62))
    _style(ax, grid_axis="x")

    ax.barh(
        positions,
        counts,
        height=0.52,
        color=[theme.CLASS_COLORS[c] for c in classes],
        zorder=3,
    )

    inversion_y = -0.4
    ax.barh(
        [inversion_y],
        [stats.inverted_reads],
        height=0.52,
        color=theme.INVERSION,
        zorder=3,
    )

    ax.set_yticks(positions + [inversion_y])
    ax.set_yticklabels(
        tick_labels,
        fontsize=_AXIS_LABEL_SIZE,
        family=theme.SANS,
    )
    for tick in ax.get_yticklabels():
        tick.set_color(theme.INK)

    largest = max(max(counts), stats.inverted_reads, 1)
    ax.set_xlim(0, largest * 1.22)
    ax.set_ylim(-1.1, 3.7)

    # Direct labels -- the relief rule for the low-contrast green, and the
    # reason no value axis is needed beyond a coarse grid.
    # strict: the three are one row each per class, and a mismatch would
    # silently drop a bar's label rather than fail.
    for y, value, share in [
        *zip(positions, counts, [stats.fraction(c) for c in classes], strict=True),
        (inversion_y, stats.inverted_reads, stats.inverted_fraction),
    ]:
        ax.text(
            value + largest * 0.015,
            y,
            f"{value:,}   {share:.1f}%",
            fontsize=_ANNOTATION_SIZE,
            family=theme.SANS,
            color=theme.INK,
            va="center",
            ha="left",
            zorder=4,
        )

    ax.axhline(0.25, color=theme.GRID, linewidth=1.0, zorder=1)
    ax.text(
        largest * 1.20,
        0.05,
        "subset of the divided reads above",
        fontsize=_ANNOTATION_SIZE - 1,
        family=theme.SANS,
        color=theme.INK_MUTED,
        ha="right",
        va="top",
    )

    _axis_label(ax, xlabel="Reads")

    side = fig.add_axes((0.78, 0.18, 0.17, 0.62))
    side.set_axis_off()
    lines = [
        ("Total reads mapped", f"{stats.total_reads:,}"),
        ("Mean coverage", f"{stats.mean_coverage:.2f}x"),
        ("Matched bases", f"{stats.matched_bases:,}"),
        ("Read N50", f"{stats.n50:,} bp"),
    ]
    for index, (name, value) in enumerate(lines):
        y = 0.92 - index * 0.19
        side.text(0, y, name, fontsize=_ANNOTATION_SIZE, family=theme.SANS,
                  color=theme.INK_SECONDARY, transform=side.transAxes)
        side.text(0, y - 0.075, value, fontsize=19, family=theme.SANS,
                  color=theme.INK, transform=side.transAxes)

    return fig


# --------------------------------------------------------------------------
# 2. Coverage profile
# --------------------------------------------------------------------------


def figure_coverage(stats: Stats, reference: Reference) -> Figure:
    fig = _page()
    shown = [c for c in reference.contigs if c.length > 0][:MAX_PANELS]
    hidden = max(0, len([c for c in reference.contigs if c.length > 0]) - len(shown))
    _heading(
        fig,
        "Coverage along the reference",
        "Aligned depth in "
        + (f"the {len(shown)} largest contigs. " if hidden else "each contig. ")
        + "Troughs are where reads stop matching -- the same places the ReadRift "
        "pages show divided reads. Every panel shares one depth scale, so they can be "
        "read against each other."
        + (f" {hidden} shorter contig(s) not shown." if hidden else ""),
    )
    if not shown or stats.total_reads == 0:
        return _empty(fig, "No coverage to plot.")

    top, bottom = 0.86, 0.09
    gap = 0.035
    height = (top - bottom - gap * (len(shown) - 1)) / len(shown)

    # One depth scale across every panel.  Scaling each panel to its own peak
    # made a 5x contig look as well covered as a 500x one, which is the opposite
    # of what this figure is for.
    peak = 1.0
    for contig in shown:
        cs = stats.per_contig[contig.name]
        contig_depth = cs.depth()
        peak = max(
            peak,
            float(contig_depth.max()) if contig_depth.size else 0.0,
            cs.mean_coverage,
        )
    depth_limit = peak * 1.25

    for index, contig in enumerate(shown):
        cs = stats.per_contig[contig.name]
        ax = fig.add_axes((0.055, top - (index + 1) * height - index * gap, 0.90, height))
        _style(ax)

        divisor, unit = _scaled_axis(contig.length)
        x = cs.positions() / divisor
        depth = cs.depth()

        ax.fill_between(x, depth, color=theme.CLASS_COLORS[ReadClass.SHORT_DIVIDED],
                        alpha=0.16, linewidth=0, zorder=2)
        ax.plot(x, depth, color=theme.CLASS_COLORS[ReadClass.SHORT_DIVIDED],
                linewidth=1.4, zorder=3)

        mean = cs.mean_coverage
        ax.axhline(mean, color=theme.INK_SECONDARY, linewidth=1.0,
                   linestyle=(0, (5, 4)), zorder=4)
        # Inside the panel, right-aligned above the line: anchored to the last
        # data point it ran off the right edge of the sheet.  `get_yaxis_transform`
        # is axes-fraction in x, data in y, which is exactly what this needs.
        ax.text(
            0.997,
            mean,
            f"mean {mean:.1f}x",
            transform=ax.get_yaxis_transform(),
            fontsize=_ANNOTATION_SIZE - 1,
            family=theme.SANS,
            color=theme.INK_SECONDARY,
            va="bottom",
            ha="right",
        )

        ax.set_xlim(0, contig.length / divisor)
        ax.set_ylim(0, depth_limit)
        ax.text(
            0.0,
            1.02,
            f"{contig.name}   {contig.length:,} bp   {cs.reads:,} reads",
            transform=ax.transAxes,
            fontsize=_ANNOTATION_SIZE,
            family=theme.SANS,
            color=theme.INK,
            va="bottom",
            ha="left",
        )
        _axis_label(ax, ylabel="depth")
        if index == len(shown) - 1:
            _axis_label(ax, xlabel=f"Position ({unit})")

    return fig


# --------------------------------------------------------------------------
# 3. Read lengths
# --------------------------------------------------------------------------


def figure_read_lengths(stats: Stats, params: Params) -> Figure:
    fig = _page()
    _heading(
        fig,
        "Read length distribution",
        "Divided reads should skew longer than undivided ones -- a longer read "
        "spans more junctions. If the distribution is cut off at the left, "
        "--min-read-length is discarding data.",
    )

    lengths = stats.all_read_lengths()
    if lengths.size == 0:
        return _empty(fig, "No reads to plot.")

    ax = fig.add_axes((0.055, 0.14, 0.90, 0.68))
    _style(ax)

    low = max(int(lengths.min()), 1)
    high = max(int(lengths.max()), low + 1)
    bins = np.logspace(np.log10(low), np.log10(high), 60)

    for cls in theme.CLASS_ORDER:
        data = np.asarray(stats.read_lengths[cls], dtype=np.int64)
        if data.size == 0:
            continue
        counts, edges = np.histogram(data, bins=bins)
        centres = np.sqrt(edges[:-1] * edges[1:])
        ax.step(
            centres,
            counts,
            where="mid",
            color=theme.CLASS_COLORS[cls],
            linewidth=2.0,
            label=theme.CLASS_LABELS[cls],
            zorder=3,
        )
        # Direct label at the mode of each curve, so identity never rests on
        # the legend alone.
        peak = int(np.argmax(counts))
        if counts[peak]:
            ax.text(
                centres[peak],
                counts[peak] * 1.06,
                theme.CLASS_LABELS[cls].split(",")[0],
                fontsize=_ANNOTATION_SIZE - 1,
                family=theme.SANS,
                color=theme.CLASS_COLORS[cls],
                ha="center",
                va="bottom",
                zorder=4,
            )

    ax.axvline(stats.n50, color=theme.INK_SECONDARY, linewidth=1.0,
               linestyle=(0, (5, 4)), zorder=5)
    ax.text(
        stats.n50,
        ax.get_ylim()[1] * 0.97,
        f" N50 = {stats.n50:,} bp",
        fontsize=_ANNOTATION_SIZE,
        family=theme.SANS,
        color=theme.INK_SECONDARY,
        ha="left",
        va="top",
    )

    if params.min_read_length > low:
        ax.axvline(params.min_read_length, color=theme.INK_MUTED, linewidth=1.0,
                   linestyle=(0, (2, 3)), zorder=5)
        ax.text(
            params.min_read_length,
            ax.get_ylim()[1] * 0.80,
            f" --min-read-length = {params.min_read_length:,}",
            fontsize=_ANNOTATION_SIZE - 1,
            family=theme.SANS,
            color=theme.INK_MUTED,
            ha="left",
            va="top",
        )

    ax.set_xscale("log")
    _axis_label(ax, xlabel="Read length (bp, log scale)", ylabel="Reads")
    legend = ax.legend(frameon=False, fontsize=_ANNOTATION_SIZE, loc="upper right")
    for text in legend.get_texts():
        text.set_color(theme.INK_SECONDARY)
        text.set_family(theme.SANS)

    return fig


# --------------------------------------------------------------------------
# 4. Division distance
# --------------------------------------------------------------------------


def figure_division_distance(stats: Stats, params: Params) -> Figure:
    fig = _page()
    # Clamped: the threshold is drawn on a log axis, where zero has no place.
    cut = max(params.division_cut, 1)
    _heading(
        fig,
        "Where the short/long cut-off falls",
        "For every divided read: how much more reference it spans than its own "
        "length accounts for. Two clear modes mean the cut-off separates real "
        "populations; one broad mode means it is arbitrary.",
    )

    distances = np.asarray(stats.distances, dtype=np.int64)
    if distances.size == 0:
        return _empty(fig, "No divided reads to plot.")

    ax = fig.add_axes((0.055, 0.14, 0.90, 0.68))
    _style(ax)

    clipped = np.clip(distances, 1, None)
    high = max(int(clipped.max()), 2)
    bins = np.logspace(0, np.log10(high), 60)

    ax.hist(
        clipped,
        bins=bins,
        color=theme.CLASS_COLORS[ReadClass.SHORT_DIVIDED],
        edgecolor=theme.SURFACE,
        linewidth=0.5,
        zorder=3,
    )

    ax.axvline(cut, color=theme.INK, linewidth=1.6, zorder=5)
    short = int((distances < cut).sum())
    long_ = int(distances.size - short)
    ax.text(
        cut,
        ax.get_ylim()[1] * 0.97,
        f"  cut-off {cut:,} bp",
        fontsize=_ANNOTATION_SIZE,
        family=theme.SANS,
        color=theme.INK,
        ha="left",
        va="top",
    )
    ax.text(
        cut * 0.85,
        ax.get_ylim()[1] * 0.80,
        f"{short:,} short  ",
        fontsize=_ANNOTATION_SIZE,
        family=theme.SANS,
        color=theme.CLASS_COLORS[ReadClass.SHORT_DIVIDED],
        ha="right",
        va="top",
    )
    ax.text(
        cut * 1.18,
        ax.get_ylim()[1] * 0.80,
        f"  {long_:,} long",
        fontsize=_ANNOTATION_SIZE,
        family=theme.SANS,
        color=theme.CLASS_COLORS[ReadClass.LONG_DIVIDED],
        ha="left",
        va="top",
    )

    ax.set_xscale("log")
    _axis_label(
        ax,
        xlabel="|reference span - read span|  (bp, log scale; zero shown at 1)",
        ylabel="Divided reads",
    )
    return fig


# --------------------------------------------------------------------------
# 5. Structural event map
# --------------------------------------------------------------------------


def figure_event_map(stats: Stats, reference: Reference) -> Figure:
    fig = _page()
    _heading(
        fig,
        "Structural events",
        "Every junction, inversion breakpoint and contig or circular join, on "
        "one page. This is ReadRift's main result without the reads.",
    )

    events = stats.events()
    contigs = [c for c in reference.contigs if c.length > 0][:MAX_EVENT_ROWS]
    if not events or not contigs:
        return _empty(fig, "No structural events were found.")

    ax = fig.add_axes((0.17, 0.16, 0.78, 0.66))
    _style(ax, grid_axis="none")
    ax.spines["left"].set_visible(False)

    longest = max(c.length for c in contigs)
    divisor, unit = _scaled_axis(longest)

    kinds = ("junction", "inversion", "circle", "contig_join")
    offsets = {kind: (index - 1.5) * 0.13 for index, kind in enumerate(kinds)}

    by_contig: dict[str, list] = {c.name: [] for c in contigs}
    for event in events:
        if event.contig in by_contig:
            by_contig[event.contig].append(event)

    for row, contig in enumerate(contigs):
        y = len(contigs) - row - 1
        ax.plot(
            [0, contig.length / divisor],
            [y, y],
            color=theme.GRID,
            linewidth=3.0,
            solid_capstyle="round",
            zorder=2,
        )
        for kind in kinds:
            positions = [
                e.position / divisor for e in by_contig[contig.name] if e.kind == kind
            ]
            if not positions:
                continue
            ax.plot(
                positions,
                [y + offsets[kind]] * len(positions),
                marker="|",
                markersize=11,
                markeredgewidth=1.6,
                linestyle="none",
                color=theme.EVENT_COLORS[kind],
                zorder=4,
            )

    # Contig names as tick labels, so matplotlib reserves room for them
    # instead of letting long accessions run off the page.
    ax.set_yticks(list(range(len(contigs))))
    ax.set_yticklabels(
        [c.name for c in reversed(contigs)],
        fontsize=_ANNOTATION_SIZE,
        family=theme.SANS,
    )
    for label in ax.get_yticklabels():
        label.set_color(theme.INK)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.7, len(contigs) - 0.3)
    ax.set_xlim(0, longest / divisor * 1.02)
    _axis_label(ax, xlabel=f"Position ({unit})")

    handles = [
        Line2D(
            [], [],
            marker="|", linestyle="none", markersize=11, markeredgewidth=1.8,
            color=theme.EVENT_COLORS[kind],
            label=f"{theme.EVENT_LABELS[kind]}  ({sum(1 for e in events if e.kind == kind):,})",
        )
        for kind in kinds
    ]
    legend = ax.legend(
        handles=handles, frameon=False, fontsize=_ANNOTATION_SIZE,
        loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=4,
    )
    for text in legend.get_texts():
        text.set_color(theme.INK_SECONDARY)
        text.set_family(theme.SANS)

    return fig


# --------------------------------------------------------------------------
# 6. Alignment identity
# --------------------------------------------------------------------------


def figure_identity(stats: Stats) -> Figure:
    fig = _page()
    _heading(
        fig,
        "Alignment identity",
        "Derived from the BTOP trace column, which the PostScript version read "
        "and never used. A low-identity tail separates genuine divergence from "
        "mapping artefacts.",
    )

    values = np.asarray(stats.identities, dtype=float) * 100.0
    if values.size == 0:
        return _empty(fig, "No BTOP trace data available (was --identity set?).")

    ax = fig.add_axes((0.055, 0.14, 0.90, 0.68))
    _style(ax)

    low = float(np.floor(values.min()))
    ax.hist(
        values,
        bins=np.linspace(low, 100.0, 60),
        color=theme.CLASS_COLORS[ReadClass.UNDIVIDED],
        edgecolor=theme.SURFACE,
        linewidth=0.5,
        zorder=3,
    )

    median = float(np.median(values))
    ax.axvline(median, color=theme.INK, linewidth=1.4, zorder=5)
    ax.text(
        median,
        ax.get_ylim()[1] * 0.97,
        f"  median {median:.2f}%",
        fontsize=_ANNOTATION_SIZE,
        family=theme.SANS,
        color=theme.INK,
        ha="left",
        va="top",
    )

    _axis_label(ax, xlabel="Identity (% of aligned columns)", ylabel="Reads")
    return fig


# --------------------------------------------------------------------------
# 7. Table view
# --------------------------------------------------------------------------


def figure_table(stats: Stats, reference: Reference) -> Figure:
    fig = _page()
    _heading(
        fig,
        "Per-contig summary",
        "The same numbers as the figures above, in text -- the accessible "
        "reading of every chart in this section.",
    )

    headers = [
        "Contig", "Length (bp)", "Coverage", "Undivided",
        "Short div.", "Long div.", "Inverted", "Events",
    ]
    widths = [0.22, 0.13, 0.09, 0.11, 0.11, 0.11, 0.10, 0.09]

    rows: list[list[str]] = []
    for contig in reference.contigs:
        cs = stats.per_contig.get(contig.name)
        if cs is None:
            continue
        rows.append([
            contig.name,
            f"{contig.length:,}",
            f"{cs.mean_coverage:.2f}x",
            f"{cs.counts[ReadClass.UNDIVIDED]:,}",
            f"{cs.counts[ReadClass.SHORT_DIVIDED]:,}",
            f"{cs.counts[ReadClass.LONG_DIVIDED]:,}",
            f"{cs.inverted:,}",
            f"{len(cs.events):,}",
        ])

    rows.append([
        "All contigs",
        f"{stats.reference_length:,}",
        f"{stats.mean_coverage:.2f}x",
        f"{stats.counts[ReadClass.UNDIVIDED]:,}",
        f"{stats.counts[ReadClass.SHORT_DIVIDED]:,}",
        f"{stats.counts[ReadClass.LONG_DIVIDED]:,}",
        f"{stats.inverted_reads:,}",
        f"{len(stats.events()):,}",
    ])

    ax = fig.add_axes((0.055, 0.08, 0.90, 0.76))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    xs: list[float] = []
    running = 0.0
    for width in widths:
        running += width
        xs.append(running)

    row_height = min(0.055, 0.90 / max(len(rows) + 2, 1))
    y = 0.95

    for index, header in enumerate(headers):
        ax.text(
            xs[index], y, header,
            fontsize=_ANNOTATION_SIZE, family=theme.SANS,
            color=theme.INK_SECONDARY, ha="right", va="bottom",
        )
    ax.plot([0, 1], [y - 0.012, y - 0.012], color=theme.AXIS, linewidth=0.9)

    for row_index, row in enumerate(rows):
        y -= row_height
        last = row_index == len(rows) - 1
        if last:
            ax.plot([0, 1], [y + row_height - 0.012, y + row_height - 0.012],
                    color=theme.AXIS, linewidth=0.9)
        for index, cell in enumerate(row):
            ax.text(
                xs[index], y, cell,
                fontsize=_ANNOTATION_SIZE, family=theme.SANS,
                color=theme.INK if last else theme.INK_SECONDARY,
                ha="right", va="bottom",
            )

    return fig


# --------------------------------------------------------------------------
# Section assembly
# --------------------------------------------------------------------------


def render_plots(
    stats: Stats, reference: Reference, params: Params
) -> Iterator[Figure]:
    """Yield the plot pages, in reading order."""
    yield figure_classification(stats)
    yield figure_coverage(stats, reference)
    yield figure_read_lengths(stats, params)
    yield figure_division_distance(stats, params)
    yield figure_event_map(stats, reference)
    if params.identity:
        yield figure_identity(stats)
    yield figure_table(stats, reference)
