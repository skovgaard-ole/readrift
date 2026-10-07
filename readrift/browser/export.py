"""Rendering one visible region to PDF or PNG.

The point of the button is that a screenshot of the browser is not publishable
but a figure of the same region is.  It draws from the same
:class:`~readrift.browser.region.RegionView` the canvas drew from, in the same
colours, so "export" means what it says rather than "draw something similar".

The figure is built through ``Figure`` plus an explicit canvas class, the
pattern already used in :mod:`readrift.render.mapfig`: it needs no backend
selection and no global matplotlib state, which matters in a threaded server.
"""

from __future__ import annotations

import io
from itertools import pairwise

from matplotlib.collections import LineCollection
from matplotlib.figure import Figure

from readrift.browser.region import RegionView, continuation_direction
from readrift.render import theme

#: A3 landscape at 72 dpi, the same sheet the front page and plots use (D5).
FIGSIZE = (theme.A3_LANDSCAPE[0] / 72.0, theme.A3_LANDSCAPE[1] / 72.0)

_PNG_DPI = 200

#: Vertical share of the axes given to each band.
_COVERAGE_TOP = 1.00
_COVERAGE_BOTTOM = 0.86
_ANNOTATION_TOP = 0.83
_ANNOTATION_BOTTOM = 0.72
_PILEUP_TOP = 0.68
_PILEUP_BOTTOM = 0.02

#: Line width fed to :func:`readrift.render.theme.arrow_points`, and an overall
#: shrink applied to the result.
#:
#: The shapes are the printed map's, unchanged.  Their fixed offsets -- a 20 pt
#: shaft, a second chevron at 44 pt, the circular ring at 62 pt -- are sized for
#: a 10 750 pt sheet, where they are a rounding error.  On A3 they would be
#: three quarters of an inch, so the whole polyline is scaled down; the
#: proportions, and therefore the symbols, are identical.
_ARROW_LINE_WIDTH = 2.0
_ARROW_SCALE = 0.42

#: Stroke of the decorations, in points.  Kept well under a third of the
#: arrowhead length (``_ARROW_LINE_WIDTH * 2.5 * _ARROW_SCALE``, about 2.1 pt):
#: at this size a heavier stroke closes the chevron up into a solid blob, which
#: is what the canvas used to do before the same ratio was fixed there.
_ARROW_STROKE = 0.7


def render_region(
    view: RegionView,
    sample: str,
    fmt: str = "pdf",
    depth: dict | None = None,
) -> bytes:
    """Draw *view* and return the encoded document.

    *depth* is the coverage profile for the same window, as returned by
    :meth:`readrift.browser.store.BrowserStore.depth`.  It is optional only so
    the renderer stays testable without a store; the server always passes it,
    because a coverage track on screen and none in the export would make the
    button a lie.
    """
    fig = Figure(figsize=FIGSIZE, facecolor=theme.SURFACE)
    ax = fig.add_axes((0.045, 0.055, 0.93, 0.87))
    ax.set_facecolor(theme.SURFACE)
    ax.set_xlim(view.start, view.end)
    ax.set_ylim(0.0, 1.0)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(theme.AXIS)
    ax.tick_params(axis="x", colors=theme.AXIS, labelsize=9)
    ax.xaxis.set_major_formatter(lambda value, _pos: theme.thousands(value))
    ax.set_xlabel(f"Position on {view.contig} (bp)", color=theme.INK_SECONDARY, fontsize=10)

    _draw_title(fig, view, sample)
    _draw_coverage(ax, view, depth)
    _draw_annotations(ax, view)
    _draw_reads(ax, view)
    _draw_legend(ax, view)

    buffer = io.BytesIO()
    if fmt == "png":
        from matplotlib.backends.backend_agg import FigureCanvasAgg

        FigureCanvasAgg(fig)
        fig.savefig(buffer, format="png", dpi=_PNG_DPI, facecolor=theme.SURFACE)
    else:
        from matplotlib.backends.backend_pdf import FigureCanvasPdf

        FigureCanvasPdf(fig)
        fig.savefig(buffer, format="pdf", facecolor=theme.SURFACE)

    fig.clear()
    return buffer.getvalue()


# --------------------------------------------------------------------------


def _draw_title(fig: Figure, view: RegionView, sample: str) -> None:
    span = view.end - view.start
    shown = sum(view.shown.values())
    total = sum(view.counts.values())
    detail = f"{shown:,} reads drawn"
    if shown != total:
        detail = f"{shown:,} of {total:,} reads drawn (density limit)"

    fig.text(
        0.045,
        0.965,
        f"{sample} — {view.contig}:{theme.thousands(view.start)}–{theme.thousands(view.end)}",
        fontsize=15,
        family=theme.MONO,
        color=theme.TITLE,
        ha="left",
        va="baseline",
    )
    fig.text(
        0.045,
        0.938,
        f"{theme.thousands(span)} bp    {detail}",
        fontsize=9.5,
        family=theme.MONO,
        color=theme.INK_MUTED,
        ha="left",
        va="baseline",
    )


def _draw_coverage(ax, view: RegionView, depth: dict | None) -> None:
    """Depth of matched bases: the mean as a line, min-to-max as a band.

    The band is the point.  Averaged over a whole contig on screen, a narrow
    coverage trough disappears into the mean -- and a trough is exactly where a
    deletion or an excised element lives.
    """
    band_top, band_bottom = _COVERAGE_TOP, _COVERAGE_BOTTOM
    ax.axhline(band_bottom, color=theme.GRID, linewidth=0.8, zorder=1)

    positions = (depth or {}).get("positions") or []
    if not positions:
        return

    lows = depth["min"]
    means = depth["mean"]
    highs = depth["max"]
    peak = max(highs) or 1.0
    height = band_top - band_bottom

    def scale(value: float) -> float:
        return band_bottom + height * min(1.0, value / peak)

    ax.fill_between(
        positions,
        [scale(v) for v in lows],
        [scale(v) for v in highs],
        color=theme.CLASS_COLORS[theme.CLASS_ORDER[1]],
        alpha=0.22,
        linewidth=0.0,
        zorder=2,
    )
    ax.plot(
        positions,
        [scale(v) for v in means],
        color=theme.CLASS_COLORS[theme.CLASS_ORDER[1]],
        linewidth=1.0,
        zorder=3,
    )
    ax.text(
        view.start,
        band_top - height * 0.14,
        f"coverage  0–{theme.thousands(peak)}×",
        fontsize=8,
        family=theme.MONO,
        color=theme.INK_MUTED,
        ha="left",
        va="center",
    )


def _draw_annotations(ax, view: RegionView) -> None:
    if not view.annotations:
        return

    rows = max(a.row for a in view.annotations) + 1
    height = (_ANNOTATION_TOP - _ANNOTATION_BOTTOM) / max(2 * rows, 1)
    span = max(1, view.end - view.start)
    label_room = span * 0.012

    lines: dict[str, list] = {}
    for annotation in view.annotations:
        # + strand above the mid-line of the band, - strand below, matching the
        # printed map's convention.
        if annotation.strand > 0:
            y = _ANNOTATION_TOP - height * (annotation.row + 0.5)
        else:
            y = _ANNOTATION_BOTTOM + height * (annotation.row + 0.5)
        colour = theme.ANNOTATION_COLORS.get(
            annotation.group, theme.ANNOTATION_COLORS[""]
        )
        lines.setdefault(colour, []).append(
            [(annotation.start, y), (annotation.end, y)]
        )

        if annotation.name and annotation.end - annotation.start > label_room:
            symbol = theme.ANNOTATION_SYMBOLS.get(annotation.group, "")
            ax.text(
                (annotation.start + annotation.end) / 2.0,
                y + height * 0.32,
                f"{annotation.name}{symbol}",
                fontsize=5.5,
                family=theme.MONO,
                color=colour,
                ha="center",
                va="bottom",
                clip_on=True,
            )

    for colour, segments in lines.items():
        ax.add_collection(
            LineCollection(segments, colors=colour, linewidths=2.4, zorder=3)
        )

    ax.axhline(
        (_ANNOTATION_TOP + _ANNOTATION_BOTTOM) / 2.0,
        color=theme.GRID,
        linewidth=0.8,
        zorder=1,
    )


def _point_scales(view: RegionView) -> tuple[float, float]:
    """How many data units one PDF point is worth, on each axis.

    The arrow polylines are in points, because that is what the printed map
    works in.  Converting per axis -- rather than picking one factor -- means a
    chevron drawn here occupies exactly the same number of points on the page
    as the corresponding chevron in the report, despite x being measured in
    bases and y in a 0..1 band.
    """
    width_pt = FIGSIZE[0] * 72.0 * 0.93
    height_pt = FIGSIZE[1] * 72.0 * 0.87
    return (view.end - view.start) / width_pt, 1.0 / height_pt


def _mirrors(style: str) -> bool:
    """Whether *style* has to be flipped to point the other way.

    ``arrow_fwd`` and ``arrow_rev`` are already two shapes, one per direction --
    that is the whole difference between them -- so mirroring either would point
    it back the way it came.  The junction and inversion arrows are one shape
    used both ways, and do need it.
    """
    return continuation_direction(style) == 0


def _marker_polyline(
    style: str, x: float, y: float, direction: int, sx: float, sy: float
) -> list[tuple[float, float]]:
    """One decoration, in data coordinates, from the printed map's own shape."""
    if style == "inversion_arrow":
        points = theme.transform(
            theme.arrow_points("arrow_plain", _ARROW_LINE_WIDTH),
            0.0,
            0.0,
            scale=theme.INVERSION_ARROW_SCALE,
            rotation=theme.INVERSION_ARROW_ROTATION,
        )
        # Put the tip, not the tail, at the anchor -- the inverted piece's left
        # end -- nudged up and left so the arrow points at the piece without
        # covering it, as the printed map's placement in layout.py does.
        tip_x, tip_y = points[1]
        points = [(px - tip_x - 3.0, py - tip_y + 3.0) for px, py in points]
    else:
        points = theme.arrow_points(style, _ARROW_LINE_WIDTH)

    flip = -1.0 if (direction < 0 and _mirrors(style)) else 1.0
    return [
        (x + px * _ARROW_SCALE * flip * sx, y + py * _ARROW_SCALE * sy)
        for px, py in points
    ]


def _draw_markers(ax, view: RegionView, lane_height: float, y_of) -> None:
    """Arrowheads, junction ticks and circular-join rings.

    ``at_edge`` markers are drawn here without further thought, unlike in the
    canvas: the exporter asks for exactly the range it draws, so the window the
    region was built for and the window on the page are the same one.
    """
    sx, sy = _point_scales(view)
    tick_half = max(lane_height * 0.45, 0.25 * sy * 6.0)

    chevrons: dict[str, list] = {}
    ticks: dict[str, list] = {}
    rings: list[tuple[float, float, str]] = []

    for read in view.reads:
        if not read.markers:
            continue
        y = y_of(read)
        colour = theme.CLASS_COLORS[_class_of(read.cls)]

        for marker in read.markers:
            x = float(marker.x)
            if marker.style == "tick":
                ticks.setdefault(colour, []).append(
                    [(x, y - tick_half), (x, y + tick_half)]
                )
                continue

            ink = theme.INVERSION if marker.style == "inversion_arrow" else colour
            chevrons.setdefault(ink, []).append(
                _marker_polyline(marker.style, x, y, marker.direction, sx, sy)
            )
            if marker.style == "arrow_circular":
                flip = -1.0 if (marker.direction < 0 and _mirrors(marker.style)) else 1.0
                rings.append(
                    (x + theme.CIRCULAR_RING_X * _ARROW_SCALE * flip * sx, y, colour)
                )

    for colour, segments in ticks.items():
        ax.add_collection(
            LineCollection(segments, colors=colour, linewidths=_ARROW_STROKE, zorder=7)
        )
    for colour, polylines in chevrons.items():
        ax.add_collection(
            LineCollection(
                polylines,
                colors=colour,
                linewidths=_ARROW_STROKE,
                joinstyle="round",
                capstyle="round",
                zorder=8,
            )
        )
    for x, y, colour in rings:
        ax.plot(
            [x],
            [y],
            marker="o",
            markersize=4.0,
            markerfacecolor=theme.SURFACE,
            markeredgecolor=colour,
            markeredgewidth=_ARROW_STROKE,
            linestyle="none",
            zorder=9,
        )


def _draw_reads(ax, view: RegionView) -> None:
    if not view.reads:
        ax.text(
            (view.start + view.end) / 2.0,
            (_PILEUP_TOP + _PILEUP_BOTTOM) / 2.0,
            "no reads in this region",
            fontsize=11,
            color=theme.INK_MUTED,
            ha="center",
            va="center",
        )
        return

    plus_lanes = view.lanes.get("1", 0)
    minus_lanes = view.lanes.get("-1", 0)
    total_lanes = max(1, plus_lanes + minus_lanes + 1)
    band = _PILEUP_TOP - _PILEUP_BOTTOM
    lane_height = band / total_lanes
    divider = _PILEUP_TOP - lane_height * (plus_lanes + 0.5)

    width = max(0.5, min(3.0, lane_height * 260.0))

    def y_of(read) -> float:
        if read.strand > 0:
            return _PILEUP_TOP - lane_height * (read.lane - 0.5)
        return divider - lane_height * (read.lane + 0.5)

    by_colour: dict[str, list] = {}
    inversions: list = []
    joins: list[tuple[float, float, str]] = []

    for read in view.reads:
        y = y_of(read)
        colour = theme.CLASS_COLORS[_class_of(read.cls)]

        for segment in read.segments:
            line = [(segment.start, y), (segment.end, y)]
            if segment.inverted:
                inversions.append(line)
            else:
                by_colour.setdefault(colour, []).append(line)

        # Divided reads: a thin connector across the gap the read jumps over,
        # so the two pieces read as one molecule rather than two reads.
        ordered = sorted(read.segments, key=lambda s: s.start)
        for left, right in pairwise(ordered):
            if right.start > left.end:
                by_colour.setdefault(theme.GRID, []).append(
                    [(left.end, y), (right.start, y)]
                )

        if read.joins:
            joins.append((read.end, y, ", ".join(read.joins)))

    for colour, segments in by_colour.items():
        ax.add_collection(
            LineCollection(
                segments,
                colors=colour,
                linewidths=width if colour != theme.GRID else max(0.3, width * 0.35),
                capstyle="butt",
                zorder=4 if colour != theme.GRID else 3,
            )
        )

    if inversions:
        ax.add_collection(
            LineCollection(
                inversions,
                colors=theme.INVERSION,
                linewidths=width * 1.8,
                capstyle="butt",
                zorder=6,
            )
        )

    _draw_markers(ax, view, lane_height, y_of)

    if plus_lanes and minus_lanes:
        ax.axhline(divider, color=theme.GRID, linewidth=0.8, zorder=2)
    ax.text(
        view.start,
        _PILEUP_TOP + 0.012,
        "+ strand",
        fontsize=8,
        family=theme.MONO,
        color=theme.INK_MUTED,
        ha="left",
        va="bottom",
    )
    if minus_lanes:
        ax.text(
            view.start,
            divider - 0.012,
            "− strand",
            fontsize=8,
            family=theme.MONO,
            color=theme.INK_MUTED,
            ha="left",
            va="top",
        )

    # Only annotate joins when there are few enough for the labels to be read.
    if len(joins) <= 40:
        for x, y, text in joins:
            ax.text(
                x,
                y,
                f" →{text}",
                fontsize=5.5,
                family=theme.MONO,
                color=theme.EVENT_COLORS["contig_join"],
                ha="left",
                va="center",
                clip_on=True,
            )


def _draw_legend(ax, view: RegionView) -> None:
    from matplotlib.lines import Line2D

    handles = []
    for cls in theme.CLASS_ORDER:
        name = cls.value
        count = view.shown.get(name, 0)
        total = view.counts.get(name, 0)
        text = theme.CLASS_LABELS[cls]
        text += f"  {count:,}" if count == total else f"  {count:,} of {total:,}"
        handles.append(
            Line2D([], [], color=theme.CLASS_COLORS[cls], linewidth=3, label=text)
        )
    handles.append(
        Line2D([], [], color=theme.INVERSION, linewidth=3, label="Inverted segment")
    )

    legend = ax.legend(
        handles=handles,
        loc="lower right",
        bbox_to_anchor=(1.0, 1.005),
        ncol=4,
        frameon=False,
        fontsize=8.5,
        handlelength=1.6,
        columnspacing=1.6,
    )
    for text in legend.get_texts():
        text.set_color(theme.INK_SECONDARY)


def _class_of(name: str):
    from readrift.models import ReadClass

    return ReadClass(name)
