"""The ReadRift pages.

One figure per contig page, in the original poster format: the page is
``10750 x 7600`` points (times ``--page-scale``) and the axes fill it exactly,
so one data unit is one point.  That keeps every offset, font size and line
width from the PostScript meaningful without rescaling.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

from matplotlib.backends.backend_pdf import FigureCanvasPdf
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from readrift.labels import gene_label
from readrift.models import Annotation, ContigLayout, SegmentKind
from readrift.params import Params
from readrift.render import theme

#: Horizontal distance, in drawing units, within which two annotations are
#: considered crowded enough to stack on separate rows.
_ANNOTATION_CROWDING = 10.0


def _new_page(params: Params) -> tuple[Figure, Any]:
    fig = Figure(
        figsize=(params.page_width / 72.0, params.page_height / 72.0),
        facecolor=theme.SURFACE,
    )
    FigureCanvasPdf(fig)
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(-250.0 * params.page_scale, 10_500.0 * params.page_scale)
    ax.set_ylim(-0.5 * params.page_height, 0.5 * params.page_height)
    ax.set_axis_off()
    ax.set_facecolor(theme.SURFACE)
    return fig, ax


def _draw_frame(ax, params: Params) -> None:
    ax.add_patch(
        Rectangle(
            (-250.0 * params.page_scale, -0.5 * params.page_height),
            10_750.0 * params.page_scale,
            params.page_height,
            fill=False,
            edgecolor=theme.FRAME,
            linewidth=params.line_width * 3.0,
            zorder=0,
        )
    )


def _draw_axis(ax, layout: ContigLayout, page: int, params: Params) -> None:
    """The reference axis with its three tiers of tick marks."""
    axis_width = params.axis_width
    length_units = layout.contig.length / params.map_scale
    axis_end = min(axis_width, length_units - page * axis_width)
    if axis_end <= 0:
        return

    ax.add_collection(
        LineCollection(
            [[(0.0, 0.0), (axis_end, 0.0)]],
            colors=theme.AXIS,
            linewidths=params.line_width,
            capstyle="round",
            zorder=2,
        )
    )

    tiers = (
        (1_000.0 / params.map_scale, 0.0, 3.0, params.line_width / 6.0),
        (10_000.0 / params.map_scale, -4.0, 4.0, params.line_width / 2.0),
        (50_000.0 / params.map_scale, -6.0, 6.0, params.line_width),
    )
    for step, y0, y1, width in tiers:
        if step <= 0:
            continue
        count = int(axis_end / step)
        if count > 20_000:
            continue
        marks = [
            [(j * step, y0), (j * step, y1)] for j in range(1, count + 1)
        ]
        if marks:
            ax.add_collection(
                LineCollection(
                    marks, colors=theme.AXIS, linewidths=width, zorder=2
                )
            )

    # Coordinate labels on the coarse tier.
    coarse = 50_000.0 / params.map_scale
    if coarse > 0:
        for j in range(int(axis_end / coarse) + 1):
            position = page * params.map_scale * axis_width + j * 50_000.0
            ax.text(
                j * coarse,
                20.0,
                theme.thousands(position),
                fontsize=theme.AXIS_FONT_SIZE,
                family=theme.MONO,
                color=theme.AXIS,
                ha="center",
                va="baseline",
                zorder=3,
            )

    # Where the contig stops mid-page, mark its true end.
    if axis_end < axis_width:
        ax.text(
            axis_end + 30.0,
            -6.0,
            theme.thousands(layout.contig.length),
            fontsize=theme.AXIS_FONT_SIZE,
            family=theme.MONO,
            color=theme.AXIS,
            ha="left",
            va="baseline",
            zorder=3,
        )


def _draw_annotations(
    ax, annotations: list[Annotation], page: int, params: Params
) -> None:
    """The GenBank feature track, in up to three rows either side of the axis.

    Unlike the Perl, a feature crossing a page boundary is drawn on both pages
    instead of running off the edge (finding B23), and the row counters reset
    per contig instead of carrying over between them.
    """
    axis_width = params.axis_width
    rows = {1: 0, -1: 0}
    last_end = {1: 0.0, -1: 0.0}

    lines: dict[str, list[list[tuple[float, float]]]] = {}

    for annotation in annotations:
        begin = annotation.start / params.map_scale
        end = annotation.end / params.map_scale
        strand = 1 if annotation.strand >= 0 else -1

        if begin - last_end[strand] > _ANNOTATION_CROWDING:
            rows[strand] = 0
        base = theme.ANNOTATION_Y_PLUS if strand > 0 else theme.ANNOTATION_Y_MINUS
        y = base - theme.ANNOTATION_ROW_HEIGHT * rows[strand]
        if end - begin < _ANNOTATION_CROWDING:
            rows[strand] = (rows[strand] + 1) % theme.ANNOTATION_ROWS
            last_end[strand] = end

        page_start = page * axis_width
        page_end = page_start + axis_width
        if end < page_start or begin > page_end:
            continue

        x0 = max(begin, page_start) - page_start
        x1 = min(end, page_end) - page_start

        color = theme.ANNOTATION_COLORS.get(annotation.group, theme.ANNOTATION_COLORS[""])
        lines.setdefault(color, []).append([(x0, y), (x1, y)])

        if begin >= page_start:
            # Exactly what the browser draws -- one rule for a gene's name on
            # every surface.  The Perl instead cut every name to its last five
            # characters (`substr($name,-5)`, line 1141, uncommented, in the
            # middle of its drawing loop). That is right for a locus tag only
            # by accident, because a tag ends in its number; on anything else
            # it eats the name from the front, turning `fadE16_1` into `E16_1`
            # and stripping the series letter off `t00010`.
            name = gene_label(annotation.name)
            name += theme.ANNOTATION_SYMBOLS.get(annotation.group, "")
            if name:
                ax.text(
                    x0,
                    y + (2.0 if strand > 0 else -3.5),
                    name,
                    fontsize=theme.ANNOTATION_FONT_SIZE,
                    family=theme.MONO,
                    color=color,
                    ha="left",
                    va="baseline",
                    zorder=3,
                )

    for color, segments in lines.items():
        ax.add_collection(
            LineCollection(
                segments,
                colors=color,
                linewidths=params.line_width / 3.0,
                zorder=2,
            )
        )


def _draw_reads(ax, layout: ContigLayout, page: int, params: Params) -> None:
    inversion_width = min(params.line_width * 4.0, params.line_space * 0.6)

    read_lines: dict[str, list[list[tuple[float, float]]]] = {}
    inversion_lines: list[list[tuple[float, float]]] = []
    arrow_lines: dict[str, list[list[tuple[float, float]]]] = {}
    tick_lines: dict[str, list[list[tuple[float, float]]]] = {}
    rings: list[tuple[float, float, str]] = []
    texts: list[tuple[float, float, str, str]] = []

    tick_height = 2.0 * params.line_width

    for placed in layout.placed:
        color = theme.CLASS_COLORS[placed.group.cls]

        for segment in placed.segments:
            if segment.page != page:
                continue
            line = [(segment.x0, segment.y), (segment.x1, segment.y)]
            if segment.kind is SegmentKind.INVERSION:
                inversion_lines.append(line)
            else:
                read_lines.setdefault(color, []).append(line)

        for marker in placed.markers:
            if marker.page != page:
                continue
            if marker.style == "tick":
                tick_lines.setdefault(color, []).append(
                    [
                        (marker.x, marker.y - tick_height),
                        (marker.x, marker.y + tick_height),
                    ]
                )
                continue
            if marker.style == "inversion_arrow":
                points = theme.transform(
                    theme.arrow_points("arrow_plain", params.line_width),
                    marker.x,
                    marker.y,
                    scale=theme.INVERSION_ARROW_SCALE,
                    rotation=theme.INVERSION_ARROW_ROTATION,
                )
                arrow_lines.setdefault(theme.INVERSION, []).append(points)
                continue

            points = theme.transform(
                theme.arrow_points(marker.style, params.line_width),
                marker.x,
                marker.y,
            )
            arrow_lines.setdefault(color, []).append(points)
            if marker.style == "arrow_circular":
                rings.append((marker.x + theme.CIRCULAR_RING_X, marker.y, color))

        for label in placed.labels:
            if label.page != page:
                continue
            texts.append((label.x, label.y, label.text, color))

    for color, segments in read_lines.items():
        ax.add_collection(
            LineCollection(
                segments,
                colors=color,
                linewidths=params.line_width,
                capstyle="butt",
                zorder=4,
            )
        )

    if inversion_lines:
        ax.add_collection(
            LineCollection(
                inversion_lines,
                colors=theme.INVERSION,
                linewidths=inversion_width,
                capstyle="butt",
                zorder=5,
            )
        )

    for color, segments in tick_lines.items():
        ax.add_collection(
            LineCollection(
                segments,
                colors=color,
                linewidths=params.line_width / 3.0,
                zorder=6,
            )
        )

    for color, polylines in arrow_lines.items():
        ax.add_collection(
            LineCollection(
                polylines,
                colors=color,
                linewidths=params.line_width,
                joinstyle="round",
                zorder=6,
            )
        )

    for x, y, color in rings:
        ax.plot(
            [x],
            [y],
            marker="o",
            markersize=params.line_width * 4.0,
            markerfacecolor=theme.SURFACE,
            markeredgecolor=color,
            markeredgewidth=params.line_width,
            linestyle="none",
            zorder=7,
        )

    for x, y, text, color in texts:
        ax.text(
            x,
            y,
            text,
            fontsize=theme.LABEL_FONT_SIZE,
            family=theme.MONO,
            color=color,
            ha="left",
            va="baseline",
            zorder=8,
        )


def render_contig(
    layout: ContigLayout,
    annotations: list[Annotation],
    params: Params,
    sample: str,
) -> Iterator[Figure]:
    """Yield one figure per page of *layout*."""
    for page in range(layout.pages):
        fig, ax = _new_page(params)
        _draw_frame(ax, params)

        first_base = int(page * params.axis_width * params.map_scale)
        last_base = min(
            int((page + 1) * params.axis_width * params.map_scale),
            layout.contig.length,
        )
        ax.text(
            0.0,
            0.5 * params.page_height - 70.0,
            f"ReadRift of {sample} on {layout.contig.name}."
            f"      Page {page + 1} of {layout.pages}."
            f"      Bases {theme.thousands(first_base)} to {theme.thousands(last_base)}.",
            fontsize=theme.TITLE_FONT_SIZE,
            family=theme.MONO,
            color=theme.TITLE,
            ha="left",
            va="baseline",
            zorder=9,
        )

        _draw_axis(ax, layout, page, params)
        if annotations:
            _draw_annotations(ax, annotations, page, params)
        _draw_reads(ax, layout, page, params)

        yield fig
