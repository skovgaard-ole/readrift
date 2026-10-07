"""The front page: parameters, legend, project information, statistics, stamp.

Laid out on A3 landscape rather than the original 149 x 105 inch sheet, which
carried the same amount of text in one corner (decision D5).

Every parameter shown here comes from :data:`readrift.params.OPTIONS`, the same
table that builds the argument parser, so the two cannot drift apart the way
``--space2reads`` did in the Perl (finding B18).
"""

from __future__ import annotations

import platform
import socket
from datetime import datetime
from pathlib import Path

from matplotlib.backends.backend_pdf import FigureCanvasPdf
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure

from readrift.inputs.reference import Reference
from readrift.models import ReadClass
from readrift.params import Params
from readrift.render import theme
from readrift.stats import Stats

_LEFT = 58.0
_RIGHT = 630.0
_COLUMN_WIDTH = 500.0
_LINE = 15.5
_BODY_SIZE = 9.5
_SECTION_SIZE = 13.0


class _Cursor:
    """A downward text cursor, in page points."""

    __slots__ = ("ax", "x", "y", "_lines")

    def __init__(self, ax, x: float, y: float) -> None:
        self.ax = ax
        self.x = x
        self.y = y
        self._lines: list[list[tuple[float, float]]] = []

    def nl(self, n: float = 1.0) -> None:
        self.y -= _LINE * n

    def text(
        self,
        content: str,
        indent: float = 0.0,
        color: str = theme.INK_SECONDARY,
        size: float = _BODY_SIZE,
        family: str = theme.MONO,
    ) -> None:
        self.ax.text(
            self.x + indent, self.y, content,
            fontsize=size, family=family, color=color, ha="left", va="baseline",
        )
        self.nl()

    def row(self, left: str, right: str, indent: float = 12.0, right_x: float = 300.0,
            color: str = theme.INK_SECONDARY) -> None:
        self.ax.text(self.x + indent, self.y, left, fontsize=_BODY_SIZE,
                     family=theme.MONO, color=color, ha="left", va="baseline")
        self.ax.text(self.x + right_x, self.y, right, fontsize=_BODY_SIZE,
                     family=theme.MONO, color=color, ha="left", va="baseline")
        self.nl()

    def section(self, title: str) -> None:
        self.nl(1.0)
        self.ax.text(self.x, self.y, title, fontsize=_SECTION_SIZE,
                     family=theme.SANS, color=theme.INK, ha="left", va="baseline")
        self.nl(0.9)
        self._lines.append([(self.x, self.y + 4), (self.x + _COLUMN_WIDTH, self.y + 4)])
        self.nl(0.5)

    def flush(self) -> None:
        if self._lines:
            self.ax.add_collection(
                LineCollection(self._lines, colors=theme.GRID, linewidths=0.9)
            )
            self._lines = []


def _arrow(ax, style: str, x: float, y: float, color: str, scale: float = 0.75,
           rotation: float = 0.0) -> None:
    points = theme.transform(
        theme.arrow_points(style, 2.0), x, y, scale=scale, rotation=rotation
    )
    ax.add_collection(
        LineCollection([points], colors=color, linewidths=1.4, joinstyle="round")
    )
    if style == "arrow_circular":
        # The ring goes through the same transform as the polyline, so it stays
        # on the arrow's axis whatever the rotation.
        (ring_x, ring_y), = theme.transform(
            [(theme.CIRCULAR_RING_X, 0.0)], x, y, scale=scale, rotation=rotation
        )
        ax.plot([ring_x], [ring_y], marker="o", markersize=7,
                markerfacecolor=theme.SURFACE, markeredgecolor=color,
                markeredgewidth=1.6, linestyle="none")


def render_frontpage(
    stats: Stats,
    reference: Reference,
    params: Params,
    sample: str,
    notes: list[str],
) -> Figure:
    width, height = theme.A3_LANDSCAPE
    fig = Figure(figsize=(width / 72.0, height / 72.0), facecolor=theme.SURFACE)
    FigureCanvasPdf(fig)
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.set_axis_off()

    ax.text(_LEFT, height - 52, f"ReadRift of {sample}", fontsize=25,
            family=theme.SANS, color=theme.TITLE, ha="left", va="baseline")
    ax.text(_LEFT, height - 72, "Long sequence reads mapped by BLAST, classified by "
            "how their alignment divides", fontsize=11, family=theme.SANS,
            color=theme.INK_SECONDARY, ha="left", va="baseline")

    # ---- left column -----------------------------------------------------

    left = _Cursor(ax, _LEFT, height - 110)

    left.section("Parameters")
    for flag, value, default in params.as_frontpage_rows():
        left.row(flag, f"{value:>10}   (default {default})")

    left.section("Read lines")
    lines: list[list[tuple[float, float]]] = []
    colors: list[str] = []
    for cls in theme.CLASS_ORDER:
        lines.append([(left.x + 12, left.y + 3), (left.x + 92, left.y + 3)])
        colors.append(theme.CLASS_COLORS[cls])
        left.text(theme.CLASS_LABELS[cls], indent=110)
    lines.append([(left.x + 12, left.y + 3), (left.x + 92, left.y + 3)])
    colors.append(theme.INVERSION)
    left.text("Segment running opposite to the rest of its read (inversion)", indent=110)
    ax.add_collection(LineCollection(lines, colors=colors, linewidths=3.0))

    left.section("Arrows")
    for style, description in (
        ("arrow_fwd", "This read continues, forward direction"),
        ("arrow_rev", "This read continues, reverse direction"),
        ("arrow_plain", "This read joins an unexpected position"),
        ("arrow_ct", "This read joins the end of another contig"),
        ("arrow_circular", "This read joins the start and end of the same contig "
                           "-- indicates circularity"),
    ):
        _arrow(ax, style, left.x + 12, left.y + 3, theme.CLASS_COLORS[ReadClass.LONG_DIVIDED])
        left.text(description, indent=110)

    # The rotated arrow that marks an inverted segment.  `mapfig` draws this
    # shape at this angle on every map page, so the key has to be diagonal too
    # or it does not describe what is on the page.  Anchored higher than the
    # arrows above because rotating the shape carries it below its anchor.
    _arrow(ax, "inversion_arrow", left.x + 12, left.y + 9, theme.INVERSION,
           rotation=theme.INVERSION_ARROW_ROTATION)
    left.text("Indicator for inverted sequence", indent=110)

    if reference.is_genbank:
        left.section("Annotations")
        left.text("Coloured by their likelihood of causing genomic instability:", indent=12)
        for group, description in theme.ANNOTATION_LEGEND:
            left.text(description, indent=30, color=theme.ANNOTATION_COLORS[group])

    left.flush()

    # ---- right column ----------------------------------------------------

    right = _Cursor(ax, _RIGHT, height - 110)

    right.section("Project and assembly")
    rows = reference.metadata.rows()
    if rows:
        for name, value in rows:
            right.row(f"{name}:", _truncate(value, 58), right_x=190.0)
    else:
        right.text("Reference is FASTA; first line reads:", indent=12)
        right.text(_truncate(reference.metadata.first_line, 78), indent=30)

    right.section("Reference sequences")
    for contig in reference.contigs[:14]:
        right.row(contig.name, f"{contig.length:>15,} bases", right_x=250.0)
    if len(reference.contigs) > 14:
        right.text(f"... and {len(reference.contigs) - 14} more", indent=12)

    right.section("Reads mapped")
    right.row("Total reads mapped:", f"{stats.total_reads:>14,}", right_x=250.0,
              color=theme.INK)
    for cls in theme.CLASS_ORDER:
        right.row(f"{theme.CLASS_LABELS[cls]}:",
                  f"{stats.counts[cls]:>14,}    {stats.fraction(cls):6.2f}%",
                  right_x=250.0, color=theme.CLASS_COLORS[cls])
    right.row("Reads indicating inversions:",
              f"{stats.inverted_reads:>14,}    {stats.inverted_fraction:6.2f}%",
              right_x=250.0, color=theme.INVERSION)
    right.row("Mean coverage (matched bases):", f"{stats.mean_coverage:>13.2f}x",
              right_x=250.0)

    if notes:
        right.section("Notes")
        for note in notes:
            for chunk in _wrap(note, 78):
                right.text(chunk, indent=12, color=theme.INK_SECONDARY)

    right.section("Stamp")
    stamp = [
        ("Sample", sample),
        ("Analysis time", datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        ("Reference", params.seq_file),
        ("Reads", params.btop_file),
        ("Command", " ".join(params.argv) if params.argv else ""),
        ("Working directory", str(Path.cwd())),
        ("Host", socket.gethostname()),
        ("Platform", f"{platform.system()} {platform.release()} / Python "
                     f"{platform.python_version()}"),
    ]
    for name, value in stamp:
        if value:
            right.row(f"{name}:", _truncate(value, 62), right_x=150.0,
                      color=theme.INK_MUTED)

    right.flush()
    return fig


def _wrap(text: str, width: int) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if len(candidate) > width and current:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    return lines


def _truncate(text: str, width: int) -> str:
    return text if len(text) <= width else text[: width - 3] + "..."
