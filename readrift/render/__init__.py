"""PDF assembly.

The Perl wrote one PostScript file per contig page, concatenated them with a
``while (<>)`` loop, then shelled out to ``ps2pdf``.  That required Ghostscript,
littered the working directory with ``.ps`` intermediates that were never
cleaned up, and -- because no label was ever escaped -- produced a broken
document if any read name, product or file path contained a parenthesis
(findings B09, B29).

Here everything goes straight into one PDF.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from readrift.inputs.reference import Reference
from readrift.models import Annotation, ContigLayout
from readrift.params import Params
from readrift.stats import Stats

__all__ = ["write_pdf"]

# matplotlib and the page modules are imported inside write_pdf, not here.
# `from readrift.render import theme` executes this file, and the interactive
# browser does exactly that just to read the palette -- it should not pay a
# second of matplotlib import to find out what colour an inversion is.


def _annotations_for(annotations: Iterable[Annotation], contig: str) -> list[Annotation]:
    selected = [a for a in annotations if a.contig == contig]
    selected.sort(key=lambda a: (a.start, a.end))
    return selected


def write_pdf(
    path: str | Path,
    stats: Stats,
    reference: Reference,
    layouts: dict[str, ContigLayout],
    params: Params,
    sample: str,
    notes: list[str],
) -> int:
    """Write the whole document.  Returns the number of pages."""
    from matplotlib.backends.backend_pdf import PdfPages

    from readrift.render.frontpage import render_frontpage
    from readrift.render.mapfig import render_contig
    from readrift.render.plots import render_plots

    pages = 0
    with PdfPages(path) as pdf:
        info = pdf.infodict()
        info["Title"] = f"ReadRift of {sample}"
        info["Subject"] = (
            f"{stats.total_reads:,} long reads mapped onto "
            f"{len(reference.contigs)} contig(s), {stats.reference_length:,} bp"
        )
        info["Creator"] = "readrift (Python port of read_print_23.pl)"

        front = render_frontpage(stats, reference, params, sample, notes)
        pdf.savefig(front)
        front.clear()
        pages += 1

        if not params.no_plots:
            for figure in render_plots(stats, reference, params):
                pdf.savefig(figure)
                figure.clear()
                pages += 1

        for contig in reference.contigs:
            layout = layouts.get(contig.name)
            if layout is None:
                continue
            annotations = _annotations_for(reference.annotations, contig.name)
            for figure in render_contig(layout, annotations, params, sample):
                pdf.savefig(figure)
                figure.clear()
                pages += 1

    return pages
