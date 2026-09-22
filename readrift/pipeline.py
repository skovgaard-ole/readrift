"""Run orchestration: input -> classification -> layout -> PDF (and cache).

The parsing front half is factored out as :func:`analyse` because two entry
points need it: the PDF run and the interactive browser.  On a real dataset
that half is the expensive one -- reading and classifying a multi-gigabyte BTOP
file takes minutes -- so it must never be written twice or run twice.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path

from readrift import extract, layout, report
from readrift import stats as stats_module
from readrift.classify import ClassifyStats, classify_stream
from readrift.inputs import btop as btop_reader
from readrift.inputs.reference import Reference, load_reference
from readrift.labels import LabelAssigner
from readrift.models import ReadClass, ReadGroup
from readrift.params import Params
from readrift.stats import Stats

_BANNER = (
    "\n" + "=" * 78 + "\n"
    "  readrift - visualising long sequence reads mapped by BLAST (BTOP)\n"
    + "=" * 78 + "\n"
)

#: Suffix of the browser cache written beside the PDF.
CACHE_SUFFIX = ".readriftdb.npz"


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def cache_path_for(params: Params) -> Path:
    """Where the browser cache for this run lives."""
    return Path(f"{params.out}{CACHE_SUFFIX}")


# --------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------


@dataclass(slots=True)
class Analysis:
    """The result of reading and classifying the inputs.

    Everything downstream -- layout, the PDF, the TSV report, the browser
    cache -- is derived from this and nothing else.
    """

    reference: Reference
    stats: Stats
    reads: list[ReadGroup]
    notes: list[str] = field(default_factory=list)
    index: extract.ReadIndex | None = None
    """Only built when ``-e`` was given (finding B25)."""


class AnalysisFailed(Exception):
    """Raised when the inputs cannot produce a meaningful result.

    Carries the process exit code the caller should return.
    """

    def __init__(self, code: int) -> None:
        super().__init__(f"analysis failed with exit code {code}")
        self.code = code


def analyse(params: Params) -> Analysis:
    """Read the reference and the alignments, and classify every read."""
    # ---- reference -------------------------------------------------------

    log(f"Reading reference {params.seq_file}")
    reference = load_reference(params.seq_file)
    print(
        f"    format: {'GenBank' if reference.is_genbank else 'FASTA'}   "
        f"contigs: {len(reference.contigs)}   "
        f"total: {reference.total_length:,} bp"
    )
    for contig in reference.contigs[:10]:
        print(f"      {contig.name:<28} {contig.length:>14,} bases")
    if len(reference.contigs) > 10:
        print(f"      ... and {len(reference.contigs) - 10} more")
    if reference.annotations:
        print(f"    annotations: {len(reference.annotations):,}")

    # ---- reads -----------------------------------------------------------

    log(f"Reading alignments {params.btop_file}")
    btop_stats = btop_reader.BtopStats()
    groups = btop_reader.iter_read_groups(
        params.btop_file, params.min_read_length, btop_stats
    )

    index: extract.ReadIndex | None = None
    if params.extract_reads:
        index = extract.ReadIndex(params.index_window, params.min_match_length)
        groups = extract.index_groups(groups, index)

    labeller = LabelAssigner()
    classify_stats = ClassifyStats()
    classified = classify_stream(
        groups, reference, params, labeller=labeller, stats=classify_stats
    )

    summary, kept = stats_module.collect(
        classified, reference.contigs, with_identity=params.identity
    )

    print(
        f"    BTOP lines: {btop_stats.lines:,}   "
        f"reads seen: {btop_stats.reads:,}   "
        f"below --min-read-length: {btop_stats.short_reads:,}"
    )

    # A file with no usable line -- a wrong -outfmt, or the wrong file -- used to
    # finish with exit 0 and an empty map: the accession check below never sees
    # an accession, so it has nothing to fail on.
    if btop_stats.nothing_parsed:
        explanation = btop_stats.explain_nothing_parsed(params.btop_file)
        print("\nERROR: " + explanation.replace("\n", "\n       ") + "\n")
        raise AnalysisFailed(2)

    # A total accession mismatch between the BTOP file and the reference is the
    # commonest way to get an empty result, and the Perl gave no hint that it
    # had happened (finding B11).
    if classify_stats.unknown_contigs and not summary.total_reads:
        found = ", ".join(sorted(classify_stats.unknown_contigs)[:5])
        expected = ", ".join(reference.names[:5])
        print(
            f"\nERROR: none of the accessions in {params.btop_file} match the "
            f"reference.\n"
            f"       BTOP file has:  {found}\n"
            f"       Reference has:  {expected}\n"
            f"       These usually differ only by a version suffix; rebuild the "
            f"BLAST database\n"
            f"       from the same file, or pass --allow-unknown-contigs to "
            f"continue anyway.\n"
        )
        if not params.allow_unknown_contigs:
            raise AnalysisFailed(2)

    notes: list[str] = []
    notes.extend(labeller.warnings())
    notes.extend(btop_stats.warnings())
    notes.extend(classify_stats.warnings(reference))

    return Analysis(
        reference=reference,
        stats=summary,
        reads=kept,
        notes=notes,
        index=index,
    )


def write_cache(analysis: Analysis, params: Params) -> Path | None:
    """Write the browser cache.  Never fatal -- warns and returns ``None``.

    Called after the PDF so that a failure here cannot cost the user the
    document they were waiting minutes for.
    """
    from readrift.browser.store import build_store

    path = cache_path_for(params)
    try:
        build_store(
            path,
            analysis.reference,
            analysis.reads,
            analysis.stats,
            params,
            notes=analysis.notes,
        )
    except Exception as exc:
        # Deliberately broad.  On a real dataset the analysis above cost
        # minutes; nothing that goes wrong writing an auxiliary file is worth
        # discarding it for.  Run with -d to see the traceback instead.
        if params.debug:
            raise
        print(
            f"    WARNING: could not write the browser cache {path}\n"
            f"             {type(exc).__name__}: {exc}\n"
            f"             The PDF is unaffected. Re-run with -d for the traceback."
        )
        return None
    return path


# --------------------------------------------------------------------------
# The PDF run
# --------------------------------------------------------------------------


def run(params: Params) -> int:
    """Execute a full run.  Returns a process exit code."""
    started = time.perf_counter()
    print(_BANNER)

    try:
        analysis = analyse(params)
    except AnalysisFailed as exc:
        return exc.code

    reference = analysis.reference
    summary = analysis.stats
    notes = analysis.notes

    # ---- layout ----------------------------------------------------------

    pdf_path: Path | None = None
    if not params.no_pdf:
        log("Placing reads")
        layouts = layout.build_layout(analysis.reads, reference.contigs, params)
        overflow = sum(item.overflow for item in layouts.values())
        total_pages = sum(item.pages for item in layouts.values())

        if overflow:
            notes.append(
                f"{overflow:,} read(s) found no free lane and were drawn on the "
                f"outermost one. Raise --page-scale or lower --cov-max to see them "
                f"separately. (The Perl crashed at this point -- finding B03.)"
            )

        # ---- output ------------------------------------------------------

        from readrift.render import write_pdf

        sample = Path(params.out).name
        pdf_path = Path(f"{params.out}.pdf")
        log(f"Writing {pdf_path} ({total_pages} map page(s))")
        try:
            pages = write_pdf(
                pdf_path, summary, reference, layouts, params, sample, notes
            )
        except PermissionError:
            print(
                f"\nERROR: cannot write {pdf_path} -- it is probably open in another "
                f"program. Close it and run again, or pass -o with another prefix.\n"
            )
            return 1

    if params.report_analysis:
        sample = Path(params.out).name
        tsv = report.write_report(
            f"{params.out}_Analysis.tsv", summary, reference, params, sample
        )
        log(f"Wrote {tsv}")

    for spec in params.extract_reads:
        assert analysis.index is not None
        requested, written, path = extract.run_extraction(spec, analysis.index, params)
        log(
            f"Extracted {written:,} of {requested:,} reads for "
            f"{spec.contig}:{spec.start:,}-{spec.end:,} -> {path}"
        )
        if written < requested:
            print(
                f"    note: {requested - written:,} indexed read(s) were not found "
                f"in {spec.reads_file}"
            )

    cache = None
    if not params.no_cache:
        log("Writing the browser cache")
        cache = write_cache(analysis, params)

    # ---- summary ---------------------------------------------------------

    print("\nReads mapped:")
    for cls, name in (
        (ReadClass.UNDIVIDED, "Undivided reads"),
        (ReadClass.SHORT_DIVIDED, "Divided reads, short distance"),
        (ReadClass.LONG_DIVIDED, "Divided reads, long distance"),
    ):
        print(
            f"    {name:<32}{summary.counts[cls]:>12,}"
            f"{summary.fraction(cls):>9.2f}%"
        )
    print(
        f"    {'Reads indicating inversions':<32}"
        f"{summary.inverted_reads:>12,}{summary.inverted_fraction:>9.2f}%"
    )
    print(f"    {'Total':<32}{summary.total_reads:>12,}")
    print(f"\n    Mean coverage (matched bases): {summary.mean_coverage:.2f}x")
    print(f"    Read N50: {summary.n50:,} bp")

    if notes:
        print("\nNotes:")
        for note in notes:
            print(f"  - {note}")

    elapsed = time.perf_counter() - started
    if pdf_path is not None:
        print(f"\nWrote {pdf_path} ({pages} pages) in {elapsed:.1f}s")
    else:
        print(f"\nFinished in {elapsed:.1f}s")

    if cache is not None:
        print(f"Wrote {cache}\n  Browse it with:  python -m readrift browse {cache}\n")
    else:
        print()

    return 0


# --------------------------------------------------------------------------
# The interactive browser
# --------------------------------------------------------------------------


def browse(args) -> int:
    """Serve the interactive browser.  Takes a :class:`readrift.cli.BrowseArgs`."""
    from readrift.browser.server import serve
    from readrift.browser.store import BrowserStore, StoreFormatError

    store: BrowserStore | None = None

    if args.cache is not None:
        try:
            store = BrowserStore.open(args.cache)
        except StoreFormatError as exc:
            # A cache path on its own carries no way to rebuild -- the inputs
            # are not on the command line -- so say what would.
            print(f"\nERROR: {exc}")
            print(
                "       A cache is rebuilt from its inputs, which this form of "
                "the command\n"
                "       does not name. Re-run with them:\n"
                "           readrift browse <reference> <btop>\n"
            )
            return 1
    else:
        params = args.params
        assert params is not None
        path = cache_path_for(params)

        if path.exists() and not args.rebuild:
            try:
                candidate = BrowserStore.open(path)
            except StoreFormatError as exc:
                print(f"  {path} cannot be reused ({exc}); rebuilding.")
            else:
                stale, reason = candidate.is_stale(params)
                if stale:
                    print(f"  {path} is out of date ({reason}); rebuilding.")
                else:
                    log(f"Reusing {path}")
                    store = candidate

        if store is None:
            print(_BANNER)
            try:
                analysis = analyse(params)
            except AnalysisFailed as exc:
                return exc.code
            log("Writing the browser cache")
            written = write_cache(analysis, params)
            if written is None:
                return 1
            try:
                store = BrowserStore.open(written)
            except StoreFormatError as exc:
                print(f"\nERROR: {exc}\n")
                return 1

    return serve(
        store,
        host=args.host,
        port=args.port,
        open_browser=args.open_browser,
        verbose=args.verbose,
    )
