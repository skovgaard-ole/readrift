"""Command-line interface.

The parser is generated from :data:`readrift.params.OPTIONS`, the same table the
PDF front page reads, so a documented option cannot fail to exist -- which is
what happened to ``--space2reads`` in the Perl (finding B18).

Every legacy short flag keeps its original letter and meaning.

There is one subcommand, ``browse``.  It is dispatched by looking at the first
argument rather than by an ``argparse`` subparser, because a subparser would
change the meaning of the existing positional interface: ``readrift ref.gb
reads.btop`` must keep working exactly as it does today.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from readrift.inputs.btop import OUTFMT
from readrift.params import OPTIONS, ExtractSpec, Params

#: Words that are treated as a subcommand when they appear first.  A reference
#: file is never called this, so the ambiguity is theoretical.
SUBCOMMANDS = frozenset({"browse"})

_EPILOG = f"""\
examples:
  readrift ref.gb reads.btop
  readrift ref.fa reads.btop -r 4000 -m 1000 -x 20
      include reads over 4 kb, each match over 1 kb, stop at 20x coverage
  readrift ref.fa reads.btop -o result --identity
      set the output prefix and add the alignment-identity figure
  readrift ref.fa reads.btop -e reads.fastq,chr1,50000,100000
      also write out the reads covering chr1:50000-100000
  readrift ref.fa reads.btop -e reads.fastq.gz,chr1,50000,100000 -e reads.fastq.gz,chr2,25000,75000

Every run also writes <prefix>.readriftdb.npz, the cache the interactive browser
reads.  Open it with:
  readrift browse result.readriftdb.npz

The BTOP file is made by NCBI BLAST+, which ReadRift needs but does not
install or run -- it reads the table BLAST writes:
  makeblastdb -in ref.fa -dbtype nucl -out refdb
  blastn -db refdb -query reads.fa -out reads.btop \\
    -outfmt "{OUTFMT}"
Copy the -outfmt string exactly. Do not add delim=\\t: tab is already BLAST's
separator, and BLAST would write the two characters \\t between the columns.
"""

_BROWSE_EPILOG = """\
examples:
  readrift browse result.readriftdb.npz
      open a cache written by an earlier run -- starts in about a second

  readrift browse ref.gb reads.btop -x 10
      build the cache first if it does not exist, then open it. The options are
      the same ones the normal run takes, and must match: a cache built with
      different filtering holds different reads, so it is rebuilt when they
      differ.

The server binds 127.0.0.1 and serves a read-only view of the cache. It has no
authentication, so do not point --host at a network interface.
"""

#: Suffixes stripped when deriving the default output prefix.
_STRIPPABLE = {".gz", ".btop", ".tsv", ".txt", ".out"}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="readrift",
        description=(
            "Visualise long sequence reads mapped onto a reference by BLAST, "
            "classified by how their alignment divides."
        ),
        epilog=_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "seq_file",
        help="reference sequence, GenBank or FASTA (may be gzipped)",
    )
    parser.add_argument(
        "btop_file",
        help="BLAST alignment in BTOP tabular format (may be gzipped)",
    )

    for option in OPTIONS:
        kwargs: dict = {
            "dest": option.key,
            "help": f"{option.help} (default: {option.default})",
            "default": option.default,
        }
        if option.kind == "int":
            kwargs["type"] = int
        elif option.kind == "float":
            kwargs["type"] = float
        if option.metavar:
            kwargs["metavar"] = option.metavar
        parser.add_argument(*option.flags, **kwargs)

    parser.add_argument(
        "-o", "--out", dest="out", default="",
        help="output prefix (default: derived from the BTOP filename)",
    )
    parser.add_argument(
        "-e", "--extract-reads", dest="extract_reads", action="append", default=[],
        metavar="FILE,CONTIG,START,END",
        help="write out the reads covering a region; may be repeated",
    )
    parser.add_argument(
        "-a", "--report-analysis", "--report_analysis",
        dest="report_analysis", action="store_true",
        help="also write a tab-separated summary for batch analysis",
    )
    parser.add_argument(
        "--identity", dest="identity", action="store_true",
        help="parse the BTOP trace column and add the alignment-identity figure",
    )
    parser.add_argument(
        "--no-plots", dest="no_plots", action="store_true",
        help="skip the summary figures and emit only the front page and maps",
    )
    parser.add_argument(
        "--no-cache", dest="no_cache", action="store_true",
        help="skip writing <prefix>.readriftdb.npz, the interactive browser's cache",
    )
    parser.add_argument(
        "--no-pdf", dest="no_pdf", action="store_true",
        help="write only the cache (and the -a report); skip the map document",
    )
    parser.add_argument(
        "--allow-unknown-contigs", dest="allow_unknown_contigs", action="store_true",
        help="warn instead of counting reads whose accession is not in the reference",
    )
    parser.add_argument(
        "--keep-fold-back", dest="keep_fold_back", action="store_true",
        help="keep HSPs that re-read reference the same read already mapped, in "
             "the opposite direction; they are dropped by default as an "
             "end-ligation artefact (see CHANGES.md)",
    )
    parser.add_argument(
        "--compat", dest="compat", action="store_true",
        help="restore the Perl's ambiguous behaviour where it differed from the "
             "documented options (see CHANGES.md)",
    )
    parser.add_argument(
        "-d", "--debug", dest="debug", action="store_true",
        help="print per-stage detail while running",
    )

    return parser


def default_prefix(btop_file: str) -> str:
    path = Path(btop_file)
    name = path.name
    while True:
        stem, dot, suffix = name.rpartition(".")
        if dot and f".{suffix}".lower() in _STRIPPABLE:
            name = stem
            continue
        break
    return name or path.stem or "readrift"


def parse_args(argv: Sequence[str] | None = None) -> Params:
    raw = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(raw)

    for name, path in (("reference", args.seq_file), ("BTOP", args.btop_file)):
        if not Path(path).exists():
            raise SystemExit(f"ERROR: {name} file not found: {path}")

    try:
        specs = tuple(ExtractSpec.parse(item) for item in args.extract_reads)
    except ValueError as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

    values = {
        option.key: getattr(args, option.key) for option in OPTIONS
    }

    params = Params(
        seq_file=args.seq_file,
        btop_file=args.btop_file,
        out=args.out or default_prefix(args.btop_file),
        extract_reads=specs,
        report_analysis=args.report_analysis,
        identity=args.identity,
        no_plots=args.no_plots,
        no_cache=args.no_cache,
        no_pdf=args.no_pdf,
        allow_unknown_contigs=args.allow_unknown_contigs,
        keep_fold_back=args.keep_fold_back,
        compat=args.compat,
        debug=args.debug,
        argv=tuple(raw),
        **values,
    )

    try:
        params.validate()
    except ValueError as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

    return params


# --------------------------------------------------------------------------
# browse
# --------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class BrowseArgs:
    """A resolved ``readrift browse`` invocation.

    Exactly one of *cache* and *params* describes where the data comes from:
    an existing cache file, or the inputs to build one from.
    """

    cache: Path | None
    params: Params | None
    host: str
    port: int
    open_browser: bool
    rebuild: bool
    verbose: bool


#: Recognised by extension, so `browse something.readriftdb.npz` needs no flag.
CACHE_SUFFIX = ".readriftdb.npz"


def build_browse_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="readrift browse",
        description=(
            "Open an interactive genome browser over a readrift run: pan and "
            "zoom the reference, with the reads coloured by how their alignment "
            "divides."
        ),
        epilog=_BROWSE_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="CACHE | REFERENCE BTOP",
        help=f"a {CACHE_SUFFIX} cache, or the reference and BTOP files to build one from",
    )

    for option in OPTIONS:
        kwargs: dict = {
            "dest": option.key,
            "help": f"{option.help} (default: {option.default})",
            "default": option.default,
        }
        if option.kind == "int":
            kwargs["type"] = int
        elif option.kind == "float":
            kwargs["type"] = float
        if option.metavar:
            kwargs["metavar"] = option.metavar
        parser.add_argument(*option.flags, **kwargs)

    parser.add_argument(
        "-o", "--out", dest="out", default="",
        help="output prefix, which decides the cache filename",
    )
    parser.add_argument(
        "--identity", dest="identity", action="store_true",
        help="parse the BTOP trace column so the detail panel can show identity",
    )
    parser.add_argument(
        "--allow-unknown-contigs", dest="allow_unknown_contigs", action="store_true",
        help="warn instead of failing when a BTOP accession is not in the reference",
    )
    parser.add_argument(
        "--keep-fold-back", dest="keep_fold_back", action="store_true",
        help="keep the fold-back tail of an end-ligation artefact; dropped by "
             "default (see CHANGES.md)",
    )
    parser.add_argument(
        "--compat", dest="compat", action="store_true",
        help="restore the Perl's ambiguous behaviour (see CHANGES.md)",
    )
    parser.add_argument(
        "--port", dest="port", type=int, default=0,
        help="port to listen on (default: any free port)",
    )
    parser.add_argument(
        "--host", dest="host", default="127.0.0.1",
        help="interface to bind (default: 127.0.0.1; anything else is unauthenticated)",
    )
    parser.add_argument(
        "--no-open", dest="open_browser", action="store_false",
        help="print the URL instead of opening a browser window",
    )
    parser.add_argument(
        "--rebuild", dest="rebuild", action="store_true",
        help="rebuild the cache even if the existing one is usable",
    )
    # No short flag: -v is already --space-to-reads in the shared option table,
    # and every legacy letter keeps its legacy meaning.
    parser.add_argument(
        "--verbose", dest="verbose", action="store_true",
        help="log every HTTP request",
    )
    return parser


def parse_browse_args(argv: Sequence[str]) -> BrowseArgs:
    parser = build_browse_parser()
    args = parser.parse_args(list(argv))
    inputs = [str(item) for item in args.inputs]

    common = {
        "host": args.host,
        "port": args.port,
        "open_browser": args.open_browser,
        "rebuild": args.rebuild,
        "verbose": args.verbose,
    }

    if len(inputs) == 1:
        cache = Path(inputs[0])
        if not cache.exists():
            raise SystemExit(
                f"ERROR: file not found: {cache}\n"
                f"       Give a {CACHE_SUFFIX} cache, or the reference and BTOP "
                f"files to build one from."
            )
        if not cache.name.endswith(CACHE_SUFFIX):
            raise SystemExit(
                f"ERROR: {cache} does not look like a readrift cache "
                f"(expected a name ending in {CACHE_SUFFIX}).\n"
                f"       To build one, pass the reference and BTOP files instead:\n"
                f"           readrift browse REFERENCE BTOP"
            )
        return BrowseArgs(cache=cache, params=None, **common)

    if len(inputs) > 2:
        raise SystemExit(
            f"ERROR: expected either one cache file or two input files, got "
            f"{len(inputs)}: {', '.join(inputs)}"
        )

    seq_file, btop_file = inputs
    for name, path in (("reference", seq_file), ("BTOP", btop_file)):
        if not Path(path).exists():
            raise SystemExit(f"ERROR: {name} file not found: {path}")

    values = {option.key: getattr(args, option.key) for option in OPTIONS}
    params = Params(
        seq_file=seq_file,
        btop_file=btop_file,
        out=args.out or default_prefix(btop_file),
        identity=args.identity,
        allow_unknown_contigs=args.allow_unknown_contigs,
        keep_fold_back=args.keep_fold_back,
        compat=args.compat,
        # Browsing never needs the poster; the cache is the whole point.
        no_pdf=True,
        no_cache=False,
        argv=tuple(argv),
        **values,
    )
    try:
        params.validate()
    except ValueError as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

    return BrowseArgs(cache=None, params=params, **common)
