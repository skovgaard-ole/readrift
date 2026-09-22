"""Layout invariants and a full end-to-end run."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from readrift import layout as layout_module
from readrift import stats as stats_module
from readrift.classify import ClassifyStats, classify_stream
from readrift.inputs.btop import BtopStats, iter_read_groups
from readrift.inputs.reference import load_reference
from readrift.labels import LabelAssigner
from readrift.layout import split_across_pages
from readrift.models import ReadClass
from readrift.params import ExtractSpec, Params
from readrift.pipeline import run


def _pipeline_inputs(params: Params):
    reference = load_reference(params.seq_file)
    groups = iter_read_groups(params.btop_file, params.min_read_length, BtopStats())
    classified = classify_stream(
        groups, reference, params, LabelAssigner(), ClassifyStats()
    )
    summary, kept = stats_module.collect(classified, reference.contigs)
    return reference, summary, kept


# --------------------------------------------------------------------------
# Page splitting
# --------------------------------------------------------------------------


def test_split_within_one_page() -> None:
    assert split_across_pages(10.0, 90.0, 100.0) == [(0, 10.0, 90.0)]


def test_split_across_two_pages() -> None:
    assert split_across_pages(90.0, 150.0, 100.0) == [
        (0, 90.0, 100.0),
        (1, 0.0, 50.0),
    ]


def test_split_across_many_pages() -> None:
    """The Perl handled exactly one page break; a longer read drew wrong."""
    pieces = split_across_pages(50.0, 320.0, 100.0)
    assert [p[0] for p in pieces] == [0, 1, 2, 3]
    assert pieces[0] == (0, 50.0, 100.0)
    assert pieces[-1] == (3, 0.0, 20.0)


# --------------------------------------------------------------------------
# Lane packing
# --------------------------------------------------------------------------


def test_lanes_never_overlap(params: Params) -> None:
    reference, _summary, kept = _pipeline_inputs(params)
    layouts = layout_module.build_layout(kept, reference.contigs, params)

    for contig_layout in layouts.values():
        occupied: dict[tuple[int, float], list[tuple[float, float]]] = {}
        for placed in contig_layout.placed:
            for segment in placed.segments:
                key = (segment.page, segment.y)
                for x0, x1 in occupied.setdefault(key, []):
                    assert segment.x1 <= x0 or segment.x0 >= x1, (
                        f"overlap on page {segment.page} lane y={segment.y}"
                    )
                occupied[key].append((segment.x0, segment.x1))


def test_contig_joining_read_is_drawn_on_both_contigs(params: Params) -> None:
    reference, _summary, kept = _pipeline_inputs(params)
    layouts = layout_module.build_layout(kept, reference.contigs, params)

    drawn_on = {
        name
        for name, contig_layout in layouts.items()
        for placed in contig_layout.placed
        if placed.group.read == "read.6"
    }
    assert drawn_on == {"ctgA", "ctgB"}


def test_high_coverage_overflows_instead_of_crashing(params: Params) -> None:
    """Finding B03: the Perl walked off its lane array and died."""
    reference, _summary, kept = _pipeline_inputs(params)
    cramped = replace(params, line_space=370)  # k_max == 2, easily exhausted

    layouts = layout_module.build_layout(kept, reference.contigs, cramped)
    assert sum(item.overflow for item in layouts.values()) > 0


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------


def test_coverage_counts_matched_bases_only(params: Params) -> None:
    """Decision D3 / finding B10: no skipped reference counts as covered."""
    _reference, summary, kept = _pipeline_inputs(params)
    ctg_a = summary.per_contig["ctgA"]

    # Coverage is exactly the sum of aligned spans -- nothing more.
    manual = sum(
        h.sspan for group in kept for h in group.hits if h.contig == "ctgA"
    )
    assert ctg_a.matched_bases == manual

    # read.3 is short-divided across a 300 bp gap on ctgA. The Perl would have
    # added its whole 20000..24300 span (4300); only the two 2000 bp matches
    # count here, so the gap is absent from the total.
    read3 = next(g for g in kept if g.read == "read.3")
    assert read3.matched_bases == 4000
    assert read3.ref_high - read3.ref_low == 4300


def test_class_fractions_sum_to_100(params: Params) -> None:
    _reference, summary, _kept = _pipeline_inputs(params)
    total = sum(summary.fraction(cls) for cls in ReadClass)
    assert abs(total - 100.0) < 1e-6


# --------------------------------------------------------------------------
# End to end
# --------------------------------------------------------------------------


def test_full_run_writes_a_pdf(params: Params, fixture_dir: Path) -> None:
    exit_code = run(replace(params, report_analysis=True, identity=True))
    assert exit_code == 0

    pdf = Path(f"{params.out}.pdf")
    assert pdf.exists()
    assert pdf.stat().st_size > 10_000
    assert pdf.read_bytes().startswith(b"%PDF")

    tsv = Path(f"{params.out}_Analysis.tsv")
    assert tsv.exists()
    lines = tsv.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    assert header[0] == "record_type"
    # Finding B19: every row has the same width, unlike the Perl's 14-then-18.
    assert all(len(line.split("\t")) == len(header) for line in lines[1:])
    assert lines[1].split("\t")[0] == "run"


def test_full_run_against_genbank(params: Params, fixture_dir: Path) -> None:
    exit_code = run(
        replace(
            params,
            seq_file=str(fixture_dir / "ref.gb"),
            out=str(fixture_dir / "gb_out"),
        )
    )
    assert exit_code == 0
    assert Path(f"{fixture_dir / 'gb_out'}.pdf").exists()


def test_extraction_writes_a_list_and_sequences(
    params: Params, fixture_dir: Path
) -> None:
    """Finding B06: the Perl wrote the list file's header and nothing else."""
    spec = ExtractSpec(
        reads_file=str(fixture_dir / "reads.fastq.gz"),
        contig="ctgA",
        start=20_000,
        end=24_500,
    )
    exit_code = run(replace(params, extract_reads=(spec,), no_plots=True))
    assert exit_code == 0

    listing = fixture_dir / f"{spec.stem}.txt"
    body = listing.read_text(encoding="utf-8").splitlines()
    data = [line for line in body if not line.startswith("#")]
    assert data[0].split("\t") == ["read", "contig", "start", "end"]
    assert len(data) > 1, "the read list must contain the reads, not just a header"
    assert any(line.startswith("read.3") for line in data)

    sequences = fixture_dir / f"{spec.stem}.fastq"
    assert sequences.exists()
    assert sequences.read_text(encoding="utf-8").startswith("@")


def test_unknown_accessions_are_reported_not_silent(
    params: Params, fixture_dir: Path
) -> None:
    """Finding B11: a total mismatch produced a blank PDF and exit 0."""
    bad = fixture_dir / "bad.btop"
    bad.write_text(
        "someread\t1\t10\t2010\t1000\t3000\t5000\tNOT_A_CONTIG\t2000\n",
        encoding="utf-8",
    )
    exit_code = run(
        replace(params, btop_file=str(bad), out=str(fixture_dir / "bad_out"))
    )
    assert exit_code == 2


@pytest.mark.parametrize(
    "content",
    [
        pytest.param("", id="empty file"),
        pytest.param(
            "someread\\t1\\t10\\t2010\\t1000\\t3000\\t5000\\tctgA\\t2000\n",
            id="delim=\\t written literally",
        ),
    ],
)
def test_unparsable_btop_is_an_error_not_an_empty_map(
    params: Params, fixture_dir: Path, content: str, capsys: pytest.CaptureFixture[str]
) -> None:
    """A file with no usable line used to finish with exit 0 and an empty PDF."""
    bad = fixture_dir / "unparsable.btop"
    bad.write_text(content, encoding="utf-8")
    exit_code = run(
        replace(params, btop_file=str(bad), out=str(fixture_dir / "unparsable_out"))
    )
    assert exit_code == 2
    assert not (fixture_dir / "unparsable_out.pdf").exists()
    assert "blastn" in capsys.readouterr().out
