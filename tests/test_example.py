"""The bundled example (examples/): real data, with its results pinned.

The synthetic fixtures in conftest.py have expected classes derived by hand
from the algorithm. These numbers are what the classifier *says* about real
data: the first 20x of Oxford Nanopore run DRR325755 against the M. iwaonis
SS37A-Re genome. They were checked equal, read by read, to a -x 20 run over the
original full alignment.

A failure here means the science changed. That needs an entry in CHANGES.md and
a person's review of the new numbers, not only an update to this file.
"""

from __future__ import annotations

import contextlib
import io
import statistics
from pathlib import Path

import pytest

from readrift.cli import parse_args
from readrift.models import ReadClass
from readrift.pipeline import Analysis, analyse

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
REFERENCE = EXAMPLES / "AP027148.gb.gz"
ALIGNMENT = EXAMPLES / "DRR325755_20x.btop.gz"

pytestmark = pytest.mark.skipif(
    not (REFERENCE.exists() and ALIGNMENT.exists()), reason="examples/ is not present"
)


@pytest.fixture(scope="module")
def analysis() -> Analysis:
    """One run for the whole module: about ten seconds, most of it the traces."""
    params = parse_args(
        [str(REFERENCE), str(ALIGNMENT), "-x", "20", "--identity", "--no-pdf", "--no-cache"]
    )
    with contextlib.redirect_stdout(io.StringIO()):
        return analyse(params)


def test_example_classification(analysis: Analysis) -> None:
    stats = analysis.stats
    assert stats.total_reads == 3_120
    assert stats.counts[ReadClass.UNDIVIDED] == 2_663
    assert stats.counts[ReadClass.SHORT_DIVIDED] == 340
    assert stats.counts[ReadClass.LONG_DIVIDED] == 117
    assert stats.inverted_reads == 12


def test_example_coverage_length_and_events(analysis: Analysis) -> None:
    stats = analysis.stats
    assert stats.mean_coverage == pytest.approx(19.6205, abs=1e-4)
    assert stats.n50 == 34_922
    assert len(stats.events()) == 318


def test_example_notes(analysis: Analysis) -> None:
    assert any(note.startswith("156 read(s) fold back") for note in analysis.notes)
    assert any("Stopped after 20.00x" in note for note in analysis.notes)


def test_example_identity_comes_from_the_trace_column(analysis: Analysis) -> None:
    """The example keeps BLAST's btop column precisely so --identity has data."""
    identities = list(analysis.stats.identities)
    assert len(identities) == 3_120
    assert statistics.median(identities) == pytest.approx(0.940, abs=0.001)
