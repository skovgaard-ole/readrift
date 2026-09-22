"""Parser tests, including the specific things the Perl got wrong."""

from __future__ import annotations

from pathlib import Path

import pytest

from readrift.btop_trace import parse as parse_btop
from readrift.inputs.btop import BtopStats, iter_read_groups, parse_line
from readrift.inputs.fasta import accession_from_header
from readrift.inputs.reads import iter_records
from readrift.inputs.reference import load_reference
from readrift.labels import gene_label
from tests.conftest import CONTIG_LENGTHS


def test_fasta_lengths(fixture_dir: Path) -> None:
    reference = load_reference(fixture_dir / "ref.fa")
    assert {c.name: c.length for c in reference.contigs} == CONTIG_LENGTHS


def test_fasta_counts_ambiguity_codes(fixture_dir: Path) -> None:
    """Finding B22: the Perl deleted N, R, Y ... before measuring."""
    reference = load_reference(fixture_dir / "ref.fa")
    # The generated sequence contains N, R and Y; if they were stripped the
    # lengths would come out short.
    assert reference.contigs[0].length == CONTIG_LENGTHS["ctgA"]


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        (">NZ_CP012345.1 Escherichia coli", "NZ_CP012345.1"),
        (">gi|12345|ref|NC_000913.3|", "NC_000913.3"),
        (">lcl|contig_7 something", "contig_7"),
        (">plain", "plain"),
    ],
)
def test_fasta_header_accessions(header: str, expected: str) -> None:
    """Finding B21: greedy \\S+ after a pipe swallowed the whole header."""
    accession, _aliases = accession_from_header(header)
    assert accession == expected


def test_genbank_lengths_and_metadata(fixture_dir: Path) -> None:
    reference = load_reference(fixture_dir / "ref.gb")
    assert {c.name: c.length for c in reference.contigs} == CONTIG_LENGTHS

    meta = reference.metadata
    assert meta.organism == "Testus exemplaris"
    assert meta.bio_project == "PRJNA000001"
    assert meta.bio_sample == "SAMN00000001"
    assert meta.sra == "SRR0000001"
    assert meta.assembly_method == "Flye v. 2.9"
    assert meta.sequencing_technology == "Oxford Nanopore (MinION)"


def test_genbank_versioned_accessions_resolve(fixture_dir: Path) -> None:
    """Finding B11: a BLAST db built with VERSION never matched LOCUS."""
    reference = load_reference(fixture_dir / "ref.gb")
    assert reference.resolve("ctgA") == "ctgA"
    assert reference.resolve("ctgA.1") == "ctgA"
    assert reference.resolve("nonesuch") is None


def test_genbank_feature_classification(fixture_dir: Path) -> None:
    reference = load_reference(fixture_dir / "ref.gb")
    by_name = {a.name: a for a in reference.annotations}

    assert by_name["thrA"].group == ""
    assert by_name["thrA"].strand == 1
    assert by_name["thrA"].start == 1000
    assert by_name["thrA"].end == 2500

    assert by_name["TST_0002"].group == "hyp"
    assert by_name["TST_0002"].strand == -1
    assert by_name["TST_0003"].group == "RNA"
    assert by_name["TST_0003"].feature == "tRNA"
    assert by_name["TST_0004"].group == "IS"       # wrapped /product, joined
    assert by_name["TST_0005"].group == "pseudo"
    assert by_name["TST_0006"].group == "phage"

    # Finding B24: a standalone feature is no longer folded into the gene
    # before it.
    assert any(a.feature == "repeat_region" for a in reference.annotations)
    assert any(a.contig == "ctgB" and a.name == "repB" for a in reference.annotations)


def test_genbank_keeps_the_qualifiers_the_browser_shows(fixture_dir: Path) -> None:
    """/product, /locus_tag and /note survive parsing, not just grouping."""
    reference = load_reference(fixture_dir / "ref.gb")
    by_name = {a.name: a for a in reference.annotations}

    # /gene wins the name; the locus tag is still reachable beside it.
    assert by_name["thrA"].locus_tag == "TST_0001"
    assert by_name["thrA"].product == "aspartokinase"

    # No /gene at all, so the name *is* the locus tag.
    assert by_name["TST_0002"].locus_tag == "TST_0002"
    assert by_name["TST_0002"].product == "hypothetical protein"

    # A wrapped /product is joined before it is stored, as it is before it is
    # matched for the group.
    assert by_name["TST_0004"].product.startswith("IS3 family transposase")

    # Neither qualifier present is an empty string, never None.
    assert by_name["repB"].locus_tag == ""
    assert by_name["repB"].note == ""


def test_gene_label_shortens_locus_tags_only() -> None:
    """A tag is drawn by its number; a real name is left alone."""
    assert gene_label("SS37A_42600") == "42600"
    assert gene_label("SS37A_41660") == "41660"
    assert gene_label("TST_0002") == "0002"
    assert gene_label("ECDH10B-1234") == "1234"

    # Real gene names, and anything else not shaped like a tag.
    assert gene_label("recA") == "recA"
    assert gene_label("thrA") == "thrA"
    assert gene_label("b0001") == "b0001"          # no separator
    assert gene_label("sup_1234_like") == "sup_1234_like"   # digits not at the end
    assert gene_label("_42600") == "_42600"        # no letter in the prefix
    assert gene_label("") == ""


def test_gene_label_keeps_a_numbered_gene_name() -> None:
    """`ftsH_5` is the fifth ftsH, not tag number 5.

    Both shapes are <prefix><sep><digits>, so shape alone cannot tell them
    apart -- only how the number is written.  A locus tag is zero-padded to a
    fixed width across the replicon; a paralog suffix is a small integer.
    """
    assert gene_label("ftsH_5") == "ftsH_5"
    assert gene_label("rpoB_2") == "rpoB_2"
    assert gene_label("ynfM_12") == "ynfM_12"
    assert gene_label("tuf_1") == "tuf_1"

    # The boundary, stated explicitly: three digits is not yet a tag, four is.
    assert gene_label("ABC_999") == "ABC_999"
    assert gene_label("ABC_0999") == "0999"


def test_gene_label_keeps_the_rna_series_letter() -> None:
    """`SS37A_r00010` is the same tag prefix's RNA series -- keep the `r`.

    Dropping it would leave an rRNA indistinguishable from the protein-coding
    gene at the same index.  The width test still counts digits only, so a
    letter cannot smuggle a short number past it.
    """
    assert gene_label("SS37A_r00010") == "r00010"
    assert gene_label("SS37A_t00120") == "t00120"

    # The letter is not padding: three digits behind it is still not a tag.
    assert gene_label("ABC_r999") == "ABC_r999"
    assert gene_label("ftsH_a5") == "ftsH_a5"


def test_btop_grouping_is_by_real_read_name(fixture_dir: Path) -> None:
    stats = BtopStats()
    groups = list(iter_read_groups(fixture_dir / "reads.btop", 0, stats))
    sizes = {hits[0].read: len(hits) for hits in groups}

    assert sizes["read.3"] == 2
    assert sizes["read.12"] == 3
    assert sizes["0a1b2c3d-e4f5-6789-abcd-ef0123456789"] == 1
    assert stats.malformed == 0
    assert stats.out_of_order == 0


def test_btop_optional_trace_column() -> None:
    """The tool's own usage text says the trace column is optional."""
    hit = parse_line("r\t1\t1\t100\t500\t600\t1000\tctgA")
    assert hit.btop == ""
    assert hit.sspan == 100


def test_btop_rejects_wrong_outfmt() -> None:
    with pytest.raises(ValueError):
        parse_line("r\tplus\t1\t100\t500\t600\t1000\tctgA\t100")


# Byte for byte what BLAST+ 2.13.0 wrote for `-outfmt '6 delim=\t qseqid ...'`,
# the command the old --help text gave: the escape is not interpreted, so the
# "separator" is a backslash and a t.
_LITERAL_BACKSLASH_T = "probe.1\\t1\\t1\\t3000\\t1001\\t4000\\t3000\\tCTG1.1\\t3000\n"


def test_documented_outfmt_has_no_delim() -> None:
    """The one -outfmt string ReadRift prints must not make BLAST write '\\t'."""
    from readrift.cli import build_parser
    from readrift.inputs.btop import OUTFMT

    assert "delim" not in OUTFMT
    assert OUTFMT.split()[1:] == [
        "qseqid", "sframe", "qstart", "qend", "sstart", "send", "qlen", "sseqid", "btop",
    ]
    assert OUTFMT in build_parser().epilog


def test_nothing_parsed_explains_literal_backslash_t(tmp_path: Path) -> None:
    path = tmp_path / "delim.btop"
    path.write_text(_LITERAL_BACKSLASH_T * 3, encoding="utf-8")
    stats = BtopStats()
    assert list(iter_read_groups(path, 0, stats)) == []

    assert stats.nothing_parsed
    why = stats.explain_nothing_parsed(path)
    assert "delim=\\t" in why
    assert "Leave delim= out" in why


def test_nothing_parsed_explains_missing_outfmt(tmp_path: Path) -> None:
    path = tmp_path / "pairwise.txt"
    path.write_text("BLASTN 2.13.0+\n\n\nReference: Zheng Zhang ...\n", encoding="utf-8")
    stats = BtopStats()
    list(iter_read_groups(path, 0, stats))

    assert stats.nothing_parsed
    assert "without -outfmt" in stats.explain_nothing_parsed(path)


def test_short_reads_are_not_nothing_parsed(tmp_path: Path) -> None:
    """Every line dropped by --min-read-length still *parsed*: not a format error."""
    path = tmp_path / "short.btop"
    path.write_text("r\t1\t1\t100\t500\t600\t1000\tctgA\n", encoding="utf-8")
    stats = BtopStats()
    list(iter_read_groups(path, 10_000, stats))

    assert stats.short_reads == 1
    assert not stats.nothing_parsed


def test_btop_trace_counts() -> None:
    trace = parse_btop("120AG45-T7A-")
    assert trace.matches == 172
    assert trace.mismatches == 1
    assert trace.query_gaps == 1
    assert trace.subject_gaps == 1
    assert trace.identity == pytest.approx(172 / 175)


@pytest.mark.parametrize("name", ["reads.fastq", "reads.fastq.gz"])
def test_read_extraction_handles_plain_and_gzip(fixture_dir: Path, name: str) -> None:
    """Findings B07/B08: the Perl lost the first record of any .gz file."""
    records = list(iter_records(fixture_dir / name))
    ids = [read_id for read_id, _text in records]

    assert "read.1" in ids
    assert len(ids) == len(set(ids))
    # The first record must be present -- this is exactly what seek() on a
    # pipe silently destroyed.
    assert records[0][1].startswith("@")
