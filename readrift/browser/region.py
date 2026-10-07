"""One visible window, packed into lanes.

:class:`RegionView` is the single intermediate between the store and anything
that draws: the JSON the canvas front end consumes, and the matplotlib figure
:mod:`readrift.browser.export` writes.  Both read the same object, so a saved
image cannot disagree with what was on screen.

Lane assignment uses :class:`readrift.layout.LanePacker` -- the same first-fit
rule, applied to the same sort order, as the printed map.  A region therefore
looks like the corresponding slice of the poster rather than a second opinion
about it.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from readrift.browser.store import CLASS_NAMES, EVENT_KINDS, BrowserStore
from readrift.labels import gene_label
from readrift.layout import LanePacker

#: Default cap on reads returned for one window.  Above this the pileup is a
#: solid block anyway, and the browser has to draw every one of them.
DEFAULT_LIMIT = 4_000

#: Lanes available per strand.  Generous: the front end scrolls.
MAX_LANES = 4_000

#: Horizontal room reserved after each read, as a fraction of the span the
#: reader is actually looking at.
#:
#: Unlike the printed map, which reserves ``--space`` for the read label it
#: writes after every read, the browser draws no labels in the pileup.  All
#: this has to do is keep two reads sharing a lane from meeting and reading as
#: one long read, which is a handful of pixels -- so the fraction is chosen to
#: be about 3 px on a typical canvas and nothing more.  Every pixel beyond that
#: is a read that could have shared the lane and did not.
LANE_GAP_FRACTION = 0.0025

#: Annotation rows per strand, matching the printed map.
ANNOTATION_ROWS = 3

#: Lane order by read class -- long-divided (orange) against the coordinate
#: axis, then short-divided (blue), then undivided (green).
#:
#: The same order ``layout._requests`` uses for the printed map, and for the
#: same reason: the reads that carry structural evidence belong in the lanes
#: nearest the axis, where they are read first, rather than somewhere in the
#: depth of a pileup the user has to scroll through.  Packing is still
#: first-fit within that order, so a later class still takes a low lane when
#: one is horizontally free.
CLASS_LANE_ORDER: dict[str, int] = {
    "long_divided": 0,
    "short_divided": 1,
    "undivided": 2,
}

#: How close to a contig end a read must reach to earn a terminal arrow.
#:
#: The printed map uses ``layout._EDGE_TOLERANCE``, 10 *drawing units*, which is
#: 1000 bases at the default ``--map-scale`` and something else at any other.
#: A browser has no map scale, so this is stated in bases directly.  For reads
#: tens of kilobases long, landing within 100 bp of the end is "at the end".
CONTIG_EDGE_TOLERANCE = 100


@dataclass(slots=True)
class PackedSegment:
    """One aligned piece of a read, in reference coordinates."""

    start: int
    end: int
    strand: int
    inverted: bool
    """Runs opposite to the read's first HSP."""


@dataclass(slots=True)
class PackedMarker:
    """One decoration on a read: an arrowhead or a junction tick.

    The styles are the printed map's, from ``render/theme.arrow_points``:

    ``arrow_fwd`` / ``arrow_rev``
        the read runs on past this point -- off the edge of the window, or off
        the end of the contig.  Which of the two says which way it runs.
    ``arrow_ct``
        a double chevron: the read continues onto another contig.
    ``arrow_circular``
        a chevron and a ring: the read bridges the origin of a circular
        replicon.
    ``inversion_arrow``
        this piece of the read runs opposite to the rest of it.
    ``tick``
        a junction inside a short-divided read.
    """

    x: int
    """Reference coordinate."""

    style: str
    direction: int
    """``+1`` points right, ``-1`` left, ``0`` for a tick."""

    at_edge: bool = False
    """Placed at the window edge rather than at a real feature."""

    def as_json(self) -> list:
        # A list, not an object: a dense window carries thousands of these.
        return [self.x, self.style, self.direction, 1 if self.at_edge else 0]


@dataclass(slots=True)
class PackedRead:
    """One *line* of a read, placed in a lane, ready to draw.

    Usually one per read.  A long-divided read is emitted once per piece --
    see :func:`_placements` -- so several entries can share an :attr:`id`, and
    anything counting reads has to count distinct ids rather than entries.
    """

    id: int
    label: str
    cls: str
    junction: str
    inverted: bool
    strand: int
    lane: int
    start: int
    """Leftmost reference coordinate on this contig -- may be left of the view."""

    end: int
    segments: list[PackedSegment] = field(default_factory=list)

    hits: int = 1
    """Pieces the *whole read* has on this contig, across all its placements.

    Carried on every line so a tooltip on one of them can say what the read
    does, not just what this line shows."""

    joins: list[str] = field(default_factory=list)
    """Other contigs this read also maps to, for a ``Contig_Join``."""

    markers: list[PackedMarker] = field(default_factory=list)

    cont: str = "arrow_fwd"
    """Which arrow marks this read running on past the drawn range."""

    cont_dir: int = 0
    """``+1`` / ``-1`` when :attr:`cont` fixes the direction, ``0`` when the
    caller should point it outwards from whichever end it is drawn at.

    Sent because the client places its own edge arrows -- only the browser
    knows where the viewport ends, and the window fetched from the server is
    three times wider than it.  Sending the resolved style keeps the rule in
    one place instead of restating it in JavaScript."""

    def as_json(self) -> dict:
        return {
            "id": self.id,
            "label": self.label,
            "cls": self.cls,
            "junction": self.junction,
            "inverted": self.inverted,
            "strand": self.strand,
            "lane": self.lane,
            "start": self.start,
            "end": self.end,
            "segments": [
                [s.start, s.end, 1 if s.inverted else 0] for s in self.segments
            ],
            "hits": self.hits,
            "joins": self.joins,
            "markers": [m.as_json() for m in self.markers],
            "cont": self.cont,
            "cont_dir": self.cont_dir,
        }


#: Longest ``/product`` carried in a region payload, for the hover tooltip.
#:
#: A whole contig on screen is every feature it has, and a product runs to a
#: sentence.  Enough to recognise a gene by, with the full text one click away
#: from ``/api/annotation``.
PRODUCT_PREVIEW = 140


@dataclass(slots=True)
class PackedAnnotation:
    id: int
    """Index into the store's annotation arrays -- what ``/api/annotation`` takes."""

    start: int
    end: int
    strand: int
    row: int
    group: str
    feature: str
    name: str
    """Verbatim, as GenBank wrote it."""

    label: str
    """What to draw: :func:`readrift.labels.gene_label` of *name*."""

    product: str = ""
    """``/product``, truncated to :data:`PRODUCT_PREVIEW` for the tooltip."""

    def as_json(self) -> dict:
        return {
            "id": self.id,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "row": self.row,
            "group": self.group,
            "feature": self.feature,
            "name": self.name,
            "label": self.label,
            "product": self.product,
        }


@dataclass(slots=True)
class RegionView:
    """Everything visible in one window."""

    contig: str
    contig_length: int
    start: int
    end: int
    reads: list[PackedRead] = field(default_factory=list)
    """One entry per drawn *line*, not per read -- see :func:`_placements`."""

    annotations: list[PackedAnnotation] = field(default_factory=list)
    events: list[dict] = field(default_factory=list)
    counts: dict[str, int] = field(default_factory=dict)
    """Reads of each class overlapping the window, before the density cap."""

    shown: dict[str, int] = field(default_factory=dict)
    """Reads of each class actually returned."""

    lanes: dict[str, int] = field(default_factory=dict)
    """Highest lane used, per strand (``"1"`` / ``"-1"``)."""

    truncated: bool = False

    def as_json(self) -> dict:
        return {
            "contig": self.contig,
            "contig_length": self.contig_length,
            "start": self.start,
            "end": self.end,
            "reads": [r.as_json() for r in self.reads],
            "annotations": [a.as_json() for a in self.annotations],
            "events": self.events,
            "counts": self.counts,
            "shown": self.shown,
            "lanes": self.lanes,
            "truncated": self.truncated,
        }


# --------------------------------------------------------------------------
# Density cap
# --------------------------------------------------------------------------


def _priority(cls: str, junction: str, inverted: bool) -> int:
    """Lower sorts first.  Structural evidence outranks background pileup.

    When a window holds more reads than can be drawn, dropping them at random
    would hide exactly what the program exists to show.  Long-divided reads,
    inversions and joins are kept before anything else; undivided reads are the
    ones that get sampled, and the caller is told how many were left out.
    """
    if cls == "long_divided":
        return 0
    if inverted:
        return 1
    if junction != "C":
        return 2
    if cls == "short_divided":
        return 3
    return 4


def drawn_segment_minimum(store: BrowserStore) -> int:
    """``--min-match-length`` from the run that built *store*, or 0.

    The printed map does not draw the short anchors either side of a
    long-distance junction: :func:`readrift.classify.drawable_hits` drops HSPs
    with ``sspan <= --min-match-length`` from a ``LONG_DIVIDED`` read so the
    map is not cluttered by the fragments a structural event leaves behind.
    Applying the same cut here is what keeps a browser region the *same
    picture* as the matching slice of the poster.

    It is a drawing rule and only a drawing rule.  Those bases are still
    matched bases: they are in the coverage profile, in the statistics and in
    the read's alignment table under ``/api/read``, exactly as in the report
    (decision D3).  So the filter lives here, in the region builder, and never
    in the cache.

    ``0`` when the cache predates the parameter being recorded, which keeps
    every segment -- the old behaviour.
    """
    try:
        return int(store.meta.get("params", {}).get("min_match_length", 0) or 0)
    except (TypeError, ValueError):  # pragma: no cover - hand-edited meta
        return 0


def continuation_style(junction: str, strand: int) -> str:
    """Which arrow marks a read running on past the point it is drawn to.

    Identical to :func:`readrift.layout._continuation_style`, which decides the
    same thing for the printed map -- a circular join outranks a contig join,
    which outranks the plain direction arrow.
    """
    if junction == "Circle":
        return "arrow_circular"
    if junction == "Contig_Join":
        return "arrow_ct"
    return "arrow_fwd" if strand > 0 else "arrow_rev"


def continuation_direction(style: str) -> int:
    """Which way *style* points, or ``0`` if it depends on where it is drawn.

    ``arrow_fwd`` and ``arrow_rev`` encode the strand, and the printed map
    draws the same one at both ends of a read -- the chevron says which way the
    molecule runs, not which end you are looking at.  The junction arrows carry
    no strand, so they point outwards instead.
    """
    if style == "arrow_fwd":
        return +1
    if style == "arrow_rev":
        return -1
    return 0


def _markers_for(
    packed: PackedRead,
    contig_length: int,
    view_start: int,
    view_end: int,
) -> list[PackedMarker]:
    """Every decoration one read earns.

    The printed map draws its continuation arrows at page folds; a browser has
    no folds, so the window edge takes their place -- the meaning is the same,
    "this line carries on past here".  Contig-end, circular, contig-join and
    inversion marks follow the map's rules exactly.
    """
    style = packed.cont
    fixed = packed.cont_dir
    markers: list[PackedMarker] = []

    def arrow(x: int, side: int, at_edge: bool = False) -> PackedMarker:
        return PackedMarker(x, style, fixed or side, at_edge=at_edge)

    # Runs on past the left or right edge of what is being drawn.  These are
    # tagged at_edge because they are only right for a caller drawing exactly
    # this range -- the image exporter.  The interactive client discards them
    # and places its own, since the window it fetched is wider than the one it
    # shows.
    if packed.start < view_start:
        markers.append(arrow(view_start, -1, at_edge=True))
    if packed.end > view_end:
        markers.append(arrow(view_end, +1, at_edge=True))

    # Reaches the very beginning or the very end of the contig.  For a Circle
    # read this fires at both ends, which is exactly what makes the join
    # readable: the same read leaves one end and arrives at the other.
    if packed.start <= CONTIG_EDGE_TOLERANCE:
        markers.append(arrow(packed.start, -1))
    if contig_length - packed.end <= CONTIG_EDGE_TOLERANCE:
        markers.append(arrow(packed.end, +1))

    # A piece running against the rest of the read.  As on the printed map
    # (layout.py), the marker is the piece's left end: the arrow comes in
    # diagonally from above-left and stops short of it, so it points at the
    # inverted piece without covering it -- a short inversion may be only a
    # few pixels wide.
    for segment in packed.segments:
        if segment.inverted:
            markers.append(PackedMarker(segment.start, "inversion_arrow", +1))

    # Junctions inside a short-divided read, drawn as several pieces on one lane.
    if packed.cls == "short_divided" and len(packed.segments) > 1:
        for segment in packed.segments:
            markers.append(PackedMarker(segment.start, "tick", 0))
            markers.append(PackedMarker(segment.end, "tick", 0))

    return markers


def _placements(read: PackedRead) -> list[PackedRead]:
    """Split *read* into the lines it should occupy -- the printed map's rule.

    ``layout._requests`` gives a short-divided read **one** placement covering
    all its pieces, because a small indel is one event and its pieces belong on
    one line; every other class gets **one placement per HSP**, so the two ends
    of a structural event are independent lines carrying the same label.

    The browser used to place every read once, across the hull of its pieces.
    For an ordinary long-divided read that wastes the lane between the two
    ends; for a read that bridges the origin of a circular replicon it is
    ruinous.  Such a read has a piece at each end of the contig::

        AP027142        1-    20,357  -   read      1- 20,079
        AP027142  3,729,139- 3,733,940  -   read 20,087- 24,835

    so its hull is the **whole replicon**: one read reserved a lane across
    3.7 Mb and was drawn joined by a connector spanning reference it never
    crosses -- it goes the other way, around the origin.  Placed per piece,
    each end reserves only itself, the lane between them is free for other
    reads, and no connector is drawn because each line has a single segment.
    The circular arrow at each end still says what happened; that is what the
    arrow vocabulary is for.
    """
    if read.cls == "short_divided" or len(read.segments) < 2:
        return [read]
    return [
        replace(read, segments=[segment], start=segment.start, end=segment.end)
        for segment in read.segments
    ]


def _apply_limit(rows: list[tuple], limit: int) -> tuple[list[tuple], bool]:
    """Keep the most informative *limit* rows, evenly sampling the remainder."""
    if len(rows) <= limit:
        return rows, False

    ranked = sorted(range(len(rows)), key=lambda i: rows[i][0])
    keep_ranked: list[int] = []
    overflow: list[int] = []
    for position, index in enumerate(ranked):
        if rows[index][0] <= 2 and position < limit:
            keep_ranked.append(index)
        else:
            overflow.append(index)

    room = limit - len(keep_ranked)
    if room > 0 and overflow:
        # Sample the background reads evenly *along the contig*, not in rank
        # order, so thinning the pileout does not carve a fake coverage hole
        # into one end of the window.
        overflow.sort(key=lambda i: rows[i][1])
        step = len(overflow) / room
        keep_ranked.extend(overflow[int(i * step)] for i in range(room))

    chosen = set(keep_ranked)
    return [row for i, row in enumerate(rows) if i in chosen], True


# --------------------------------------------------------------------------
# Building
# --------------------------------------------------------------------------


def build_region(
    store: BrowserStore,
    contig: str,
    start: int,
    end: int,
    classes: set[str] | None = None,
    inverted_only: bool = False,
    limit: int = DEFAULT_LIMIT,
    with_events: bool = True,
    view_span: int | None = None,
) -> RegionView:
    """Collect and pack everything visible on *contig* between *start* and *end*.

    *view_span* is how many bases the caller will actually show, when that is
    narrower than ``start..end``.  The interactive client fetches three screens
    around the view so panning does not hit the network, and the lane gap is a
    pixel quantity expressed as a fraction of a span -- so without this it would
    be computed from three times the span on screen and come out three times
    too wide, wasting lanes, and would not match what ``export.py`` gets by
    asking for exactly the range it draws.  Defaults to the fetched range.
    """
    info = store.require_contig(contig)
    start = max(0, int(start))
    end = min(info.length, int(end))
    if end <= start:
        end = min(info.length, start + 1)

    view = RegionView(
        contig=info.name,
        contig_length=info.length,
        start=start,
        end=end,
        counts={name: 0 for name in CLASS_NAMES},
        shown={name: 0 for name in CLASS_NAMES},
    )

    # Annotations and events first: a stretch with no reads over it still has
    # genes under it, and that is exactly the stretch worth looking at.
    view.lanes = {"1": 0, "-1": 0}
    view.annotations = _pack_annotations(store, info.name, start, end)
    if with_events:
        view.events = store.events_in(info.name, start, end, EVENT_KINDS)

    d = store.column
    hit_indices = store.hits_in(info.name, start, end)
    if hit_indices.size == 0:
        return view

    # One read may contribute several hits to the window; pack it once.
    read_ids = np.unique(d("hit_read")[hit_indices])

    read_cls = d("read_cls")
    read_junction = d("read_junction")
    read_inverted = d("read_inverted")
    read_strand = d("read_strand")
    hit_contig = d("hit_contig")
    hit_slow = d("hit_slow")
    hit_shigh = d("hit_shigh")
    hit_strand = d("hit_strand")

    # None means "no filter given"; an empty set means "every class unticked".
    wanted = set(CLASS_NAMES) if classes is None else classes

    min_drawn = drawn_segment_minimum(store)

    rows: list[tuple] = []
    for read_id in read_ids.tolist():
        cls = CLASS_NAMES[int(read_cls[read_id])]
        junction = store.junction_name(read_junction[read_id])
        is_inverted = bool(read_inverted[read_id])

        # The read's whole footprint on this contig, not just the part in view:
        # a lane must stay reserved across the pieces of a divided read, and the
        # front end needs to know the line continues past the edge.
        own = store.read_hits_on(read_id, info.index)
        if cls == "long_divided" and min_drawn > 0:
            own = own[(hit_shigh[own] - hit_slow[own]) > min_drawn]

        # Nothing left to draw for this read here.  Counted before the class
        # filter but after this one, because `counts` answers "is there
        # anything here the filter is hiding from me" -- and a read with no
        # drawable segment is not being hidden by the filter.
        if own.size == 0:
            continue

        view.counts[cls] += 1
        if cls not in wanted:
            continue
        if inverted_only and not is_inverted:
            continue

        # Measured against the read's first HSP in BTOP order (`read_strand`),
        # which is the same reference `ReadGroup.inverted` uses -- so the
        # segments drawn in the inversion colour are exactly the ones that made
        # the read count as inverted.
        #
        # Only a short-distance inversion is drawn as one, though: thick, red
        # and arrowed.  It is the mark of phase variation (Klemm 1986, PMID
        # 2874022), typically a few hundred bases that would vanish otherwise.
        # A long-divided read's pieces are separate lines anyway, already
        # above and below the axis by strand, and the printed map lays them
        # out one piece at a time, so it never marks them either.  The read
        # still counts as inverted -- badge, tooltip and the inverted-only
        # filter go by `is_inverted`, not by its segments.
        strand = int(read_strand[read_id])
        emphasise = cls != "long_divided"
        segments = [
            PackedSegment(
                start=int(hit_slow[h]),
                end=int(hit_shigh[h]),
                strand=int(hit_strand[h]),
                inverted=emphasise and int(hit_strand[h]) != strand,
            )
            for h in own.tolist()
        ]
        segments.sort(key=lambda s: s.start)

        joins: list[str] = []
        if junction == "Contig_Join":
            all_hits = store.read_hits(read_id)
            joins = sorted(
                {
                    store.contigs[int(hit_contig[h])].name
                    for h in all_hits.tolist()
                    if int(hit_contig[h]) != info.index
                }
            )

        cont = continuation_style(junction, strand)
        packed = PackedRead(
            id=int(read_id),
            label=store.labels[int(read_id)],
            cls=cls,
            junction=junction,
            inverted=is_inverted,
            cont=cont,
            cont_dir=continuation_direction(cont),
            strand=strand,
            lane=0,
            start=segments[0].start,
            end=segments[-1].end,
            segments=segments,
            hits=len(segments),
            joins=joins,
        )
        rows.append((_priority(cls, junction, is_inverted), packed.start, packed))

    kept, truncated = _apply_limit(rows, max(1, limit))
    view.truncated = truncated

    # Split into drawn lines first, then sort *those*.  `shown` counts reads;
    # `view.reads` carries one entry per line.
    lines: list[PackedRead] = []
    for _priority_value, _start, packed in kept:
        view.shown[packed.cls] += 1
        lines.extend(_placements(packed))

    # Pack by class, then in ascending start order, one packer per strand --
    # `layout.build_layout`'s rule applied to `layout._requests`' output, so the
    # same read lands in the same relative place as on the poster.  See
    # CLASS_LANE_ORDER for why the class comes first.
    #
    # Sorting by *line* and not by read is what makes the split worth anything.
    # LanePacker keeps one "last occupied x" per lane and cannot represent a
    # hole, so packing a read's far piece straight after its near one drags that
    # lane's end across the gap between them and blocks every read in between
    # from a lane that is in fact free -- exactly the space the split was meant
    # to release.  Sorted by line, the far piece is packed among the reads it
    # actually sits next to.
    lines.sort(
        key=lambda line: (
            CLASS_LANE_ORDER.get(line.cls, len(CLASS_LANE_ORDER)),
            line.start,
            line.end,
        )
    )

    shown_span = end - start if not view_span or view_span <= 0 else view_span
    gap = max(1.0, min(shown_span, end - start) * LANE_GAP_FRACTION)
    packers = {1: LanePacker(MAX_LANES), -1: LanePacker(MAX_LANES)}
    highest = {1: 0, -1: 0}

    for line in lines:
        key = 1 if line.strand > 0 else -1
        lane = packers[key].place(float(line.start), float(line.end) + gap)
        line.lane = lane
        line.markers = _markers_for(line, info.length, start, end)
        highest[key] = max(highest[key], lane)
        view.reads.append(line)

    view.lanes = {"1": highest[1], "-1": highest[-1]}
    return view


def _pack_annotations(
    store: BrowserStore, contig: str, start: int, end: int
) -> list[PackedAnnotation]:
    """Assign annotations to stacked rows so neighbouring labels stay legible.

    The printed map advances a row counter modulo three; here the rows are
    packed first-fit instead, which is the same idea done properly -- an
    annotation only moves down when the row above is genuinely occupied, so a
    sparse region stays on one line and a dense one fans out.

    "Occupied" includes room for the label, which is a pixel quantity, so the
    reserved gap is a fraction of the visible span rather than a fixed number
    of bases.
    """
    indices = store.annotations_in(contig, start, end)
    if indices.size == 0:
        return []

    d = store.column
    groups = store.meta.get("anno_groups", [])
    features = store.meta.get("anno_features", [])
    label_gap = max(1.0, (end - start) * 0.02)

    row_end: dict[int, list[float]] = {
        1: [float("-inf")] * ANNOTATION_ROWS,
        -1: [float("-inf")] * ANNOTATION_ROWS,
    }
    out: list[PackedAnnotation] = []

    for index in indices.tolist():
        begin = int(d("anno_start")[index])
        finish = int(d("anno_end")[index])
        strand = 1 if int(d("anno_strand")[index]) >= 0 else -1

        ends = row_end[strand]
        assigned = next(
            (r for r in range(ANNOTATION_ROWS) if ends[r] < begin),
            min(range(ANNOTATION_ROWS), key=lambda r: ends[r]),
        )
        ends[assigned] = max(ends[assigned], finish + label_gap)

        gid = int(d("anno_group")[index])
        fid = int(d("anno_feature")[index])
        name = store.anno_names[index]
        product = store.anno_products[index]
        if len(product) > PRODUCT_PREVIEW:
            product = product[:PRODUCT_PREVIEW].rstrip() + "…"
        out.append(
            PackedAnnotation(
                id=index,
                start=begin,
                end=finish,
                strand=strand,
                row=assigned,
                group=groups[gid] if gid < len(groups) else "",
                feature=features[fid] if fid < len(features) else "",
                name=name,
                label=gene_label(name),
                product=product,
            )
        )

    return out
