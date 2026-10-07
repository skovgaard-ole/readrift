/* readrift browser -- interactive view over one cached run.
 *
 * Plain ES2020, no framework, no bundler, nothing loaded from a network.
 *
 * Two canvases: a fixed header (ruler, coverage, annotations, event marks) and
 * a scrollable pileup below it.  Both share one horizontal coordinate system
 * -- `state.start` .. `state.end` in reference bases -- so panning and zooming
 * only ever change two numbers.
 *
 * Colours come from /api/meta, which reads them from readrift/render/theme.py.
 * Nothing meaningful is hard-coded here; see style.css for why.
 */

'use strict';

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

// The view is laid out like the printed map: forward-strand reads above the
// coordinate axis, reverse-strand reads below it, mirrored about the axis.
// The header carries everything that belongs to the reference itself and is
// the only fixed row; the two read panes scroll away from the axis.
//
//        + strand reads      lane 1 against the axis, growing upward
//        coverage            standing on the axis
//   ==== coordinate axis ====
//        coordinates, event marks
//        annotations
//        − strand reads      lane 1 against the axis, growing downward

const COVERAGE_H = 74;      // coverage band; its baseline *is* the axis
const AXIS_Y = COVERAGE_H;  // the coordinate axis, in header-canvas pixels
const RULER_H = 28;         // ticks, event marks and coordinates, under it

// Structural-event marks: a triangle under the axis, apex on it, pointing back
// up at the coordinate it marks.  Stated once because drawEventMarks draws
// them and eventAt has to be able to find them again under the pointer.
const EVENT_MARK_TOP = AXIS_Y + 1;
const EVENT_MARK_H = 7;
const EVENT_MARK_W = 4;     // half-width at the base
const ANNO_ROW_H = 13;
const ANNO_ROWS = 3;
const ANNO_H = ANNO_ROWS * 2 * ANNO_ROW_H + 14;
const HEADER_H = COVERAGE_H + RULER_H + ANNO_H;

const PILEUP_PAD = 10;
const MIN_LANE = 1.5;
const MAX_LANE = 11;
const MAX_CANVAS_H = 16000;

const MIN_SPAN = 60;          // never zoom past ~60 bases across the window
const FETCH_FACTOR = 3;       // fetch this many screens' worth around the view
const REFETCH_ZOOM = 4;       // refetch once the view is this much smaller
const DEBOUNCE_MS = 60;

// Arrowheads.  The symbols are the printed map's -- single chevron for a read
// running on, double for a contig join, chevron plus ring for a circular join,
// a chevron against the run for an inverted piece, a bar for a junction -- but
// the proportions are the browser's own. The PostScript shapes in theme.py are
// tuned for a 10 750 pt sheet with a 20 pt shaft; scaled into a lane a few
// pixels tall, that shaft would swamp the read and the head would vanish.
// browser/export.py draws the real theme.py geometry, so an exported figure is
// still literally the report's arrows.
const ARROW_MIN_LANE = 4.5;      // below this the head is illegible, so skip it
const ARROW_LANE_FRACTION = 1.0; // head length, and height, as a share of a lane
const ARROW_MIN_PX = 4;
const ARROW_MAX_PX = 11;

// Decorations -- arrowheads, junction ticks, the circular ring -- are a line
// drawing, and their stroke has to stay thin whatever the read bars do.  Taking
// it from the bar width, as this used to, put a 5 px stroke on a 9 px chevron:
// the open V filled in solid and every arrow read as a blob.  The printed map
// has the same relation the other way round -- a 2 pt stroke on a 5 pt chevron
// in a 15 pt lane -- because its bars are a fraction of their lane, where the
// browser's fill most of it.  So the two widths cannot share a source.
const MARKER_STROKE_FRACTION = 0.26;   // of the arrowhead length
const MARKER_STROKE_MIN = 0.9;
const MARKER_STROKE_MAX = 1.8;

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

const state = {
  meta: null,
  contig: null,          // {name, length, ...}
  start: 0,
  end: 1,
  classes: new Set(),
  invertedOnly: false,
  region: null,          // last RegionView payload
  regionKey: null,       // what the cached region covers
  depth: null,
  selected: null,        // read id
  selectedAnno: null,    // annotation id
  laneH: MAX_LANE,
  hitboxes: { '1': [], '-1': [] },   // per pane, in that pane's own pixels
  annoBoxes: [],         // annotation targets, in header-canvas pixels
  busy: 0,
};

let regionAbort = null;
let depthAbort = null;
let refreshTimer = null;
let searchTimer = null;
let suggestions = [];
let suggestionIndex = -1;

const el = (id) => document.getElementById(id);

const headerCanvas = el('header-canvas');
const plusCanvas = el('plus-canvas');
const minusCanvas = el('minus-canvas');
const plusScroll = el('plus-scroll');
const minusScroll = el('minus-scroll');
const tooltip = el('tooltip');

const strandKey = (strand) => (strand > 0 ? '1' : '-1');
const canvasFor = (strand) => (strand > 0 ? plusCanvas : minusCanvas);
const scrollFor = (strand) => (strand > 0 ? plusScroll : minusScroll);

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function commas(n) {
  return Math.round(n).toLocaleString('en-US');
}

// With a step (ruler ticks), every label gets enough decimals to resolve that
// step, so 1,130,000 and 1,135,000 read 1.130 and 1.135 Mb, never 1.13 twice.
function shortBp(n, step) {
  const v = Math.abs(n);
  const unit = v >= 1e6 ? 1e6 : (v >= 1e3 ? 1e3 : 1);
  if (unit === 1) return commas(n) + ' bp';
  const digits = step === undefined
    ? (v % unit === 0 ? 0 : (unit === 1e6 ? 2 : 1))
    : Math.max(0, Math.ceil(-Math.log10(step / unit) - 1e-9));
  return (n / unit).toFixed(digits) + (unit === 1e6 ? ' Mb' : ' kb');
}

function clamp(v, lo, hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

function niceStep(span, target) {
  const raw = span / Math.max(1, target);
  const mag = Math.pow(10, Math.floor(Math.log10(Math.max(raw, 1))));
  for (const m of [1, 2, 5]) {
    if (raw <= mag * m) return mag * m;
  }
  return mag * 10;
}

function setBusy(delta) {
  state.busy = Math.max(0, state.busy + delta);
  el('status-busy').hidden = state.busy === 0;
}

async function getJSON(url, signal) {
  const response = await fetch(url, { signal });
  const payload = await response.json();
  if (!response.ok) throw new Error(payload.error || response.statusText);
  return payload;
}

// ---------------------------------------------------------------------------
// Canvas sizing
// ---------------------------------------------------------------------------

// All three canvases must map the same base to the same pixel, so the ruler
// lines up with the reads above and below it.  The read panes can carry a
// scrollbar, so their width is the authority and the header is sized to match;
// both .strand-scroll panes reserve the scrollbar permanently (overflow-y:
// scroll) so this width never oscillates.
function plotWidth() {
  return plusCanvas.clientWidth || minusCanvas.clientWidth ||
         headerCanvas.clientWidth || 800;
}

function sizeCanvas(canvas, cssHeight) {
  const dpr = window.devicePixelRatio || 1;
  const cssWidth = plotWidth();
  if (canvas === headerCanvas) canvas.style.width = cssWidth + 'px';
  canvas.style.height = cssHeight + 'px';
  const w = Math.round(cssWidth * dpr);
  const h = Math.round(cssHeight * dpr);
  if (canvas.width !== w || canvas.height !== h) {
    canvas.width = w;
    canvas.height = h;
  }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssWidth, cssHeight);
  return ctx;
}

function xOf(position) {
  const span = state.end - state.start || 1;
  return (position - state.start) / span * plotWidth();
}

function positionOf(x) {
  const span = state.end - state.start || 1;
  return state.start + x / plotWidth() * span;
}

// ---------------------------------------------------------------------------
// Drawing: header
// ---------------------------------------------------------------------------

function drawHeader() {
  const width = plotWidth();
  const ctx = sizeCanvas(headerCanvas, HEADER_H);
  const palette = state.meta.palette;

  ctx.fillStyle = palette.surface;
  ctx.fillRect(0, 0, width, HEADER_H);

  // Coverage first: the ruler draws the axis line it stands on, and the axis
  // has to land on top of the profile rather than under it.
  drawCoverage(ctx, width, palette);
  drawRuler(ctx, width, palette);
  drawAnnotations(ctx, width, palette);
  drawEventMarks(ctx, width, palette);
}

// The coordinate axis, with its ticks and coordinates hanging below it -- the
// reference's own line, with forward reads and the coverage profile above and
// reverse reads and the annotation tracks below.
function drawRuler(ctx, width, palette) {
  const span = state.end - state.start;
  const step = niceStep(span, Math.max(4, Math.floor(width / 130)));
  const minor = step / 5;
  const axis = AXIS_Y + 0.5;

  ctx.strokeStyle = palette.axis;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(0, axis);
  ctx.lineTo(width, axis);
  ctx.stroke();

  ctx.beginPath();
  ctx.strokeStyle = palette.grid;
  for (let p = Math.ceil(state.start / minor) * minor; p <= state.end; p += minor) {
    const x = Math.round(xOf(p)) + 0.5;
    ctx.moveTo(x, axis);
    ctx.lineTo(x, axis + 5);
  }
  ctx.stroke();

  ctx.font = '11px ui-monospace, Consolas, monospace';
  ctx.textBaseline = 'alphabetic';
  ctx.fillStyle = palette.ink_secondary;
  ctx.strokeStyle = palette.axis;
  ctx.beginPath();
  for (let p = Math.ceil(state.start / step) * step; p <= state.end; p += step) {
    const x = Math.round(xOf(p)) + 0.5;
    ctx.moveTo(x, axis);
    ctx.lineTo(x, axis + 9);
    const label = step >= 1000 ? shortBp(p, step) : commas(p);
    ctx.textAlign = x < 34 ? 'left' : (x > width - 34 ? 'right' : 'center');
    ctx.fillText(label, x, AXIS_Y + 24);
  }
  ctx.stroke();
}

// Sits directly on the coordinate axis, growing upward towards the forward
// reads: the axis line drawn by drawRuler is this profile's baseline, so
// coverage and the reference it was measured against share one zero.
function drawCoverage(ctx, width, palette) {
  const top = 6;
  const bottom = AXIS_Y;
  const height = bottom - top;

  const depth = state.depth;
  ctx.font = '10px ui-monospace, Consolas, monospace';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';

  if (!depth || !depth.positions.length) {
    ctx.fillStyle = palette.ink_muted;
    ctx.fillText('coverage', 4, top);
    return;
  }

  const peak = Math.max(1, ...depth.max);
  const scale = (v) => bottom - height * Math.min(1, v / peak);
  const colour = palette.classes.short_divided;

  // min-to-max band, then the mean.  The band is what keeps a narrow trough
  // visible when a whole contig is on screen.
  ctx.beginPath();
  for (let i = 0; i < depth.positions.length; i++) {
    const x = xOf(depth.positions[i]);
    if (i === 0) ctx.moveTo(x, scale(depth.max[i]));
    else ctx.lineTo(x, scale(depth.max[i]));
  }
  for (let i = depth.positions.length - 1; i >= 0; i--) {
    ctx.lineTo(xOf(depth.positions[i]), scale(depth.min[i]));
  }
  ctx.closePath();
  ctx.fillStyle = colour;
  ctx.globalAlpha = 0.22;
  ctx.fill();
  ctx.globalAlpha = 1;

  ctx.beginPath();
  for (let i = 0; i < depth.positions.length; i++) {
    const x = xOf(depth.positions[i]);
    const y = scale(depth.mean[i]);
    if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.strokeStyle = colour;
  ctx.lineWidth = 1.2;
  ctx.stroke();

  ctx.fillStyle = palette.ink_muted;
  ctx.fillText('coverage  0–' + commas(peak) + '×', 4, top);
}

function drawAnnotations(ctx, width, palette) {
  const centre = AXIS_Y + RULER_H + ANNO_H / 2;
  ctx.strokeStyle = palette.grid;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(0, centre + 0.5);
  ctx.lineTo(width, centre + 0.5);
  ctx.stroke();

  const annotations = state.region ? state.region.annotations : [];
  state.annoBoxes = [];
  if (!annotations.length) {
    ctx.font = '10px ui-monospace, Consolas, monospace';
    ctx.fillStyle = palette.ink_muted;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.fillText(
      state.meta.is_genbank ? 'no features here' : 'no annotations (FASTA reference)',
      4, centre - ANNO_ROW_H);
    return;
  }

  ctx.lineWidth = 4;
  ctx.font = '9px ui-monospace, Consolas, monospace';
  ctx.textBaseline = 'alphabetic';

  for (const a of annotations) {
    const x0 = xOf(a.start);
    const x1 = xOf(a.end);
    if (x1 < -40 || x0 > width + 40) continue;
    const colour = palette.annotations[a.group] !== undefined
      ? palette.annotations[a.group]
      : palette.annotations[''];
    const y = a.strand > 0
      ? centre - 8 - a.row * ANNO_ROW_H
      : centre + 8 + a.row * ANNO_ROW_H;

    const right = Math.max(x1, x0 + 1);
    ctx.strokeStyle = colour;
    ctx.beginPath();
    ctx.moveTo(x0, y);
    ctx.lineTo(right, y);
    ctx.stroke();

    // A bracket around the selected feature: the bar is 4 px of colour that
    // already means something, so the selection cannot be a colour change.
    if (a.id === state.selectedAnno) {
      ctx.save();
      ctx.strokeStyle = palette.ink;
      ctx.lineWidth = 1;
      ctx.setLineDash([3, 2]);
      ctx.strokeRect(x0 - 2, y - 5, (right - x0) + 4, 10);
      ctx.setLineDash([]);
      ctx.restore();
    }

    state.annoBoxes.push({ id: a.id, x0: x0, x1: right, y: y, anno: a });

    // `label`, not `name`: a locus tag is the same prefix on every feature of
    // the assembly, so only its number distinguishes this one.  The full name
    // is in the tooltip and the panel.
    const text = a.label || a.name;
    if (text && x1 - x0 > 26) {
      const symbol = palette.annotation_symbols[a.group] || '';
      ctx.fillStyle = colour;
      ctx.textAlign = 'center';
      const cx = clamp((x0 + x1) / 2, 26, width - 26);
      ctx.fillText(text + symbol, cx, a.strand > 0 ? y - 4 : y + 10);
    }
  }
}

// The feature bar under the pointer, if any.  The narrowest wins, so a small
// gene inside a long `repeat_region` is still reachable.
function annotationAt(x, y) {
  let best = null;
  for (const box of state.annoBoxes) {
    if (Math.abs(box.y - y) > 5) continue;
    if (x < box.x0 - 2 || x > box.x1 + 2) continue;
    if (!best || (box.x1 - box.x0) < (best.x1 - best.x0)) best = box;
  }
  return best;
}

// Below the axis, pointing back up at the coordinate they mark.  Drawn after
// the ruler so a mark covers the tick it lands on rather than the other way
// round -- an event is worth more than one grey tick.
//
// These are the run's structural evidence, one triangle per event recorded by
// stats._record_events: an orange junction where a read jumps to a distant
// position, a red inversion breakpoint, an indigo circular join, a black
// contig join.  Hover one for what it means; `n` and `p` tour them.
function drawEventMarks(ctx, width, palette) {
  const events = state.region ? state.region.events : [];
  if (!events || !events.length) return;
  for (const event of events) {
    const x = xOf(event.position);
    if (x < -6 || x > width + 6) continue;
    ctx.fillStyle = palette.events[event.kind] || palette.ink;
    ctx.beginPath();
    ctx.moveTo(x, EVENT_MARK_TOP);
    ctx.lineTo(x - EVENT_MARK_W, EVENT_MARK_TOP + EVENT_MARK_H);
    ctx.lineTo(x + EVENT_MARK_W, EVENT_MARK_TOP + EVENT_MARK_H);
    ctx.closePath();
    ctx.fill();
  }
}

// The event mark under the pointer, if any.  Generous by a couple of pixels in
// both directions: the triangle is 8 px wide and the reader is aiming at a
// coordinate, not at a shape.
function eventAt(x, y) {
  const events = state.region ? state.region.events : [];
  if (!events || !events.length) return null;
  if (y < EVENT_MARK_TOP - 2 || y > EVENT_MARK_TOP + EVENT_MARK_H + 2) return null;

  let best = null;
  let bestDistance = EVENT_MARK_W + 2;
  for (const event of events) {
    const distance = Math.abs(xOf(event.position) - x);
    if (distance <= bestDistance) { best = event; bestDistance = distance; }
  }
  return best;
}

// ---------------------------------------------------------------------------
// Drawing: pileup
// ---------------------------------------------------------------------------

// One lane height for both panes, taken from whichever pane is under more
// pressure: a read has to be the same thickness above and below the axis or
// the two halves stop reading as one picture.
function pileupLayout() {
  const region = state.region;
  const lanes = {
    1: region ? (region.lanes['1'] || 0) : 0,
    '-1': region ? (region.lanes['-1'] || 0) : 0,
  };
  const room = (strand) => Math.max(60, scrollFor(strand).clientHeight);
  const fit = (strand) =>
    (room(strand) - PILEUP_PAD * 2) / Math.max(1, lanes[strandKey(strand)]);
  const laneH = clamp(Math.min(fit(1), fit(-1)), MIN_LANE, MAX_LANE);

  const height = (strand) => Math.max(
    room(strand),
    Math.min(MAX_CANVAS_H,
             PILEUP_PAD * 2 + lanes[strandKey(strand)] * laneH));

  return { lanes, laneH, heights: { 1: height(1), '-1': height(-1) } };
}

function drawPileup() {
  const layout = pileupLayout();
  state.laneH = layout.laneH;
  drawStrand(1, layout);
  drawStrand(-1, layout);
}

// One pane: the forward reads above the axis or the reverse reads below it.
// Lane 1 is always the row against the axis, so the forward pane counts its
// lanes up from its own bottom edge and the reverse pane down from its top.
function drawStrand(strand, layout) {
  const width = plotWidth();
  const canvas = canvasFor(strand);
  const scroll = scrollFor(strand);
  const height = layout.heights[strandKey(strand)];
  const palette = state.meta.palette;

  // Keep the forward pane against the axis.  Only re-pin when it was already
  // there, so scrolling up into the deep lanes is not undone by the next pan.
  const pinned = strand > 0 &&
    scroll.scrollHeight - scroll.scrollTop - scroll.clientHeight <= 2;
  const restore = () => {
    if (pinned) scroll.scrollTop = scroll.scrollHeight;
  };

  const ctx = sizeCanvas(canvas, height);
  ctx.fillStyle = palette.surface;
  ctx.fillRect(0, 0, width, height);

  const boxes = [];
  state.hitboxes[strandKey(strand)] = boxes;

  const yOf = strand > 0
    ? (lane) => height - PILEUP_PAD - (lane - 0.5) * layout.laneH
    : (lane) => PILEUP_PAD + (lane - 0.5) * layout.laneH;

  const region = state.region;
  const reads = region
    ? region.reads.filter((read) => (read.strand > 0) === (strand > 0))
    : [];

  if (!reads.length) {
    // The "nothing here" message belongs to the region, not to one strand, so
    // it is written once -- in the forward pane, where the eye starts.
    if (strand > 0 && (!region || !region.reads.length)) {
      let message = 'no reads in this region';
      if (state.busy) message = 'loading…';
      else if (!state.classes.size) message = 'every read class is switched off';
      else if (region && Object.values(region.counts).some((n) => n > 0)) {
        message = 'reads here, but none match the current filter';
      }
      ctx.font = '13px ui-sans-serif, system-ui, sans-serif';
      ctx.fillStyle = palette.ink_muted;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText(message, width / 2, height / 2);
    }
    drawStrandLabel(ctx, strand, width, height, palette);
    restore();
    return;
  }

  const lineWidth = clamp(layout.laneH * 0.62, 1, 6);
  // Fatter than a plain read so an inverted piece stands out, but never wider
  // than its own lane: 1.7x the bar exceeds MAX_LANE and bled over the reads
  // above and below it.  render/mapfig.py caps the same way, against
  // --line-space.
  const inversionWidth = Math.min(Math.max(lineWidth * 1.7, 1.5),
                                  layout.laneH * 0.95);
  const arrowSize = clamp(layout.laneH * ARROW_LANE_FRACTION,
                          ARROW_MIN_PX, ARROW_MAX_PX);
  const markerStroke = clamp(arrowSize * MARKER_STROKE_FRACTION,
                             MARKER_STROKE_MIN, MARKER_STROKE_MAX);
  const paths = new Map();     // colour -> Path2D
  const chevronPaths = new Map();
  const tickPaths = new Map();
  const rings = [];
  const gapPath = new Path2D();
  const inversionPath = new Path2D();
  const joinMarks = [];

  const pathFor = (colour) => {
    let path = paths.get(colour);
    if (!path) { path = new Path2D(); paths.set(colour, path); }
    return path;
  };

  for (const read of reads) {
    const y = yOf(read.lane);
    const colour = palette.classes[read.cls];

    let left = Infinity;
    let right = -Infinity;
    let previousEnd = null;

    for (const segment of read.segments) {
      const x0 = xOf(segment[0]);
      const x1 = xOf(segment[1]);
      if (x1 < -20 || x0 > width + 20) { previousEnd = x1; continue; }
      const target = segment[2] ? inversionPath : pathFor(colour);
      target.moveTo(x0, y);
      target.lineTo(Math.max(x1, x0 + 0.7), y);
      if (previousEnd !== null && x0 - previousEnd > 1) {
        gapPath.moveTo(previousEnd, y);
        gapPath.lineTo(x0, y);
      }
      previousEnd = x1;
      left = Math.min(left, x0);
      right = Math.max(right, x1);
    }

    if (left === Infinity) continue;
    if (read.joins.length) joinMarks.push([right, y, read.joins.join(', ')]);

    boxes.push({
      id: read.id, x0: xOf(read.start), x1: xOf(read.end),
      y: y, read: read,
    });
  }

  // The skipped reference inside a divided read, drawn faintly so the pieces
  // read as one molecule without competing with the alignments themselves.
  ctx.strokeStyle = palette.grid;
  ctx.lineWidth = Math.max(0.5, lineWidth * 0.4);
  ctx.stroke(gapPath);

  ctx.lineWidth = lineWidth;
  ctx.lineCap = 'butt';
  for (const [colour, path] of paths) {
    ctx.strokeStyle = colour;
    ctx.stroke(path);
  }

  ctx.strokeStyle = palette.inversion;
  ctx.lineWidth = inversionWidth;
  ctx.stroke(inversionPath);

  // Arrowheads, junction ticks and circular-join rings, on top of the lines.
  const inversionArrows = new Path2D();
  drawMarkers(reads, yOf, layout, arrowSize, chevronPaths, tickPaths, rings, palette,
              inversionArrows);

  ctx.lineWidth = markerStroke;
  for (const [colour, path] of tickPaths) {
    ctx.strokeStyle = colour;
    ctx.stroke(path);
  }

  ctx.lineJoin = 'round';
  ctx.lineCap = 'round';
  for (const [colour, path] of chevronPaths) {
    ctx.strokeStyle = colour;
    ctx.stroke(path);
  }

  if (rings.length) {
    const radius = clamp(layout.laneH * 0.42, 2, 4.5);
    ctx.lineWidth = markerStroke;
    for (const [x, y, colour] of rings) {
      ctx.strokeStyle = colour;
      ctx.fillStyle = palette.surface;
      ctx.beginPath();
      ctx.arc(x, y, radius, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    }
  }

  // Inversion arrows last, over everything: a pale halo first so a red arrow
  // stays legible where it crosses the lanes of other reads.
  ctx.strokeStyle = palette.surface;
  ctx.lineWidth = 4;
  ctx.stroke(inversionArrows);
  ctx.strokeStyle = palette.inversion;
  ctx.lineWidth = 2;
  ctx.stroke(inversionArrows);
  ctx.lineCap = 'butt';

  if (joinMarks.length <= 60 && layout.laneH >= 5) {
    ctx.font = '9px ui-monospace, Consolas, monospace';
    ctx.fillStyle = palette.events.contig_join;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    for (const [x, y, text] of joinMarks) ctx.fillText(' →' + text, x, y);
  }

  drawSelection(ctx, boxes, palette);
  drawStrandLabel(ctx, strand, width, height, palette);
  restore();
}

// ---------------------------------------------------------------------------
// Arrowheads
// ---------------------------------------------------------------------------

// One open chevron, tip at (x, y), opening away from `dir`.  Half-height is
// half the head length, the same proportion theme.py's _chevron uses.
function chevron(path, x, y, dir, size) {
  const h = size * 0.5;
  path.moveTo(x - dir * size, y - h);
  path.lineTo(x, y);
  path.lineTo(x - dir * size, y + h);
}

// The printed map's inversion indicator: a shafted arrow at -30 degrees, coming
// in from above-left with its tip just short of the inverted piece's left end
// (render/theme.py INVERSION_ARROW_*, placed by layout.py).  It points at the
// piece rather than sitting on it, so even a piece a pixel or two wide stays
// visible.  Sized in pixels, not lanes, so it stays readable when lanes are thin.
const INVERSION_ANGLE = Math.PI / 6;
function inversionArrow(path, x, y, laneH) {
  const length = clamp(laneH * 2.4, 14, 24);
  const head = clamp(length * 0.4, 6, 8);
  const ux = Math.cos(INVERSION_ANGLE);
  const uy = Math.sin(INVERSION_ANGLE);   // canvas y grows downwards
  const tipX = x - 1.5;
  const tipY = y - Math.max(1.5, laneH * 0.4);
  path.moveTo(tipX - ux * length, tipY - uy * length);
  path.lineTo(tipX, tipY);
  // The head: two barbs, each 28 degrees off the shaft.
  for (const side of [-1, 1]) {
    const a = INVERSION_ANGLE + Math.PI + side * 0.49;
    path.moveTo(tipX, tipY);
    path.lineTo(tipX + Math.cos(a) * head, tipY + Math.sin(a) * head);
  }
}

function drawMarkers(reads, yOf, layout, size, chevronPaths, tickPaths, rings,
                     palette, inversionArrows) {
  const width = plotWidth();
  const showArrows = layout.laneH >= ARROW_MIN_LANE;

  const pathFor = (map, colour) => {
    let path = map.get(colour);
    if (!path) { path = new Path2D(); map.set(colour, path); }
    return path;
  };

  const emit = (style, dir, x, y, colour) => {
    if (style === 'inversion_arrow') {
      inversionArrow(inversionArrows, x, y, layout.laneH);
      return;
    }
    const path = pathFor(chevronPaths, colour);
    chevron(path, x, y, dir, size);
    if (style === 'arrow_ct') {
      // The double chevron: this read carries on onto another contig.
      chevron(path, x - dir * size * 0.85, y, dir, size);
    } else if (style === 'arrow_circular') {
      rings.push([x + dir * size * 1.5, y, colour]);
    }
  };

  for (const read of reads) {
    const y = yOf(read.lane);
    const colour = palette.classes[read.cls];

    for (const [position, style, dir, atEdge] of (read.markers || [])) {
      // atEdge markers sit at the edge of the *fetched* window, which is three
      // times wider than the view. They are for the image exporter, which asks
      // for exactly the range it draws; here they are replaced below.
      if (atEdge) continue;

      const x = xOf(position);
      if (x < -40 || x > width + 40) continue;

      if (style === 'tick') {
        const h = Math.max(1.5, layout.laneH * 0.45);
        const path = pathFor(tickPaths, colour);
        path.moveTo(x, y - h);
        path.lineTo(x, y + h);
        continue;
      }
      // An inversion is the one mark worth seeing at every zoom, and it sits
      // beside the line rather than on it, so it does not need a tall lane.
      if (showArrows || style === 'inversion_arrow') emit(style, dir, x, y, colour);
    }

    // "This read carries on past the edge of what you are looking at" -- the
    // browser's equivalent of the printed map's page-fold arrows. Only the
    // client knows where its viewport ends, so only the client can place them.
    // Reads that lie wholly outside the view are in the payload too (the
    // fetched window is wider); they must not leave a chevron on the edge.
    if (!showArrows) continue;
    if (read.end < state.start || read.start > state.end) continue;
    if (read.start < state.start) {
      emit(read.cont, read.cont_dir || -1, 0, y, colour);
    }
    if (read.end > state.end) {
      emit(read.cont, read.cont_dir || 1, width, y, colour);
    }
  }
}

// Every line the selected read occupies, not just the first: both ends of a
// structural event light up together, which is the whole point of selecting it.
function drawSelection(ctx, boxes, palette) {
  if (state.selected === null) return;
  const found = boxes.filter((b) => b.id === state.selected);
  if (!found.length) return;

  ctx.strokeStyle = palette.ink;
  ctx.lineWidth = 1;
  ctx.setLineDash([3, 2]);
  const pad = Math.max(3, state.laneH * 0.7);
  for (const box of found) {
    ctx.strokeRect(box.x0 - 2, box.y - pad, (box.x1 - box.x0) + 4, pad * 2);
  }
  ctx.setLineDash([]);
}

// Against the axis in both panes -- bottom of the forward one, top of the
// reverse one -- so the label stays next to the reference line it qualifies
// however far the pane is scrolled.  On its own patch of surface: a dense
// pileup runs right to the edge and would otherwise swallow it.
function drawStrandLabel(ctx, strand, width, height, palette) {
  const text = strand > 0 ? '+ strand  (forward)' : '− strand  (reverse)';
  ctx.font = '10px ui-monospace, Consolas, monospace';
  const boxWidth = ctx.measureText(text).width + 10;
  const top = strand > 0 ? height - 15 : 2;

  ctx.fillStyle = palette.surface;
  ctx.fillRect(width - boxWidth - 3, top, boxWidth, 13);
  ctx.fillStyle = palette.ink_muted;
  ctx.textAlign = 'right';
  ctx.textBaseline = 'top';
  ctx.fillText(text, width - 6, top + 1);
}

function redraw() {
  if (!state.meta || !state.contig) return;
  drawHeader();
  drawPileup();
  updateStatus();
}

// ---------------------------------------------------------------------------
// Data
// ---------------------------------------------------------------------------

function filterKey() {
  const classes = [...state.classes].sort().join(',');
  return state.contig.name + '|' + classes + '|' + (state.invertedOnly ? '1' : '0');
}

function regionIsUsable() {
  const key = state.regionKey;
  if (!key || key.filters !== filterKey()) return false;
  if (state.start < key.start || state.end > key.end) return false;
  const visible = state.end - state.start;
  return visible * REFETCH_ZOOM >= (key.end - key.start);
}

function scheduleRefresh() {
  syncHash();
  // Draw straight away from whatever is cached -- panning inside the fetched
  // window is then instant, and the debounced fetch only fills in the edges.
  redraw();
  clearTimeout(refreshTimer);
  refreshTimer = setTimeout(fetchAll, DEBOUNCE_MS);
}

async function fetchAll() {
  if (!state.contig) return;
  await Promise.all([fetchRegion(), fetchDepth()]);
  redraw();
}

async function fetchRegion(force) {
  if (!state.contig) return;
  if (!force && regionIsUsable()) return;

  const visible = state.end - state.start;
  const pad = visible * (FETCH_FACTOR - 1) / 2;
  const start = Math.max(0, Math.floor(state.start - pad));
  const end = Math.min(state.contig.length, Math.ceil(state.end + pad));

  if (regionAbort) regionAbort.abort();
  regionAbort = new AbortController();
  setBusy(1);
  try {
    const query = new URLSearchParams({
      contig: state.contig.name,
      start: String(start),
      end: String(end),
      // What is on screen, as opposed to the three screens being fetched: the
      // lane gap is a pixel quantity and has to be measured against the view.
      span: String(Math.max(1, Math.round(visible))),
      classes: [...state.classes].join(','),
    });
    if (state.invertedOnly) query.set('inverted', '1');
    state.region = await getJSON('/api/region?' + query, regionAbort.signal);
    state.regionKey = { start, end, filters: filterKey() };
  } catch (error) {
    if (error.name !== 'AbortError') showError(error);
  } finally {
    setBusy(-1);
  }
}

async function fetchDepth() {
  if (!state.contig) return;
  if (depthAbort) depthAbort.abort();
  depthAbort = new AbortController();
  setBusy(1);
  try {
    const query = new URLSearchParams({
      contig: state.contig.name,
      start: String(Math.floor(state.start)),
      end: String(Math.ceil(state.end)),
      bins: String(Math.max(100, Math.round(plotWidth()))),
    });
    state.depth = await getJSON('/api/depth?' + query, depthAbort.signal);
  } catch (error) {
    if (error.name !== 'AbortError') showError(error);
  } finally {
    setBusy(-1);
  }
}

function showError(error) {
  const notes = el('notes');
  const text = escapeHTML(String(error.message || error));
  if (notes) notes.innerHTML = '<p><strong>Error:</strong> ' + text + '</p>';
  else setStatusMessage('error: ' + (error.message || error));
}

// ---------------------------------------------------------------------------
// Navigation
// ---------------------------------------------------------------------------

function setView(start, end, options) {
  const length = state.contig.length || 1;
  let span = clamp(end - start, MIN_SPAN, length);
  let lo = start;
  if (span >= length) { lo = 0; span = length; }
  lo = clamp(lo, 0, length - span);
  state.start = Math.round(lo);
  state.end = Math.round(lo + span);
  if (!options || !options.quiet) scheduleRefresh();
}

function zoomAt(factor, anchorX) {
  const span = state.end - state.start;
  const anchor = anchorX === undefined
    ? (state.start + state.end) / 2
    : positionOf(anchorX);
  const fraction = (anchor - state.start) / span;
  const next = clamp(span * factor, MIN_SPAN, state.contig.length);
  setView(anchor - fraction * next, anchor - fraction * next + next);
}

function panBy(fraction) {
  const span = state.end - state.start;
  setView(state.start + span * fraction, state.end + span * fraction);
}

function selectContig(name, keepView) {
  const contig = state.meta.contigs.find((c) => c.name === name);
  if (!contig) return;
  state.contig = contig;
  el('contig').value = name;
  state.region = null;
  state.regionKey = null;
  state.depth = null;
  state.selected = null;
  if (!keepView) setView(0, contig.length, { quiet: true });
  scheduleRefresh();
}

function goTo(result) {
  if (result.contig !== state.contig.name) selectContig(result.contig, true);
  setView(result.start, result.end);
  if (result.read !== undefined && result.read !== null) selectRead(result.read);
}

async function jumpEvent(direction) {
  const centre = Math.round((state.start + state.end) / 2);
  const query = new URLSearchParams({
    contig: state.contig.name,
    from: String(centre),
    dir: direction > 0 ? 'next' : 'prev',
  });
  try {
    const payload = await getJSON('/api/events?' + query);
    if (!payload.event) {
      setStatusMessage('no structural events in this run');
      return;
    }
    const event = payload.event;
    if (event.contig !== state.contig.name) selectContig(event.contig, true);
    const span = Math.max(MIN_SPAN, Math.min(state.end - state.start, 40000));
    setView(event.position - span / 2, event.position + span / 2);
    setStatusMessage(
      (state.meta.palette.event_labels[event.kind] || event.kind) +
      ' at ' + commas(event.position) +
      (event.read_label ? '  (' + event.read_label + ')' : ''));
    if (event.read !== null && event.read !== undefined) selectRead(event.read);
  } catch (error) {
    showError(error);
  }
}

// ---------------------------------------------------------------------------
// Status
// ---------------------------------------------------------------------------

let statusMessageTimer = null;

function setStatusMessage(text) {
  const node = el('status-counts');
  node.textContent = text;
  clearTimeout(statusMessageTimer);
  statusMessageTimer = setTimeout(updateStatus, 4000);
}

function updateStatus() {
  el('status-region').textContent =
    state.contig.name + ':' + commas(state.start) + '–' + commas(state.end) +
    '   (' + shortBp(state.end - state.start) + ')';

  const region = state.region;
  if (!region) { el('status-counts').textContent = ''; return; }

  // Count what is actually on screen.  region.counts covers the fetched window,
  // which is three times wider than the view, so quoting it here would put a
  // number in front of the user that does not match what they can see.
  //
  // By id: region.reads carries one entry per drawn line, and a long-divided
  // read is several of them.  Counting entries would report one read twice for
  // showing both ends of the same event.
  const visible = {};
  const counted = new Set();
  for (const read of region.reads) {
    if (read.end < state.start || read.start > state.end) continue;
    if (counted.has(read.id)) continue;
    counted.add(read.id);
    visible[read.cls] = (visible[read.cls] || 0) + 1;
  }

  const parts = [];
  for (const cls of state.meta.classes) {
    const count = visible[cls] || 0;
    if (count) parts.push(shortClassName(cls) + ': ' + commas(count));
  }
  if (!parts.length) parts.push('no reads in view');
  if (region.truncated) {
    const shown = Object.values(region.shown).reduce((a, b) => a + b, 0);
    const total = Object.values(region.counts).reduce((a, b) => a + b, 0);
    parts.push('density limit: ' + commas(shown) + ' of ' + commas(total) +
               ' loaded — zoom in to see them all');
  }
  el('status-counts').textContent = parts.join('    ');
}

// ---------------------------------------------------------------------------
// Detail panel
// ---------------------------------------------------------------------------

function escapeHTML(text) {
  return String(text).replace(/[&<>"']/g, (c) => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
  })[c]);
}

async function selectRead(id) {
  state.selected = id;
  state.selectedAnno = null;
  drawPileup();
  drawHeader();
  try {
    const detail = await getJSON('/api/read?id=' + encodeURIComponent(id));
    renderDetail(detail);
  } catch (error) {
    showError(error);
  }
}

async function selectAnnotation(id) {
  state.selectedAnno = id;
  state.selected = null;
  drawHeader();
  drawPileup();
  try {
    const detail = await getJSON('/api/annotation?id=' + encodeURIComponent(id));
    renderAnnotationDetail(detail);
  } catch (error) {
    showError(error);
  }
}

function renderAnnotationDetail(anno) {
  const palette = state.meta.palette;
  const colour = palette.annotations[anno.group] !== undefined
    ? palette.annotations[anno.group]
    : palette.annotations[''];
  const symbol = palette.annotation_symbols[anno.group] || '';

  const badges = ['<span class="badge" style="background:' + colour + '">' +
                  escapeHTML(anno.feature || 'feature') + symbol + '</span>'];
  const groupLabel = (palette.annotation_legend.find((e) => e.group === anno.group) || {}).label;
  if (anno.group && groupLabel) {
    badges.push('<span class="badge" style="background:' + colour + '">' +
                escapeHTML(groupLabel) + '</span>');
  }

  // Only the rows GenBank actually filled: an empty /note is not information,
  // and a column of em-dashes reads as though the parser lost something.
  const rows = [
    ['name', anno.name],
    ['locus tag', anno.locus_tag && anno.locus_tag !== anno.name ? anno.locus_tag : ''],
    ['product', anno.product],
    ['note', anno.note],
    ['contig', anno.contig],
    ['position', commas(anno.start) + '–' + commas(anno.end)],
    ['length', shortBp(anno.length)],
    ['strand', anno.strand > 0 ? '+' : '−'],
  ].filter(([, value]) => value !== '' && value !== null && value !== undefined);

  el('panel-body').innerHTML =
    '<h2>' + escapeHTML(anno.label || anno.name || anno.feature) + '</h2>' +
    '<p>' + badges.join(' ') + '</p>' +
    '<h3>GenBank feature</h3>' +
    '<dl>' + rows.map(([k, v]) =>
      '<dt>' + escapeHTML(k) + '</dt><dd>' + escapeHTML(v) + '</dd>').join('') +
    '</dl>' +
    '<h3>Actions</h3><div class="buttons">' +
    '<button id="zoom-gene">zoom to this feature</button>' +
    '<button id="clear-gene">clear</button></div>' +
    '<div id="notes"></div>';

  el('zoom-gene').onclick = () => {
    const pad = Math.max(500, (anno.end - anno.start) * 0.25);
    if (anno.contig !== state.contig.name) selectContig(anno.contig, true);
    setView(anno.start - pad, anno.end + pad);
  };
  el('clear-gene').onclick = clearSelection;
}

function renderDetail(read) {
  const palette = state.meta.palette;
  const colour = palette.classes[read.cls];
  const rows = read.hits.map((hit) => {
    const identity = hit.identity === null ? '—'
      : (hit.identity * 100).toFixed(2) + '%';
    return '<tr>' +
      '<td>' + escapeHTML(hit.contig) + '</td>' +
      '<td class="num">' + commas(hit.ref_start) + '</td>' +
      '<td class="num">' + commas(hit.ref_end) + '</td>' +
      '<td class="num">' + (hit.strand > 0 ? '+' : '−') + '</td>' +
      '<td class="num">' + commas(hit.qstart) + '–' + commas(hit.qend) + '</td>' +
      '<td class="num">' + identity + '</td>' +
      '</tr>';
  }).join('');

  const badges = ['<span class="badge" style="background:' + colour + '">' +
                  escapeHTML(palette.class_labels[read.cls]) + '</span>'];
  if (read.inverted) {
    badges.push('<span class="badge" style="background:' + palette.inversion +
                '">inversion</span>');
  }
  if (read.junction !== 'C') {
    badges.push('<span class="badge" style="background:' +
                (read.junction === 'Circle' ? palette.events.circle
                                            : palette.events.contig_join) +
                '">' + escapeHTML(read.junction) + '</span>');
  }

  const identityNote = read.has_identity ? ''
    : '<p class="muted">Alignment identity needs <code>--identity</code> on the ' +
      'run that built this cache.</p>';

  el('panel-body').innerHTML =
    '<h2>' + escapeHTML(read.label) + '</h2>' +
    '<p>' + badges.join(' ') + '</p>' +
    '<h3>Read</h3>' +
    '<dl>' +
    '<dt>name</dt><dd>' + escapeHTML(read.name) + '</dd>' +
    '<dt>length</dt><dd>' + commas(read.qlen) + ' bp</dd>' +
    '<dt>matched</dt><dd>' + commas(read.matched_bases) + ' bp</dd>' +
    '<dt>alignments</dt><dd>' + read.hits.length + '</dd>' +
    (read.cls === 'undivided' ? ''
      : '<dt>division</dt><dd>' + commas(read.distance) + ' bp</dd>') +
    '</dl>' +
    '<h3>Alignments</h3>' +
    '<table class="hits"><thead><tr>' +
    '<th>contig</th><th>start</th><th>end</th><th>str</th><th>in read</th><th>ident</th>' +
    '</tr></thead><tbody>' + rows + '</tbody></table>' +
    identityNote +
    '<h3>Actions</h3><div class="buttons">' +
    '<button id="zoom-read">zoom to this read</button>' +
    '<button id="clear-read">clear</button></div>' +
    '<div id="notes"></div>';

  el('zoom-read').onclick = () => {
    const here = read.hits.filter((h) => h.contig === state.contig.name);
    const lo = Math.min(...here.map((h) => h.ref_start));
    const hi = Math.max(...here.map((h) => h.ref_end));
    if (Number.isFinite(lo) && Number.isFinite(hi)) {
      const pad = Math.max(500, (hi - lo) * 0.1);
      setView(lo - pad, hi + pad);
    }
  };
  el('clear-read').onclick = clearSelection;
}

function clearSelection() {
  state.selected = null;
  state.selectedAnno = null;
  renderIdlePanel();
  drawPileup();
  drawHeader();
}

// A small inline drawing of one arrow style, for the legend.  Same forms the
// canvas draws, so the key and the map agree.
function arrowGlyph(style, colour) {
  const stroke = ' fill="none" stroke="' + colour +
                 '" stroke-width="1.6" stroke-linejoin="round" stroke-linecap="round"/>';
  const chev = (cx, dir) => '<polyline points="' +
    (cx - dir * 5) + ',3.5 ' + cx + ',8 ' + (cx - dir * 5) + ',12.5"' + stroke;
  const shaft = (x0, x1) => '<line x1="' + x0 + '" y1="8" x2="' + x1 + '" y2="8"' + stroke;

  let body;
  if (style === 'arrow_ct') body = shaft(2, 17) + chev(17, 1) + chev(12, 1);
  else if (style === 'arrow_circular') {
    body = shaft(2, 14) + chev(14, 1) +
      '<circle cx="21" cy="8" r="3"' + stroke;
  } else if (style === 'arrow_rev') body = shaft(9, 24) + chev(9, -1);
  else if (style === 'inversion_arrow') {
    // The diagonal arrow pointing at the left end of a short red piece.
    body = '<line x1="15" y1="13" x2="24" y2="13" stroke="' + colour +
      '" stroke-width="3"/>' +
      '<polyline points="2,3 13,9.5"' + stroke +
      '<polyline points="8.6,9.9 13,9.5 11.1,5.4"' + stroke;
  }
  else if (style === 'tick') body = shaft(2, 24) + '<line x1="13" y1="2.5" x2="13" y2="13.5"' + stroke;
  else body = shaft(2, 17) + chev(17, 1);

  return '<svg class="glyph" viewBox="0 0 26 16" width="26" height="16" ' +
         'aria-hidden="true">' + body + '</svg>';
}

// The same triangle drawEventMarks puts under the coordinate axis, so the key
// and the header are visibly the same mark rather than a dot standing in for
// one.
function eventGlyph(colour) {
  return '<svg class="glyph mark" viewBox="0 0 16 16" width="16" height="16" ' +
         'aria-hidden="true"><polygon points="8,2.5 2.5,13 13.5,13" fill="' +
         colour + '"/></svg>';
}

function renderIdlePanel() {
  const meta = state.meta;
  const summary = meta.summary || {};
  const counts = summary.counts || {};
  const palette = meta.palette;

  const runRows = [
    ['sample', meta.sample],
    ['reads', commas(summary.total_reads || 0)],
    ['coverage', (summary.mean_coverage || 0).toFixed(2) + '×'],
    ['read N50', commas(summary.n50 || 0) + ' bp'],
    ['events', commas(summary.events || 0)],
    ['built', meta.created],
  ];
  if (meta.project && meta.project.organism) {
    runRows.unshift(['organism', meta.project.organism]);
  }

  const legendRows = meta.classes.map((cls) =>
    '<div class="row"><span class="swatch" style="background:' +
    palette.classes[cls] + '"></span>' + escapeHTML(palette.class_labels[cls]) +
    '  <span class="muted">' + commas(counts[cls] || 0) + '</span></div>').join('') +
    '<div class="row"><span class="swatch" style="background:' + palette.inversion +
    '"></span>Inverted segment  <span class="muted">' +
    commas(summary.inverted_reads || 0) + ' reads</span></div>' +
    (meta.is_genbank
      ? palette.annotation_legend.map((entry) =>
        '<div class="row"><span class="swatch" style="background:' +
        palette.annotations[entry.group] + '"></span>' +
        escapeHTML(entry.label) + '</div>').join('')
      : '');

  // The triangles under the coordinate axis, spelled out: what each colour is
  // evidence of, not just what it is called.
  const descriptions = palette.event_descriptions || {};
  const eventRows = meta.event_kinds.map((kind) =>
    '<div class="row event">' + eventGlyph(palette.events[kind]) +
    '<span><strong>' + escapeHTML(palette.event_labels[kind]) + '</strong>' +
    (descriptions[kind]
      ? '<br><span class="muted">' + escapeHTML(descriptions[kind]) + '</span>'
      : '') +
    '</span></div>').join('');

  const ink = palette.classes.undivided;
  const arrowRows = [
    ['arrow_fwd', ink, 'Read runs on past here, left to right'],
    ['arrow_rev', ink, 'Read runs on past here, right to left'],
    ['arrow_ct', palette.events.contig_join, 'Continues onto another contig'],
    ['arrow_circular', palette.events.circle, 'Bridges the origin — circular'],
    ['inversion_arrow', palette.inversion, 'This piece runs the other way'],
    ['tick', palette.classes.short_divided, 'Junction inside a divided read'],
  ].map(([style, colour, text]) =>
    '<div class="row">' + arrowGlyph(style, colour) +
    escapeHTML(text) + '</div>').join('');

  el('panel-body').innerHTML =
    '<h2>Nothing selected</h2>' +
    '<p class="muted">Click a read to see its name, its class and every place ' +
    'it aligned. Click a gene under the axis for its GenBank entry — product, ' +
    'note and locus tag. Hover either for a summary.</p>' +
    (meta.is_genbank
      ? '<p class="muted">Features with no gene name are drawn by the number ' +
        'from their locus tag: <code>SS37A_42600</code> appears as ' +
        '<code>42600</code>. The full tag is in the panel.</p>'
      : '') +
    '<h3>Run</h3><dl>' +
    runRows.map(([k, v]) => '<dt>' + escapeHTML(k) + '</dt><dd>' +
                escapeHTML(v || '—') + '</dd>').join('') +
    '</dl>' +
    '<h3>Legend</h3>' + legendRows +
    '<h3>Marks under the axis</h3>' +
    '<p class="muted">One triangle per structural event, at the coordinate ' +
    'it was found. Hover one on the ruler for the same note; click it to open ' +
    'the read; <code>n</code> and <code>p</code> step through them.</p>' +
    eventRows +
    '<h3>Arrows</h3>' + arrowRows +
    '<h3>Keys</h3><dl class="keys">' +
    '<dt>drag</dt><dd>pan</dd>' +
    '<dt>wheel</dt><dd>zoom at the cursor</dd>' +
    '<dt>shift+wheel</dt><dd>scroll the lanes under the cursor</dd>' +
    '<dt>← →</dt><dd>pan</dd>' +
    '<dt>+ −</dt><dd>zoom</dd>' +
    '<dt>f</dt><dd>whole contig</dd>' +
    '<dt>n p</dt><dd>next / previous event</dd>' +
    '<dt>esc</dt><dd>clear the selection</dd>' +
    '</dl>' +
    '<div id="notes">' +
    (meta.notes || []).map((note) => '<p>' + escapeHTML(note) + '</p>').join('') +
    '</div>';
}

// ---------------------------------------------------------------------------
// Pointer interaction
// ---------------------------------------------------------------------------

let drag = null;

function canvasX(event, canvas) {
  return event.clientX - canvas.getBoundingClientRect().left;
}

function hitTest(x, y, strand) {
  const tolerance = Math.max(3, state.laneH * 0.6);
  let best = null;
  for (const box of (state.hitboxes[strandKey(strand)] || [])) {
    if (Math.abs(box.y - y) > tolerance) continue;
    if (x < box.x0 - 3 || x > box.x1 + 3) continue;
    if (!best || (box.x1 - box.x0) < (best.x1 - best.x0)) best = box;
  }
  return best;
}

// `strand` is +1 / -1 for a read pane and 0 for the header, which pans and
// zooms like the rest but has no reads to hit.
function installPointer(canvas, strand) {
  const scroll = strand === 0 ? null : scrollFor(strand);
  canvas.addEventListener('pointerdown', (event) => {
    if (event.button !== 0) return;
    canvas.setPointerCapture(event.pointerId);
    drag = {
      x: canvasX(event, canvas),
      start: state.start,
      end: state.end,
      moved: false,
      pointerId: event.pointerId,
    };
    el('view').classList.add('panning');
  });

  canvas.addEventListener('pointermove', (event) => {
    const x = canvasX(event, canvas);

    if (drag) {
      const span = drag.end - drag.start;
      const shift = (drag.x - x) / plotWidth() * span;
      if (Math.abs(drag.x - x) > 2) drag.moved = true;
      setView(drag.start + shift, drag.end + shift);
      return;
    }

    el('status-cursor').textContent = commas(positionOf(x)) + ' bp';
    const rect = canvas.getBoundingClientRect();

    // The header: an event mark on the axis, or a feature in the annotation
    // rows.  They sit in different bands, so there is nothing to disambiguate.
    if (!scroll) {
      const localY = event.clientY - rect.top;
      const mark = eventAt(x, localY);
      const gene = mark ? null : annotationAt(x, localY);
      el('view').classList.toggle('over-read', Boolean(mark || gene));
      if (mark) showEventTooltip(event, mark);
      else if (gene) showAnnotationTooltip(event, gene.anno);
      else hideTooltip();
      return;
    }

    const box = hitTest(x, event.clientY - rect.top, strand);
    el('view').classList.toggle('over-read', Boolean(box));
    if (box) showTooltip(event, box.read); else hideTooltip();
  });

  const finish = (event) => {
    if (!drag) return;
    const wasClick = !drag.moved;
    const x = canvasX(event, canvas);
    try { canvas.releasePointerCapture(drag.pointerId); } catch (_) { /* ignore */ }
    drag = null;
    el('view').classList.remove('panning');

    if (!wasClick) return;
    const rect = canvas.getBoundingClientRect();

    if (!scroll) {
      const localY = event.clientY - rect.top;
      // Clicking a mark selects the read that produced it, which is the
      // follow-up question every one of them raises.
      const mark = eventAt(x, localY);
      if (mark && mark.read !== null && mark.read !== undefined) {
        selectRead(mark.read);
        return;
      }
      const gene = mark ? null : annotationAt(x, localY);
      if (gene) selectAnnotation(gene.id);
      return;
    }

    const box = hitTest(x, event.clientY - rect.top, strand);
    if (box) selectRead(box.id); else clearSelection();
  };

  canvas.addEventListener('pointerup', finish);
  canvas.addEventListener('pointercancel', finish);
  canvas.addEventListener('pointerleave', () => {
    if (!drag) { hideTooltip(); el('view').classList.remove('over-read'); }
  });

  canvas.addEventListener('wheel', (event) => {
    // deltaMode 1 counts lines, not pixels; without this a mouse wheel barely
    // moves while a trackpad flies.
    const delta = event.deltaMode === 1 ? event.deltaY * 16
                : event.deltaMode === 2 ? event.deltaY * 400
                : event.deltaY;

    // Zoom owns the wheel, because that is what this view is for. When a pane
    // holds more lanes than fit, shift+wheel scrolls the one under the cursor
    // instead -- otherwise preventDefault would leave its far lanes
    // unreachable.
    if (event.shiftKey && scroll) {
      scroll.scrollTop += delta;
      event.preventDefault();
      return;
    }

    event.preventDefault();
    zoomAt(Math.pow(1.0016, clamp(delta, -600, 600)), canvasX(event, canvas));
  }, { passive: false });

  canvas.addEventListener('dblclick', (event) => {
    zoomAt(0.5, canvasX(event, canvas));
  });
}

// A GenBank feature, from what the region payload already carries -- no
// request, so it appears the instant the pointer crosses the bar.  The rest
// (full product, note, locus tag) arrives from /api/annotation on a click.
function showAnnotationTooltip(event, anno) {
  const palette = state.meta.palette;
  const symbol = palette.annotation_symbols[anno.group] || '';
  const lines = [anno.name + symbol];

  if (anno.label && anno.label !== anno.name) lines.push('drawn as ' + anno.label);
  lines.push(anno.feature + (anno.strand > 0 ? '   + strand' : '   − strand'));
  lines.push(commas(anno.start) + '–' + commas(anno.end) +
             '  (' + shortBp(anno.end - anno.start) + ')');
  if (anno.product) lines.push('', ...wrapText(anno.product, 46));
  lines.push('', 'click for the full entry');

  placeTooltip(event, lines.join('\n'));
}

// What one of the triangles under the axis is saying.  The wording comes from
// theme.EVENT_DESCRIPTIONS via /api/meta, so the tooltip, the panel key and
// any future printed legend cannot explain the same mark three ways.
function showEventTooltip(event, mark) {
  const palette = state.meta.palette;
  const description = (palette.event_descriptions || {})[mark.kind] || '';
  const lines = [
    palette.event_labels[mark.kind] || mark.kind,
    commas(mark.position) + ' bp',
  ];
  if (mark.read_label) lines.push('read ' + mark.read_label + '   (click to open)');
  if (description) lines.push('', ...wrapText(description, 46));
  placeTooltip(event, lines.join('\n'));
}

// The tooltip is `white-space: pre`, so a long sentence would run off the side
// of the window instead of wrapping.  Break it here, on whole words.
function wrapText(text, columns) {
  const lines = [];
  let line = '';
  for (const word of text.split(/\s+/)) {
    if (line && (line + ' ' + word).length > columns) { lines.push(line); line = word; }
    else line = line ? line + ' ' + word : word;
  }
  if (line) lines.push(line);
  return lines;
}

function showTooltip(event, read) {
  const palette = state.meta.palette;
  const lines = [
    read.label,
    palette.class_labels[read.cls],
  ];
  if (read.inverted) lines.push('inversion');
  if (read.junction === 'Circle') lines.push('circular join — bridges the origin');
  else if (read.junction === 'Contig_Join') lines.push('contig join');
  if (read.joins.length) lines.push('also on ' + read.joins.join(', '));
  lines.push(commas(read.start) + '–' + commas(read.end) +
             '  (' + shortBp(read.end - read.start) + ')');
  // read.hits is the whole read's count on this contig; a long-divided read is
  // drawn as one line per piece, so read.segments is only this line's share.
  const total = read.hits === undefined ? read.segments.length : read.hits;
  lines.push(total + ' alignment' + (total === 1 ? '' : 's') +
             (total > read.segments.length ? ' on this contig' : '') +
             '   ' + (read.strand > 0 ? '+ strand' : '− strand'));

  placeTooltip(event, lines.join('\n'));
}

function placeTooltip(event, text) {
  tooltip.textContent = text;
  tooltip.hidden = false;
  const host = el('view').getBoundingClientRect();
  const width = tooltip.offsetWidth;
  const height = tooltip.offsetHeight;
  let x = event.clientX - host.left + 14;
  let y = event.clientY - host.top + 14;
  if (x + width > host.width - 8) x = event.clientX - host.left - width - 14;
  if (y + height > host.height - 8) y = event.clientY - host.top - height - 14;
  tooltip.style.left = Math.max(4, x) + 'px';
  tooltip.style.top = Math.max(4, y) + 'px';
}

function hideTooltip() {
  tooltip.hidden = true;
}

// ---------------------------------------------------------------------------
// Search
// ---------------------------------------------------------------------------

function renderSuggestions() {
  const box = el('suggestions');
  if (!suggestions.length) { box.hidden = true; box.innerHTML = ''; return; }
  box.innerHTML = suggestions.map((result, index) =>
    '<div data-index="' + index + '"' +
    (index === suggestionIndex ? ' class="active"' : '') + '>' +
    '<span class="kind">' + escapeHTML(result.kind) + '</span>' +
    '<span>' + escapeHTML(result.label) + '</span></div>').join('');
  box.hidden = false;
  for (const node of box.children) {
    node.addEventListener('mousedown', (event) => {
      event.preventDefault();
      pickSuggestion(Number(node.dataset.index));
    });
  }
}

function pickSuggestion(index) {
  const result = suggestions[index];
  if (!result) return;
  goTo(result);
  suggestions = [];
  suggestionIndex = -1;
  renderSuggestions();
  el('search').blur();
}

async function runSearch(text) {
  if (!text.trim()) { suggestions = []; renderSuggestions(); return; }
  try {
    const payload = await getJSON('/api/search?q=' + encodeURIComponent(text));
    suggestions = payload.results;
    suggestionIndex = suggestions.length ? 0 : -1;
    renderSuggestions();
  } catch (error) {
    showError(error);
  }
}

// ---------------------------------------------------------------------------
// URL hash
// ---------------------------------------------------------------------------

let hashGuard = false;
let hashTimer = null;

// Throttled: a drag fires pointermove at frame rate, and browsers rate-limit
// (and warn about) replaceState called that often.
function syncHash() {
  if (hashTimer) return;
  hashTimer = setTimeout(() => {
    hashTimer = null;
    const hash = '#' + state.contig.name + ':' + state.start + '-' + state.end;
    if (location.hash === hash) return;
    hashGuard = true;
    history.replaceState(null, '', location.pathname + hash);
    setTimeout(() => { hashGuard = false; }, 0);
  }, 250);
}

function applyHash() {
  const raw = decodeURIComponent(location.hash.replace(/^#/, ''));
  if (!raw) return false;
  const cut = raw.lastIndexOf(':');
  if (cut < 0) return false;
  const name = raw.slice(0, cut);
  const bounds = raw.slice(cut + 1).split('-');
  const contig = state.meta.contigs.find((c) => c.name === name);
  if (!contig || bounds.length !== 2) return false;
  const start = Number(bounds[0]);
  const end = Number(bounds[1]);
  if (!Number.isFinite(start) || !Number.isFinite(end) || end <= start) return false;
  state.contig = contig;
  el('contig').value = contig.name;
  setView(start, end, { quiet: true });
  return true;
}

// ---------------------------------------------------------------------------
// Export
// ---------------------------------------------------------------------------

// Fetched rather than navigated to: a plain location.assign would replace the
// page with a JSON error document if the render ever failed, throwing away the
// view the user was trying to export.
async function exportView(format) {
  const query = new URLSearchParams({
    contig: state.contig.name,
    start: String(state.start),
    end: String(state.end),
    format: format,
    classes: [...state.classes].join(','),
  });
  if (state.invertedOnly) query.set('inverted', '1');

  setBusy(1);
  setStatusMessage('rendering ' + format.toUpperCase() + '…');
  try {
    const response = await fetch('/api/export?' + query);
    if (!response.ok) {
      const payload = await response.json();
      throw new Error(payload.error || response.statusText);
    }
    const url = URL.createObjectURL(await response.blob());
    const link = document.createElement('a');
    link.href = url;
    link.download = (state.meta.sample || 'readrift') + '_' + state.contig.name +
                    '_' + state.start + '-' + state.end + '.' + format;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 30000);
    setStatusMessage('saved ' + link.download);
  } catch (error) {
    showError(error);
  } finally {
    setBusy(-1);
  }
}

// ---------------------------------------------------------------------------
// Wiring
// ---------------------------------------------------------------------------

function buildFilters() {
  const palette = state.meta.palette;
  const host = el('filters');
  host.innerHTML = state.meta.classes.map((cls) =>
    '<label title="' + escapeHTML(palette.class_labels[cls]) + '">' +
    '<input type="checkbox" data-class="' + cls + '" checked>' +
    '<span class="swatch" style="background:' + palette.classes[cls] + '"></span>' +
    escapeHTML(shortClassName(cls)) + '</label>').join('') +
    '<label title="Only reads with a segment running the other way">' +
    '<input type="checkbox" id="only-inverted">' +
    '<span class="swatch" style="background:' + palette.inversion + '"></span>' +
    'inverted only</label>';

  for (const input of host.querySelectorAll('input[data-class]')) {
    input.addEventListener('change', () => {
      if (input.checked) state.classes.add(input.dataset.class);
      else state.classes.delete(input.dataset.class);
      state.regionKey = null;
      scheduleRefresh();
    });
  }
  el('only-inverted').addEventListener('change', (event) => {
    state.invertedOnly = event.target.checked;
    state.regionKey = null;
    scheduleRefresh();
  });
}

function shortClassName(cls) {
  return { undivided: 'undivided', short_divided: 'short', long_divided: 'long' }[cls] || cls;
}

function wire() {
  el('contig').addEventListener('change', (event) => selectContig(event.target.value));
  el('zoom-in').addEventListener('click', () => zoomAt(0.5));
  el('zoom-out').addEventListener('click', () => zoomAt(2));
  el('fit').addEventListener('click', () => setView(0, state.contig.length));
  el('next-event').addEventListener('click', () => jumpEvent(1));
  el('prev-event').addEventListener('click', () => jumpEvent(-1));
  el('export-pdf').addEventListener('click', () => exportView('pdf'));
  el('export-png').addEventListener('click', () => exportView('png'));

  const search = el('search');
  search.addEventListener('input', () => {
    clearTimeout(searchTimer);
    searchTimer = setTimeout(() => runSearch(search.value), 180);
  });
  search.addEventListener('keydown', (event) => {
    if (event.key === 'ArrowDown') {
      event.preventDefault();
      suggestionIndex = Math.min(suggestionIndex + 1, suggestions.length - 1);
      renderSuggestions();
    } else if (event.key === 'ArrowUp') {
      event.preventDefault();
      suggestionIndex = Math.max(suggestionIndex - 1, 0);
      renderSuggestions();
    } else if (event.key === 'Enter') {
      event.preventDefault();
      if (suggestionIndex >= 0) {
        pickSuggestion(suggestionIndex);
      } else {
        // Enter before the debounce fired: search, then go, so one keystroke
        // never leaves the user staring at a list they have to click.
        runSearch(search.value).then(() => {
          if (suggestions.length) pickSuggestion(0);
        });
      }
    } else if (event.key === 'Escape') {
      suggestions = [];
      renderSuggestions();
      search.blur();
    }
  });
  search.addEventListener('blur', () => {
    setTimeout(() => { suggestions = []; renderSuggestions(); }, 120);
  });

  document.addEventListener('keydown', (event) => {
    if (event.target.tagName === 'INPUT' || event.target.tagName === 'SELECT') return;
    if (event.ctrlKey || event.metaKey || event.altKey) return;
    switch (event.key) {
      case 'ArrowLeft': event.preventDefault(); panBy(-0.25); break;
      case 'ArrowRight': event.preventDefault(); panBy(0.25); break;
      case '+': case '=': event.preventDefault(); zoomAt(0.5); break;
      case '-': case '_': event.preventDefault(); zoomAt(2); break;
      case 'f': case 'F': setView(0, state.contig.length); break;
      case 'n': case 'N': jumpEvent(1); break;
      case 'p': case 'P': jumpEvent(-1); break;
      case 'Escape': clearSelection(); break;
      case '/': event.preventDefault(); el('search').focus(); break;
      default: break;
    }
  });

  window.addEventListener('hashchange', () => {
    if (hashGuard) return;
    if (applyHash()) scheduleRefresh();
  });

  let resizeTimer = null;
  window.addEventListener('resize', () => {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(() => { redraw(); fetchDepth().then(redraw); }, 120);
  });

  installPointer(headerCanvas, 0);
  installPointer(plusCanvas, 1);
  installPointer(minusCanvas, -1);
}

async function boot() {
  try {
    state.meta = await getJSON('/api/meta');
  } catch (error) {
    document.body.innerHTML =
      '<p style="padding:24px;font-family:monospace">Cannot reach the readrift ' +
      'server: ' + escapeHTML(String(error.message || error)) + '</p>';
    return;
  }

  const meta = state.meta;
  el('sample').textContent = meta.sample || meta.cache;
  document.title = 'readrift — ' + (meta.sample || 'browser');

  el('contig').innerHTML = meta.contigs.map((c) =>
    '<option value="' + escapeHTML(c.name) + '">' + escapeHTML(c.name) +
    '  (' + shortBp(c.length) + ')</option>').join('');

  for (const cls of meta.classes) state.classes.add(cls);
  buildFilters();
  wire();

  if (!applyHash()) {
    state.contig = meta.contigs[0];
    if (!state.contig) {
      document.body.innerHTML =
        '<p style="padding:24px;font-family:monospace">This cache has no ' +
        'contigs.</p>';
      return;
    }
    el('contig').value = state.contig.name;
    setView(0, state.contig.length, { quiet: true });
  }

  renderIdlePanel();
  await fetchAll();
  redraw();
}

boot();
