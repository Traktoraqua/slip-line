"""Ruled-table detection and reconstruction (Phase 6).

Vector rule segments from :mod:`pdf2tex.extract` are snapped into a grid: a
table region is a cluster of intersecting horizontal/vertical rules forming at
least a 2 × 2 cell grid.  Spans are assigned to cells by containment; missing
interior rules merge unit cells into ``\\multicolumn`` / ``\\multirow`` spans.
Column alignment defaults to left, switching to right for numeric-majority
columns.  ``Table`` captions (``^Table\\s+\\d``) adjacent to the region are
attached.
"""

from __future__ import annotations

import logging
import re
import statistics

from .models import Cell, Page, Rule, Table, TableGrid

log = logging.getLogger(__name__)

SNAP_TOL = 3.0           # pt; rule coordinates within this snap together
COVER_RATIO = 0.6        # interior rule must cover this fraction of a cell edge
CAPTION_ADJACENCY = 2.5  # caption within this many line heights of the region

_CAPTION_RE = re.compile(r"^(Table|Tabell)\s+[\dA-Z][-.\d]*", re.IGNORECASE)
_NUMERIC_RE = re.compile(r"^[-+(]?\d[\d\s.,%)/×x–-]*$")


# --------------------------------------------------------------------------- #
# Region clustering (union-find over rules)
# --------------------------------------------------------------------------- #
def _bbox(rule: Rule) -> tuple[float, float, float, float]:
    return (min(rule.x0, rule.x1), min(rule.y0, rule.y1),
            max(rule.x0, rule.x1), max(rule.y0, rule.y1))


def _touch(a: Rule, b: Rule) -> bool:
    ax0, ay0, ax1, ay1 = _bbox(a)
    bx0, by0, bx1, by1 = _bbox(b)
    t = SNAP_TOL
    return not (ax1 + t < bx0 or bx1 + t < ax0 or ay1 + t < by0 or by1 + t < ay0)


def _cluster_rules(rules: list[Rule]) -> list[list[Rule]]:
    parent = list(range(len(rules)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        parent[find(i)] = find(j)

    for i in range(len(rules)):
        for j in range(i + 1, len(rules)):
            if _touch(rules[i], rules[j]):
                union(i, j)

    groups: dict[int, list[Rule]] = {}
    for i, rule in enumerate(rules):
        groups.setdefault(find(i), []).append(rule)
    return list(groups.values())


# --------------------------------------------------------------------------- #
# Grid lines
# --------------------------------------------------------------------------- #
def _snap(values: list[float]) -> list[float]:
    reps: list[float] = []
    for v in sorted(values):
        if not reps or v - reps[-1] > SNAP_TOL:
            reps.append(v)
        else:
            reps[-1] = (reps[-1] + v) / 2.0  # average within tolerance
    return reps


def _cover(intervals: list[tuple[float, float]], lo: float, hi: float) -> bool:
    """True if *intervals* cover ≥ COVER_RATIO of ``[lo, hi]``."""
    span = hi - lo
    if span <= 0:
        return False
    covered = 0.0
    for a, b in sorted(intervals):
        a, b = max(a, lo), min(b, hi)
        if b > a:
            covered += b - a
    return covered >= COVER_RATIO * span


# --------------------------------------------------------------------------- #
# Cell text
# --------------------------------------------------------------------------- #
def _cell_lines(page: Page, x0: float, y0: float, x1: float, y1: float) -> list[str]:
    spans = []
    for sp in page.spans:
        cx = (sp.bbox[0] + sp.bbox[2]) / 2.0
        cy = (sp.bbox[1] + sp.bbox[3]) / 2.0
        if x0 <= cx <= x1 and y0 <= cy <= y1:
            spans.append(sp)
    if not spans:
        return []
    spans.sort(key=lambda s: (round(s.origin[1]), s.bbox[0]))
    lines: list[list] = []
    for sp in spans:
        if lines and abs(sp.origin[1] - lines[-1][0]) <= SNAP_TOL:
            lines[-1][1].append(sp)
        else:
            lines.append([sp.origin[1], [sp]])
    out = []
    for _, group in lines:
        text = " ".join(s.text.strip() for s in group if s.text.strip())
        text = re.sub(r"\s+", " ", text).strip()
        if text:
            out.append(text)
    return out


# --------------------------------------------------------------------------- #
# Grid construction
# --------------------------------------------------------------------------- #
def _build_grid(page: Page, hrules: list[Rule], vrules: list[Rule]) -> TableGrid | None:
    ys = _snap([(_bbox(r)[1] + _bbox(r)[3]) / 2.0 for r in hrules])
    xs = _snap([(_bbox(r)[0] + _bbox(r)[2]) / 2.0 for r in vrules])
    nrows, ncols = len(ys) - 1, len(xs) - 1
    if nrows < 2 or ncols < 2:
        return None

    # Vertical-rule x → list of covered y-intervals, and similarly for h-rules.
    v_at: dict[int, list[tuple[float, float]]] = {c: [] for c in range(ncols + 1)}
    for r in vrules:
        bx0, by0, bx1, by1 = _bbox(r)
        xc = (bx0 + bx1) / 2.0
        for c, gx in enumerate(xs):
            if abs(xc - gx) <= SNAP_TOL:
                v_at[c].append((by0, by1))
                break
    h_at: dict[int, list[tuple[float, float]]] = {r: [] for r in range(nrows + 1)}
    for r in hrules:
        bx0, by0, bx1, by1 = _bbox(r)
        yc = (by0 + by1) / 2.0
        for ri, gy in enumerate(ys):
            if abs(yc - gy) <= SNAP_TOL:
                h_at[ri].append((bx0, bx1))
                break

    def vline_present(c: int, r: int) -> bool:
        if c == 0 or c == ncols:
            return True
        return _cover(v_at[c], ys[r], ys[r + 1])

    def hline_present(r: int, c: int) -> bool:
        if r == 0 or r == nrows:
            return True
        return _cover(h_at[r], xs[c], xs[c + 1])

    # Union-find over unit cells, merging across missing interior rules.
    n = nrows * ncols
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        parent[find(i)] = find(j)

    def idx(r: int, c: int) -> int:
        return r * ncols + c

    for r in range(nrows):
        for c in range(ncols):
            if c + 1 < ncols and not vline_present(c + 1, r):
                union(idx(r, c), idx(r, c + 1))
            if r + 1 < nrows and not hline_present(r + 1, c):
                union(idx(r, c), idx(r + 1, c))

    comps: dict[int, list[tuple[int, int]]] = {}
    for r in range(nrows):
        for c in range(ncols):
            comps.setdefault(find(idx(r, c)), []).append((r, c))

    cells: list[Cell] = []
    for members in comps.values():
        rs = [m[0] for m in members]
        cs = [m[1] for m in members]
        r0, r1, c0, c1 = min(rs), max(rs), min(cs), max(cs)
        lines = _cell_lines(page, xs[c0], ys[r0], xs[c1 + 1], ys[r1 + 1])
        cells.append(Cell(r0=r0, c0=c0, r1=r1, c1=c1, lines=lines))

    col_aligns = _column_alignments(cells, ncols, nrows)
    for cell in cells:
        cell.align = "c" if cell.colspan > 1 else col_aligns[cell.c0]

    return TableGrid(nrows=nrows, ncols=ncols, cells=cells, col_aligns=col_aligns)


def _column_alignments(cells: list[Cell], ncols: int, nrows: int) -> list[str]:
    aligns = ["l"] * ncols
    for c in range(ncols):
        texts = [
            " ".join(cell.lines)
            for cell in cells
            if cell.c0 == c and cell.colspan == 1 and cell.r0 > 0 and cell.lines
        ]
        if texts and sum(bool(_NUMERIC_RE.match(t)) for t in texts) * 2 >= len(texts):
            aligns[c] = "r"
    return aligns


# --------------------------------------------------------------------------- #
# Captions
# --------------------------------------------------------------------------- #
def _median_line_height(page: Page) -> float:
    hs = [ln.bbox[3] - ln.bbox[1] for ln in page.lines if ln.bbox[3] > ln.bbox[1]]
    return statistics.median(hs) if hs else 12.0


def _find_caption(page: Page, region: tuple[float, float, float, float]):
    rx0, ry0, rx1, ry1 = region
    margin = CAPTION_ADJACENCY * _median_line_height(page)
    for line in page.lines:
        text = line.text.strip()
        m = _CAPTION_RE.match(text)
        if not m:
            continue
        ly0, ly1 = line.bbox[1], line.bbox[3]
        if ry0 - margin <= ly1 <= ry0 + margin or ry1 - margin <= ly0 <= ry1 + margin:
            # Strip the "Table N" label so LaTeX's \caption adds its own.
            body = text[m.end():].lstrip(" .:)-–—\t") or text
            return line, body
    return None, None


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def detect_tables(page: Page) -> list[dict]:
    """Return detected tables as dicts: table, bbox, top_y, consumed lines."""
    if not page.rules:
        return []
    results: list[dict] = []
    for cluster in _cluster_rules(page.rules):
        hrules = [r for r in cluster if r.horizontal]
        vrules = [r for r in cluster if r.vertical]
        if len(hrules) < 2 or len(vrules) < 2:
            continue
        grid = _build_grid(page, hrules, vrules)
        if grid is None:
            continue
        boxes = [_bbox(r) for r in cluster]
        region = (min(b[0] for b in boxes), min(b[1] for b in boxes),
                  max(b[2] for b in boxes), max(b[3] for b in boxes))
        cap_line, cap_text = _find_caption(page, region)

        consumed = set()
        for line in page.lines:
            cx = (line.bbox[0] + line.bbox[2]) / 2.0
            cy = (line.bbox[1] + line.bbox[3]) / 2.0
            if region[0] <= cx <= region[2] and region[1] <= cy <= region[3]:
                consumed.add(id(line))
        if cap_line is not None:
            consumed.add(id(cap_line))

        results.append({
            "table": Table(grid=grid, caption=cap_text),
            "bbox": region,
            "top_y": region[1],
            "consumed": consumed,
        })
    log.debug("Detected %d table(s) on page %d", len(results), page.number)
    return results
