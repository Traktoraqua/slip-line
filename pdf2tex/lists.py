"""List detection: bullets, numbering and indent-based nesting (Phase 5).

A *run* of consecutive lines is a list when the first line carries a marker
(``• …``, ``1. …``, ``a) …``) and following lines are either further markers
or indented continuation lines.  Marker start-x positions are clustered into
nesting levels (≤ 4, the LaTeX limit) to build nested
``itemize``/``enumerate`` blocks.

Inline assembly is injected as a callback so this module depends only on the
data model (avoids a cycle with :mod:`pdf2tex.structure`).
"""

from __future__ import annotations

import dataclasses
import re
from typing import Callable

from .models import Inline, Line, ListBlock

# Bullet glyphs that start an itemize item.
_BULLET_RE = re.compile(r"^([•◦▪‣·–—\-\*])\s+")
# Numbered/lettered/roman markers: "1.", "(a)", "iv)".
_NUMBER_RE = re.compile(r"^\(?(?:\d+|[a-zA-Z]|[ivxlcdmIVXLCDM]+)[.)]\s+")

LEVEL_TOL = 3.0          # pt; marker start-x within this is the same level
CONT_INDENT = 2.0        # pt; continuation lines indent past the marker column
CONT_GAP_FACTOR = 2.0    # continuation must be within this × leading vertically
MAX_DEPTH = 3            # 0-based; 4 nesting levels

InlineBuilder = Callable[[list[Line]], list[Inline]]


def match_marker(text: str) -> tuple[str, int] | None:
    """Return ``(kind, marker_length)`` for a list line, else ``None``.

    ``kind`` is ``"itemize"`` or ``"enumerate"``; ``marker_length`` covers the
    marker and its trailing whitespace.
    """
    stripped = text.lstrip()
    if not stripped:
        return None
    lead = len(text) - len(stripped)
    m = _BULLET_RE.match(stripped)
    if m:
        return "itemize", lead + m.end()
    m = _NUMBER_RE.match(stripped)
    if m:
        return "enumerate", lead + m.end()
    return None


def take_list(rows: list, start: int, leading: float) -> int:
    """Return the exclusive end index of the list run beginning at *start*.

    *rows* is a list of ``(line, new_block)`` tuples.
    """
    base_x = rows[start][0].bbox[0]
    last_y = rows[start][0].bbox[1]
    i = start + 1 if match_marker(rows[start][0].text) else start
    while i < len(rows):
        line = rows[i][0]
        x0, y0 = line.bbox[0], line.bbox[1]
        if match_marker(line.text):
            base_x = min(base_x, x0)
            last_y = y0
            i += 1
            continue
        # Continuation: indented past the marker column and close vertically.
        if x0 > base_x + CONT_INDENT and (y0 - last_y) <= CONT_GAP_FACTOR * leading:
            last_y = y0
            i += 1
            continue
        break
    return i


def _strip_prefix(line: Line, n: int) -> Line:
    """Drop the first *n* characters (the marker) across the line's spans."""
    spans = []
    remaining = n
    for sp in line.spans:
        if remaining <= 0:
            spans.append(sp)
        elif len(sp.text) <= remaining:
            remaining -= len(sp.text)
        else:
            spans.append(dataclasses.replace(sp, text=sp.text[remaining:]))
            remaining = 0
    if not spans:
        spans = [dataclasses.replace(line.spans[0], text="")]
    return Line(spans=spans, bbox=line.bbox)


def _clusters(xs: list[float]) -> list[float]:
    reps: list[float] = []
    for x in sorted(set(round(v, 1) for v in xs)):
        if not reps or x - reps[-1] > LEVEL_TOL:
            reps.append(x)
    return reps


def build_list_block(
    rows: list,
    start: int,
    end: int,
    build_inlines: InlineBuilder,
) -> ListBlock:
    """Build a (possibly nested) ListBlock from ``rows[start:end]``."""
    marker_xs = [
        rows[k][0].bbox[0]
        for k in range(start, end)
        if match_marker(rows[k][0].text)
    ]
    clusters = _clusters(marker_xs) or [rows[start][0].bbox[0]]

    def level_of(x: float) -> int:
        lvl = 0
        for idx, c in enumerate(clusters):
            if x >= c - LEVEL_TOL:
                lvl = idx
        return min(lvl, MAX_DEPTH)

    stack: list[ListBlock] = []
    cur_lb: ListBlock | None = None
    cur_idx = -1
    cur_lines: list[Line] = []

    for k in range(start, end):
        line = rows[k][0]
        marker = match_marker(line.text)
        if marker:
            kind, mlen = marker
            level = min(level_of(line.bbox[0]), len(stack))  # grow ≤ 1 level
            if level == len(stack):
                nb = ListBlock(kind=kind, items=[], level=level)
                if stack:
                    stack[-1].items.append(nb)
                stack.append(nb)
            else:
                del stack[level + 1:]
            lb = stack[level]
            item_line = _strip_prefix(line, mlen)
            cur_lines = [item_line]
            lb.items.append(build_inlines(cur_lines))
            cur_lb, cur_idx = lb, len(lb.items) - 1
        elif cur_lb is not None:  # continuation
            cur_lines.append(line)
            cur_lb.items[cur_idx] = build_inlines(cur_lines)

    return stack[0]
