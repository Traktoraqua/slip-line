"""Document structure assembly (Phase 3).

Turns extracted, zone-filtered :class:`~pdf2tex.models.Page` objects into a
logical :class:`~pdf2tex.models.Document`:

* **Body font** — the character-weighted dominant font size; the yardstick for
  "larger than body" comparisons.
* **Headings** — numbered lines (``^\\d+(\\.\\d+)* …``) that are larger and/or
  bold map to ``\\section`` … ``\\paragraph`` by numbering depth; unnumbered
  bold, standalone, larger lines become starred sections.
* **Paragraphs** — lines are grouped within a block, split on a vertical gap
  > 1.5 × median leading or an indent increase; lines are joined with
  dehyphenation (``word-`` + lowercase continuation), and spans become
  bold/italic/plain inline runs.
"""

from __future__ import annotations

import logging
import re
import statistics
from collections import Counter

from .lists import build_list_block, match_marker, take_list
from .tables import detect_tables
from .models import (
    Bold,
    DocMeta,
    Document,
    Heading,
    Inline,
    Italic,
    Line,
    Page,
    Paragraph,
    Text,
)

log = logging.getLogger(__name__)

# Heading numbering: "1", "1.2", "2.3.1" followed by whitespace and content.
# Note: "1. Introduction" (dot before space) is a numbered *list* item, not a
# heading, and intentionally does not match.
_NUM_RE = re.compile(r"^(\d+(?:\.\d+)*)\s+\S")

INDENT_TOL = 10.0          # pt; larger indent than the group start ⇒ new para
GAP_FACTOR = 1.5           # vertical gap > this × median leading ⇒ new para
HEADING_SIZE_RATIO = 1.05  # "larger than body" threshold for numbered headings
UNNUMBERED_RATIO = 1.15    # size threshold for unnumbered headings
BIG_HEADING_RATIO = 1.30   # promote unnumbered heading to level 1
MAX_HEADING_CHARS = 60     # unnumbered heading must be short


# --------------------------------------------------------------------------- #
# Document-level metrics
# --------------------------------------------------------------------------- #
def _body_size(pages: list[Page]) -> float:
    counter: Counter[int] = Counter()
    for page in pages:
        for span in page.spans:
            n = len(span.text.strip())
            if n:
                counter[round(span.size)] += n
    return float(counter.most_common(1)[0][0]) if counter else 11.0


def _median_leading(pages: list[Page]) -> float:
    gaps: list[float] = []
    for page in pages:
        for block in page.blocks:
            tops = [ln.bbox[1] for ln in block.lines if ln.spans]
            gaps += [b - a for a, b in zip(tops, tops[1:]) if b - a > 0]
    return statistics.median(gaps) if gaps else 14.0


# --------------------------------------------------------------------------- #
# Inline run assembly
# --------------------------------------------------------------------------- #
_STYLE_CTOR = {"text": Text, "bold": Bold, "italic": Italic}


def _classify_span(span) -> str:
    if span.is_bold and not span.is_italic:
        return "bold"
    if span.is_italic:
        return "italic"
    return "text"


def _first_alpha_lower(text: str) -> bool:
    stripped = text.lstrip()
    return bool(stripped) and stripped[0].islower()


def _build_inlines(lines: list[Line]) -> list[Inline]:
    """Join a paragraph group's spans into styled inline runs.

    Consecutive same-style spans merge; line breaks add a space unless the
    previous line ends with a hyphen and the next begins lowercase, in which
    case the hyphen is dropped (dehyphenation).
    """
    segs: list[list[str]] = []  # [style, text]
    for li, line in enumerate(lines):
        spans = [s for s in line.spans if s.text]
        if not spans:
            continue
        if segs:
            prev = segs[-1]
            if prev[1].rstrip().endswith("-") and _first_alpha_lower(spans[0].text):
                prev[1] = prev[1].rstrip()[:-1]
            elif not prev[1].endswith(" "):
                prev[1] = prev[1] + " "
        for span in spans:
            style = _classify_span(span)
            if segs and segs[-1][0] == style:
                segs[-1][1] += span.text
            else:
                segs.append([style, span.text])

    if segs:
        segs[0][1] = segs[0][1].lstrip()
        segs[-1][1] = segs[-1][1].rstrip()
    return [_STYLE_CTOR[st](text=tx) for st, tx in segs if tx]


# --------------------------------------------------------------------------- #
# Headings
# --------------------------------------------------------------------------- #
def _bold_majority(line: Line) -> bool:
    total = sum(len(s.text.strip()) for s in line.spans)
    bold = sum(len(s.text.strip()) for s in line.spans if s.is_bold)
    return total > 0 and bold * 2 >= total


def _line_size(line: Line, default: float) -> float:
    return max((s.size for s in line.spans if s.text.strip()), default=default)


def _heading_from_group(group: list[Line], body_size: float) -> Heading | None:
    if len(group) != 1:
        return None
    line = group[0]
    text = line.text.strip()
    if not text:
        return None
    size = _line_size(line, body_size)
    bold = _bold_majority(line)

    m = _NUM_RE.match(text)
    if m and (size >= body_size * HEADING_SIZE_RATIO or bold):
        number = m.group(1)
        rest = text[m.end(1):].strip()
        if rest:
            return Heading(level=number.count(".") + 1, number=number, text=rest)

    if (
        bold
        and text[0].isalpha()  # exclude numbered/list-marker lines
        and size >= body_size * UNNUMBERED_RATIO
        and len(text) <= MAX_HEADING_CHARS
        and text[-1] not in ".:;,"
    ):
        level = 1 if size >= body_size * BIG_HEADING_RATIO else 2
        return Heading(level=level, number=None, text=text, starred=True)
    return None


# --------------------------------------------------------------------------- #
# Title
# --------------------------------------------------------------------------- #
def _largest_span(pages: list[Page]):
    best = None
    for span in (pages[0].spans if pages else []):
        if span.text.strip() and (best is None or span.size > best.size):
            best = span
    return best


def detect_title(pages: list[Page]) -> str | None:
    """Heuristic title: the line carrying the largest text on page 1."""
    span = _largest_span(pages)
    if span is None:
        return None
    sizes = [s.size for p in pages for s in p.spans if s.text.strip()]
    if not sizes:
        return None
    body = sorted(sizes)[len(sizes) // 2]
    if span.size <= body * 1.15:
        return None
    for line in pages[0].lines:
        if span in line.spans:
            return line.text.strip()
    return span.text.strip()


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def _rows_for_page(page: Page) -> list[tuple[Line, bool]]:
    """Flatten a page to ``(line, new_block)`` rows in reading order."""
    rows: list[tuple[Line, bool]] = []
    for block in page.blocks:
        for li, line in enumerate(block.lines):
            if line.spans:
                rows.append((line, li == 0))
    return rows


def _segment_groups(
    rows: list[tuple[Line, bool]], leading: float
) -> list[list[Line]]:
    """Split a non-list row segment into paragraph groups.

    A new group starts on a block boundary, a vertical gap > 1.5 × leading, or
    an indent increase — mirroring per-block grouping over a flat sequence.
    """
    groups: list[list[Line]] = []
    cur: list[Line] = []
    cur_x = 0.0
    for line, new_block in rows:
        x0, y0 = line.bbox[0], line.bbox[1]
        if cur and (
            new_block
            or (y0 - cur[-1].bbox[1]) > GAP_FACTOR * leading
            or x0 > cur_x + INDENT_TOL
        ):
            groups.append(cur)
            cur = []
        if not cur:
            cur_x = x0
        cur.append(line)
    if cur:
        groups.append(cur)
    return groups


def build_document(
    pages: list[Page],
    meta: DocMeta,
    detect_title_flag: bool = True,
) -> Document:
    title = detect_title(pages) if detect_title_flag else None
    if title:
        meta.title = title

    body_size = _body_size(pages)
    leading = _median_leading(pages)
    log.debug("body size=%.1f, median leading=%.1f", body_size, leading)

    elements: list = []
    headings = 0
    lists = 0

    def emit_group(group: list[Line], page_number: int) -> None:
        nonlocal headings
        gtext = " ".join(ln.text.strip() for ln in group).strip()
        if not gtext:
            return
        if title and page_number == 1 and gtext == title:
            return  # consumed as the title
        heading = _heading_from_group(group, body_size)
        if heading is not None:
            if title and page_number == 1 and heading.text == title:
                return
            elements.append(heading)
            headings += 1
            return
        inlines = _build_inlines(group)
        if inlines:
            elements.append(Paragraph(inlines=inlines))

    def process_rows(rows: list[tuple[Line, bool]], page_number: int) -> None:
        nonlocal lists
        i = 0
        while i < len(rows):
            if match_marker(rows[i][0].text):
                end = take_list(rows, i, leading)
                elements.append(build_list_block(rows, i, end, _build_inlines))
                lists += 1
                i = end
            else:
                j = i
                while j < len(rows) and not match_marker(rows[j][0].text):
                    j += 1
                for group in _segment_groups(rows[i:j], leading):
                    emit_group(group, page_number)
                i = j

    tables = 0
    for page in pages:
        tdicts = detect_tables(page)
        consumed: set[int] = set()
        for td in tdicts:
            consumed |= td["consumed"]

        items: list[tuple] = []  # (kind, payload, new_block, y)
        for block in page.blocks:
            for li, line in enumerate(block.lines):
                if line.spans and id(line) not in consumed:
                    items.append(("line", line, li == 0, line.bbox[1]))
        for td in tdicts:
            items.append(("table", td["table"], False, td["top_y"]))
        items.sort(key=lambda it: it[3])

        buf: list[tuple[Line, bool]] = []
        for kind, payload, new_block, _ in items:
            if kind == "line":
                buf.append((payload, new_block))
            else:
                process_rows(buf, page.number)
                buf = []
                elements.append(payload)
                tables += 1
        process_rows(buf, page.number)

    log.info(
        "Built document: %d element(s) (%d heading(s), %d list(s), %d table(s))",
        len(elements),
        headings,
        lists,
        tables,
    )
    return Document(meta=meta, elements=elements)
