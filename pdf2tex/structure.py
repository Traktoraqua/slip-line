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
from .mathconv.detect import (
    build_display_math,
    has_math_indicator,
    is_equation_number,
    is_math_span,
    line_ref,
    shift_class,
    spans_to_latex,
)
from .tables import detect_tables
from .models import (
    Bold,
    DocMeta,
    Document,
    Heading,
    Inline,
    InlineMath,
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
    """Join a paragraph group's spans into styled / math inline runs.

    Consecutive same-style spans merge; line breaks add a space unless the
    previous line ends with a hyphen and the next begins lowercase (which
    dehyphenates).  Math-flagged spans form ``$…$`` runs; a leading
    sub/superscript steals its base token from the preceding text run.
    """
    segs: list[dict] = []  # {"kind": text|bold|italic|math, "text"|"items"}
    for li, line in enumerate(lines):
        spans = [s for s in line.spans if s.text]
        if not spans:
            continue
        ref_y, ref_size = line_ref(spans)
        if segs and segs[-1]["kind"] in _STYLE_CTOR:
            prev = segs[-1]
            if prev["text"].rstrip().endswith("-") and _first_alpha_lower(spans[0].text):
                prev["text"] = prev["text"].rstrip()[:-1]
            elif not prev["text"].endswith(" "):
                prev["text"] += " "
        for span in spans:
            cls = shift_class(span, ref_y, ref_size)
            if is_math_span(span) or cls != "normal":
                if segs and segs[-1]["kind"] == "math":
                    segs[-1]["items"].append((span, cls))
                else:
                    segs.append({"kind": "math", "items": [(span, cls)]})
            else:
                style = _classify_span(span)
                if segs and segs[-1]["kind"] == style:
                    segs[-1]["text"] += span.text
                else:
                    segs.append({"kind": style, "text": span.text})

    _steal_math_bases(segs)

    inlines: list[Inline] = []
    last = len(segs) - 1
    for idx, seg in enumerate(segs):
        if seg["kind"] == "math":
            latex = seg.get("base", "") + spans_to_latex(seg["items"])
            latex = latex.strip()
            if latex:
                inlines.append(InlineMath(latex=latex))
        else:
            text = seg["text"]
            if idx == 0:
                text = text.lstrip()
            if idx == last:
                text = text.rstrip()
            if text:
                inlines.append(_STYLE_CTOR[seg["kind"]](text=text))
    return inlines


def _steal_math_bases(segs: list[dict]) -> None:
    """Move the base token of a sub/superscript into its math run.

    A run like ``x`` (text) followed by a subscript ``i`` should render as
    ``$x_{i}$``, not ``x$_{i}$`` (invalid).
    """
    import re as _re

    for idx, seg in enumerate(segs):
        if seg["kind"] != "math" or idx == 0:
            continue
        if seg["items"][0][1] == "normal":
            continue  # run starts with a real math token, no base needed
        prev = segs[idx - 1]
        if prev["kind"] not in _STYLE_CTOR or not prev["text"]:
            continue
        m = _re.search(r"(\S+)$", prev["text"])
        if not m:
            continue
        from .mathconv.unicode_map import convert_math
        seg["base"] = convert_math(m.group(1))
        prev["text"] = prev["text"][: m.start()]


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


def _merge_equation_numbers(kept: list[list], page: Page) -> None:
    """Fold a lone right-margin equation number into its sibling math line.

    Equation numbers often sit in their own block on the same visual line as
    the equation; merging lets display detection emit a numbered ``equation``.
    """
    remove: set[int] = set()
    for idx, (line, _) in enumerate(kept):
        if idx in remove or len(line.spans) != 1:
            continue
        if not is_equation_number(line.text) or line.bbox[0] <= page.width * 0.7:
            continue
        cy = (line.bbox[1] + line.bbox[3]) / 2.0
        best, best_dx = None, None
        for jdx, (other, _) in enumerate(kept):
            if jdx == idx or jdx in remove:
                continue
            ocy = (other.bbox[1] + other.bbox[3]) / 2.0
            if abs(ocy - cy) <= 6.0 and other.bbox[2] < line.bbox[0]:
                if has_math_indicator(other.spans):
                    dx = line.bbox[0] - other.bbox[2]
                    if best is None or dx < best_dx:
                        best, best_dx = jdx, dx
        if best is not None:
            other = kept[best][0]
            kept[best][0] = Line(
                spans=other.spans + line.spans,
                bbox=(other.bbox[0], min(other.bbox[1], line.bbox[1]),
                      line.bbox[2], max(other.bbox[3], line.bbox[3])),
            )
            remove.add(idx)
    for idx in sorted(remove, reverse=True):
        del kept[idx]


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
    maths = 0
    eq_index = 0

    def emit_group(group: list[Line], page: Page) -> None:
        nonlocal headings, maths, eq_index
        gtext = " ".join(ln.text.strip() for ln in group).strip()
        if not gtext:
            return
        if title and page.number == 1 and gtext == title:
            return  # consumed as the title

        dm = build_display_math(group, page, page.number, eq_index)
        if dm is not None:
            if dm.number is not None:
                eq_index += 1
            elements.append(dm)
            maths += 1
            return

        heading = _heading_from_group(group, body_size)
        if heading is not None:
            if title and page.number == 1 and heading.text == title:
                return
            elements.append(heading)
            headings += 1
            return
        inlines = _build_inlines(group)
        if inlines:
            elements.append(Paragraph(inlines=inlines))

    def process_rows(rows: list[tuple[Line, bool]], page: Page) -> None:
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
                    emit_group(group, page)
                i = j

    tables = 0
    for page in pages:
        tdicts = detect_tables(page)
        consumed: set[int] = set()
        for td in tdicts:
            consumed |= td["consumed"]

        kept: list[list] = []  # [line, new_block]
        for block in page.blocks:
            for li, line in enumerate(block.lines):
                if line.spans and id(line) not in consumed:
                    kept.append([line, li == 0])
        _merge_equation_numbers(kept, page)

        items: list[tuple] = []  # (kind, payload, new_block, y)
        for line, new_block in kept:
            items.append(("line", line, new_block, line.bbox[1]))
        for td in tdicts:
            items.append(("table", td["table"], False, td["top_y"]))
        items.sort(key=lambda it: it[3])

        buf: list[tuple[Line, bool]] = []
        for kind, payload, new_block, _ in items:
            if kind == "line":
                buf.append((payload, new_block))
            else:
                process_rows(buf, page)
                buf = []
                elements.append(payload)
                tables += 1
        process_rows(buf, page)

    log.info(
        "Built document: %d element(s) "
        "(%d heading(s), %d list(s), %d table(s), %d display math)",
        len(elements),
        headings,
        lists,
        tables,
        maths,
    )
    return Document(meta=meta, elements=elements)
