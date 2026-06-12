"""Document structure assembly.

Phase 1 is deliberately naive: every text block becomes one paragraph, with
intra-block lines joined (dehyphenating ``word-`` + lowercase continuations).
A simple title heuristic promotes the largest text on page 1.

Later phases replace this with body-font analysis, headings, lists, tables
and math; the :func:`build_document` signature stays stable.
"""

from __future__ import annotations

import logging

from .models import Block, DocMeta, Document, Page, Paragraph, Text

log = logging.getLogger(__name__)


def _join_block_lines(block: Block) -> str:
    """Join a block's lines into one paragraph string.

    Trailing ``-`` followed by a lowercase continuation is dehyphenated;
    otherwise lines are joined with a single space.
    """
    out = ""
    for line in block.lines:
        piece = line.text.strip()
        if not piece:
            continue
        if not out:
            out = piece
            continue
        if out.endswith("-") and piece[:1].islower():
            out = out[:-1] + piece
        else:
            out = out + " " + piece
    return out


def _largest_span(pages: list[Page]):
    """Return the span with the maximum font size on the first page, if any."""
    if not pages:
        return None
    best = None
    for span in pages[0].spans:
        if span.text.strip() and (best is None or span.size > best.size):
            best = span
    return best


def detect_title(pages: list[Page]) -> str | None:
    """Heuristic title: the line carrying the largest text on page 1."""
    span = _largest_span(pages)
    if span is None:
        return None
    # Only treat clearly-larger text as a title.
    sizes = [s.size for p in pages for s in p.spans if s.text.strip()]
    if not sizes:
        return None
    body = sorted(sizes)[len(sizes) // 2]  # median size
    if span.size <= body * 1.15:
        return None
    for line in pages[0].lines:
        if span in line.spans:
            return line.text.strip()
    return span.text.strip()


def build_document(
    pages: list[Page],
    meta: DocMeta,
    detect_title_flag: bool = True,
) -> Document:
    """Assemble pages into a :class:`Document` of naive paragraphs."""
    title = detect_title(pages) if detect_title_flag else None
    if title:
        meta.title = title

    elements: list = []
    for page in pages:
        for block in page.blocks:
            text = _join_block_lines(block)
            if not text:
                continue
            if title and page.number == 1 and text.strip() == title:
                continue  # already consumed as the title
            elements.append(Paragraph(inlines=[Text(text)]))

    log.info("Built document with %d paragraph(s)", len(elements))
    return Document(meta=meta, elements=elements)
