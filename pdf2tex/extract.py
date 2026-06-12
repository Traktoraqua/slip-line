"""Raw extraction from a PDF using PyMuPDF (``fitz``).

Produces :class:`~pdf2tex.models.Page` objects carrying text blocks (with
fonts, sizes, flags and geometry), vector rule segments (for later table
detection) and image-block bounding boxes (recorded, never extracted).

Text is Unicode-normalised (NFC), common ligatures are expanded and
zero-width characters are stripped so downstream stages see clean text.
"""

from __future__ import annotations

import logging
import unicodedata

import fitz

from .models import Block, ImageBlock, Line, Page, Rule, Span

log = logging.getLogger(__name__)

# Ligatures PyMuPDF may hand back as single code points.
_LIGATURES = {
    "ﬀ": "ff",
    "ﬁ": "fi",
    "ﬂ": "fl",
    "ﬃ": "ffi",
    "ﬄ": "ffl",
    "ﬅ": "st",
    "ﬆ": "st",
}

# Zero-width / formatting code points to drop outright.
_ZERO_WIDTH = dict.fromkeys(
    map(ord, "​‌‍‎‏﻿"), None
)

# A segment is treated as a rule only if its off-axis extent is below this.
_RULE_THICKNESS = 2.0


def normalize_text(text: str) -> str:
    """NFC-normalise, expand ligatures and strip zero-width characters."""
    for src, dst in _LIGATURES.items():
        if src in text:
            text = text.replace(src, dst)
    text = text.translate(_ZERO_WIDTH)
    return unicodedata.normalize("NFC", text)


def _parse_pages_arg(spec: str, page_count: int) -> list[int]:
    """Parse a ``--pages`` spec like ``"1-12,15"`` into 0-based indices."""
    indices: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo_s, hi_s = part.split("-", 1)
            lo, hi = int(lo_s), int(hi_s)
        else:
            lo = hi = int(part)
        for n in range(lo, hi + 1):
            if 1 <= n <= page_count and (n - 1) not in indices:
                indices.append(n - 1)
    return indices


def _spans_from_line(line_dict: dict) -> list[Span]:
    spans: list[Span] = []
    for sp in line_dict.get("spans", []):
        text = normalize_text(sp.get("text", ""))
        if not text:
            continue
        spans.append(
            Span(
                text=text,
                font=sp.get("font", ""),
                size=float(sp.get("size", 0.0)),
                flags=int(sp.get("flags", 0)),
                bbox=tuple(sp.get("bbox", (0, 0, 0, 0))),
                origin=tuple(sp.get("origin", (0, 0))),
                color=int(sp.get("color", 0)),
            )
        )
    return spans


def _rules_from_drawings(drawings: list[dict]) -> list[Rule]:
    rules: list[Rule] = []
    for path in drawings:
        for item in path.get("items", []):
            kind = item[0]
            if kind == "l":  # line segment: ("l", p1, p2)
                p1, p2 = item[1], item[2]
                if min(abs(p2.x - p1.x), abs(p2.y - p1.y)) <= _RULE_THICKNESS:
                    rules.append(Rule(p1.x, p1.y, p2.x, p2.y))
            elif kind == "re":  # rectangle: emit its four edges
                r = item[1]
                rules.append(Rule(r.x0, r.y0, r.x1, r.y0))
                rules.append(Rule(r.x0, r.y1, r.x1, r.y1))
                rules.append(Rule(r.x0, r.y0, r.x0, r.y1))
                rules.append(Rule(r.x1, r.y0, r.x1, r.y1))
    return rules


def extract_page(page: "fitz.Page", number: int) -> Page:
    raw = page.get_text("dict")
    blocks: list[Block] = []
    images: list[ImageBlock] = []

    for block in raw.get("blocks", []):
        if block.get("type", 0) == 1:  # image block
            images.append(ImageBlock(bbox=tuple(block.get("bbox", (0, 0, 0, 0)))))
            continue
        lines: list[Line] = []
        for line_dict in block.get("lines", []):
            spans = _spans_from_line(line_dict)
            if spans:
                lines.append(Line(spans=spans, bbox=tuple(line_dict.get("bbox", (0, 0, 0, 0)))))
        if lines:
            blocks.append(Block(lines=lines, bbox=tuple(block.get("bbox", (0, 0, 0, 0)))))

    try:
        rules = _rules_from_drawings(page.get_drawings())
    except Exception as exc:  # drawings parsing is best-effort
        log.debug("get_drawings failed on page %d: %s", number, exc)
        rules = []

    return Page(
        number=number,
        width=float(raw.get("width", page.rect.width)),
        height=float(raw.get("height", page.rect.height)),
        blocks=blocks,
        rules=rules,
        images=images,
    )


def extract_document(path: str, pages: str | None = None) -> list[Page]:
    """Open *path* and extract the requested pages (default: all)."""
    doc = fitz.open(path)
    try:
        if pages:
            indices = _parse_pages_arg(pages, doc.page_count)
        else:
            indices = list(range(doc.page_count))
        result = [extract_page(doc[i], number=i + 1) for i in indices]
    finally:
        doc.close()
    log.info("Extracted %d page(s) from %s", len(result), path)
    return result
