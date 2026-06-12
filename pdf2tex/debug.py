"""Debug dumps: per-stage JSON and annotated page PNGs (Phase 8).

Enabled with ``--debug DIR``.  Writes ``extract.json`` / ``zones.json``
(page geometry per stage), ``document.json`` (the element tree) and
``page-N.png`` images with colored boxes per element class, to aid diagnosis
of mis-detected structure.
"""

from __future__ import annotations

import dataclasses
import json
import logging
import os

import fitz

from .lists import match_marker
from .mathconv.detect import has_math_indicator, is_math_span
from .models import Document, Element, Page

log = logging.getLogger(__name__)

# RGB colors (0..1) per line/element class.
_COLORS = {
    "list": (0.9, 0.5, 0.0),
    "math": (0.8, 0.0, 0.8),
    "text": (0.5, 0.5, 0.5),
    "rule": (0.85, 0.1, 0.1),
    "image": (0.0, 0.6, 0.0),
}


def ensure_dir(debug_dir: str) -> None:
    os.makedirs(debug_dir, exist_ok=True)


def _write(debug_dir: str, name: str, obj) -> None:
    path = os.path.join(debug_dir, name)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=2)
    log.debug("debug: wrote %s", path)


def dump_pages(debug_dir: str, stage: str, pages: list[Page]) -> None:
    _write(debug_dir, f"{stage}.json", [dataclasses.asdict(p) for p in pages])


def _element_dict(element: Element) -> dict:
    return {"type": type(element).__name__, **dataclasses.asdict(element)}


def dump_document(debug_dir: str, document: Document) -> None:
    obj = {
        "meta": dataclasses.asdict(document.meta),
        "elements": [_element_dict(e) for e in document.elements],
    }
    _write(debug_dir, "document.json", obj)


def _line_class(line) -> str:
    if match_marker(line.text):
        return "list"
    if any(is_math_span(s) for s in line.spans) or has_math_indicator(line.spans):
        return "math"
    return "text"


def annotate_pages(debug_dir: str, pdf_path: str, pages: list[Page]) -> None:
    """Render each page with colored boxes for lines, rules and images."""
    doc = fitz.open(pdf_path)
    try:
        for page in pages:
            fpage = doc[page.number - 1]
            for block in page.blocks:
                for line in block.lines:
                    color = _COLORS[_line_class(line)]
                    fpage.draw_rect(fitz.Rect(line.bbox), color=color, width=0.8)
            for r in page.rules:
                fpage.draw_line((r.x0, r.y0), (r.x1, r.y1),
                                color=_COLORS["rule"], width=1.0)
            for img in page.images:
                fpage.draw_rect(fitz.Rect(img.bbox), color=_COLORS["image"], width=1.5)
            pix = fpage.get_pixmap(dpi=110)
            pix.save(os.path.join(debug_dir, f"page-{page.number}.png"))
    finally:
        doc.close()
    log.debug("debug: wrote %d annotated page PNG(s)", len(pages))
