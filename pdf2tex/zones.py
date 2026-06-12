"""Zone filtering: drop running headers/footers, page numbers and figures.

Operates on the list of :class:`~pdf2tex.models.Page` objects produced by
:mod:`pdf2tex.extract`, *before* structure assembly.  Three removals:

* **Headers/footers** — lines in the top/bottom ~8 % band whose digit-masked
  text repeats on ≥ 60 % of pages (so page-number variation does not defeat
  the match).
* **Page numbers** — isolated short numeric/roman lines in the margin band.
* **Figures** — figure image blocks are dropped, and adjacent caption lines
  (``Figure 1 …``) within ~2 line heights are removed unless
  ``keep_captions`` is set.  Small inline images on a text line are kept as
  equation-image candidates for the math phase.
"""

from __future__ import annotations

import logging
import math
import re
import statistics
from dataclasses import dataclass

from .models import Block, Line, Page

log = logging.getLogger(__name__)

BAND_FRACTION = 0.08          # top/bottom band as a fraction of page height
REPEAT_FRACTION = 0.60        # share of pages a header/footer must appear on
CAPTION_ADJACENCY = 2.0       # caption within this many line heights of figure
INLINE_IMAGE_FACTOR = 1.8     # image taller than this × line height ⇒ figure

_CAPTION_RE = re.compile(r"^(Figure|Fig\.?|Figur)\s+[\dA-Z][-.\d]*")
_DIGITS_RE = re.compile(r"\d+")
_WS_RE = re.compile(r"\s+")
# Pure page-number lines: "12", "iv", "- 5 -", "Page 7", "Side 7".
_PAGENO_RE = re.compile(
    r"^(?:page|side|p\.?)?\s*[-–—]?\s*(?:\d+|[ivxlcdm]+)\s*[-–—]?$",
    re.IGNORECASE,
)


@dataclass
class ZoneStats:
    headers_footers: int = 0
    page_numbers: int = 0
    figures: int = 0
    captions: int = 0


def _line_center_y(line: Line) -> float:
    return (line.bbox[1] + line.bbox[3]) / 2.0


def _in_band(line: Line, height: float) -> bool:
    cy = _line_center_y(line)
    return cy <= height * BAND_FRACTION or cy >= height * (1.0 - BAND_FRACTION)


def _normalize(text: str) -> str:
    """Collapse whitespace and mask digit runs so page numbers compare equal."""
    return _DIGITS_RE.sub("#", _WS_RE.sub(" ", text.strip())).lower()


def _median_line_height(pages: list[Page]) -> float:
    heights = [
        line.bbox[3] - line.bbox[1]
        for page in pages
        for line in page.lines
        if line.bbox[3] > line.bbox[1]
    ]
    return statistics.median(heights) if heights else 12.0


# --------------------------------------------------------------------------- #
# Header / footer repeat analysis
# --------------------------------------------------------------------------- #
def _repeated_band_keys(pages: list[Page]) -> set[str]:
    """Digit-masked keys of band lines that repeat on enough pages."""
    pages_per_key: dict[str, set[int]] = {}
    for page in pages:
        seen_on_page: set[str] = set()
        for line in page.lines:
            if not _in_band(line, page.height):
                continue
            key = _normalize(line.text)
            if key and key not in seen_on_page:
                seen_on_page.add(key)
                pages_per_key.setdefault(key, set()).add(page.number)

    required = max(2, math.ceil(REPEAT_FRACTION * len(pages)))
    return {key for key, pgs in pages_per_key.items() if len(pgs) >= required}


def _is_page_number(line: Line, height: float) -> bool:
    return _in_band(line, height) and bool(_PAGENO_RE.match(line.text.strip()))


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def _classify_images(page: Page, line_h: float):
    """Split image blocks into (figures, inline_candidates)."""
    figures, inline = [], []
    for img in page.images:
        ih = img.bbox[3] - img.bbox[1]
        overlaps_line = any(
            not (img.bbox[3] < ln.bbox[1] or img.bbox[1] > ln.bbox[3])
            for ln in page.lines
        )
        if ih <= INLINE_IMAGE_FACTOR * line_h and overlaps_line:
            inline.append(img)
        else:
            figures.append(img)
    return figures, inline


def _caption_near_figure(line: Line, fig_bbox, line_h: float) -> bool:
    margin = CAPTION_ADJACENCY * line_h
    below = fig_bbox[3] - margin <= line.bbox[1] <= fig_bbox[3] + margin
    above = fig_bbox[1] - margin <= line.bbox[3] <= fig_bbox[1] + margin
    return (below or above) and bool(_CAPTION_RE.match(line.text.strip()))


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def filter_zones(
    pages: list[Page],
    keep_captions: bool = False,
) -> tuple[list[Page], ZoneStats]:
    """Filter headers/footers, page numbers and figures from *pages* in place."""
    stats = ZoneStats()
    line_h = _median_line_height(pages)
    repeated = _repeated_band_keys(pages)

    for page in pages:
        figures, inline = _classify_images(page, line_h)
        stats.figures += len(figures)
        page.images = inline  # keep only equation-image candidates

        kept_blocks: list[Block] = []
        for block in page.blocks:
            kept_lines: list[Line] = []
            for line in block.lines:
                if _in_band(line, page.height) and _normalize(line.text) in repeated:
                    stats.headers_footers += 1
                    continue
                if _is_page_number(line, page.height):
                    stats.page_numbers += 1
                    continue
                if (
                    not keep_captions
                    and figures
                    and any(_caption_near_figure(line, f.bbox, line_h) for f in figures)
                ):
                    stats.captions += 1
                    continue
                kept_lines.append(line)
            if kept_lines:
                kept_blocks.append(Block(lines=kept_lines, bbox=block.bbox))
        page.blocks = kept_blocks

    log.info(
        "Zones: dropped %d header/footer line(s), %d page number(s), "
        "%d figure(s), %d caption(s)",
        stats.headers_footers,
        stats.page_numbers,
        stats.figures,
        stats.captions,
    )
    return pages, stats
