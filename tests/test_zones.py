"""Tests for Phase 2 zone filtering (headers/footers/page numbers/figures)."""

from __future__ import annotations

import fitz
import pytest

from pdf2tex.extract import extract_document
from pdf2tex.zones import _normalize, filter_zones


@pytest.fixture
def multipage_pdf(tmp_path):
    """3-page PDF: repeated header, per-page footer page number, unique body."""
    path = tmp_path / "multi.pdf"
    doc = fitz.open()
    for i in range(1, 4):
        page = doc.new_page()  # default A4 ~ 595 x 842
        page.insert_text((72, 24), "CONFIDENTIAL DRAFT", fontsize=9)   # top band
        page.insert_text((72, 400), f"Unique body text for page {i} here.",
                          fontsize=11)                                  # middle
        page.insert_text((290, 820), str(i), fontsize=9)               # bottom band
    doc.save(path)
    doc.close()
    return str(path)


@pytest.fixture
def figure_pdf(tmp_path):
    """Single page with an embedded image and a figure caption beneath it."""
    path = tmp_path / "fig.pdf"
    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((72, 120), "Body paragraph above the figure.", fontsize=11)
    pix = fitz.Pixmap(fitz.csRGB, fitz.IRect(0, 0, 120, 120))
    pix.clear_with(180)
    page.insert_image(fitz.Rect(72, 300, 240, 460), stream=pix.tobytes("png"))
    page.insert_text((72, 480), "Figure 1: A grey square.", fontsize=10)
    page.insert_text((72, 520), "Body paragraph below the figure.", fontsize=11)
    doc.save(path)
    doc.close()
    return str(path)


def test_normalize_masks_digits():
    assert _normalize("Page 12") == "page #"
    assert _normalize("Page  3") == _normalize("Page 99")


def test_header_and_pagenumber_removed(multipage_pdf):
    pages = extract_document(multipage_pdf)
    pages, stats = filter_zones(pages)

    body = "\n".join(line.text for p in pages for line in p.lines)
    assert "CONFIDENTIAL DRAFT" not in body
    assert "Unique body text for page 1 here." in body
    assert stats.headers_footers >= 3
    # Page numbers caught either as repeats ("#") or by the numeric rule.
    assert stats.headers_footers + stats.page_numbers >= 6


def test_figure_caption_removed(figure_pdf):
    pages = extract_document(figure_pdf)
    pages, stats = filter_zones(pages)

    body = "\n".join(line.text for p in pages for line in p.lines)
    assert "Figure 1: A grey square." not in body
    assert "Body paragraph above the figure." in body
    assert "Body paragraph below the figure." in body
    assert stats.figures == 1
    assert stats.captions == 1


def test_keep_captions_flag(figure_pdf):
    pages = extract_document(figure_pdf)
    pages, stats = filter_zones(pages, keep_captions=True)
    body = "\n".join(line.text for p in pages for line in p.lines)
    assert "Figure 1: A grey square." in body
    assert stats.captions == 0


def test_single_page_keeps_band_text(tmp_path):
    """Repeat analysis must not strip a lone page's header-band text."""
    path = tmp_path / "one.pdf"
    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((72, 24), "Report Heading In Band", fontsize=12)
    page.insert_text((72, 400), "Body text.", fontsize=11)
    doc.save(path)
    doc.close()

    pages, stats = filter_zones(extract_document(str(path)))
    body = "\n".join(line.text for p in pages for line in p.lines)
    assert "Report Heading In Band" in body
    assert stats.headers_footers == 0
