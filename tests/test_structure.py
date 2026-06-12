"""Tests for Phase 3 structure: headings, dehyphenation, inline styles."""

from __future__ import annotations

import datetime

import fitz
import pytest

from pdf2tex import __version__
from pdf2tex.emit import render
from pdf2tex.extract import extract_document
from pdf2tex.models import Bold, DocMeta, Heading, Italic, Line, Paragraph, Span, Text
from pdf2tex.structure import (
    _body_size,
    _build_inlines,
    _heading_from_group,
    build_document,
)
from pdf2tex.zones import filter_zones


# --------------------------------------------------------------------------- #
# Unit tests on synthetic spans/lines
# --------------------------------------------------------------------------- #
def _span(text, size=11.0, bold=False, italic=False, x=72.0, y=100.0):
    flags = (1 << 4 if bold else 0) | (1 << 1 if italic else 0)
    return Span(text=text, font="Test-Bold" if bold else "Test", size=size,
                flags=flags, bbox=(x, y - 8, x + 50, y), origin=(x, y))


def _line(spans, y=100.0):
    return Line(spans=spans, bbox=(spans[0].bbox[0], y - 8, spans[-1].bbox[2], y))


def test_inline_styles_mixed_run():
    line = _line([
        _span("This is "),
        _span("bold", bold=True),
        _span(" and "),
        _span("italic", italic=True),
        _span(" text."),
    ])
    inlines = _build_inlines([line])
    kinds = [(type(i).__name__, i.text) for i in inlines]
    assert kinds == [
        ("Text", "This is "),
        ("Bold", "bold"),
        ("Text", " and "),
        ("Italic", "italic"),
        ("Text", " text."),
    ]


def test_dehyphenation_across_lines():
    l1 = _line([_span("This is infor-")], y=100.0)
    l2 = _line([_span("mation that continues")], y=114.0)
    inlines = _build_inlines([l1, l2])
    assert len(inlines) == 1
    assert isinstance(inlines[0], Text)
    assert inlines[0].text == "This is information that continues"


def test_line_join_adds_space_without_hyphen():
    l1 = _line([_span("first line")], y=100.0)
    l2 = _line([_span("second line")], y=114.0)
    inlines = _build_inlines([l1, l2])
    assert inlines[0].text == "first line second line"


def test_numbered_heading_levels():
    body = 11.0
    h1 = _heading_from_group([_line([_span("1 Introduction", size=14, bold=True)])], body)
    h2 = _heading_from_group([_line([_span("2.3.1 Scope", size=12, bold=True)])], body)
    assert isinstance(h1, Heading) and h1.level == 1 and h1.text == "Introduction"
    assert h1.number == "1"
    assert h2.level == 3 and h2.text == "Scope" and h2.number == "2.3.1"


def test_numbered_list_item_is_not_heading():
    # "1. Foo" (dot before space) must not be treated as a heading.
    assert _heading_from_group([_line([_span("1. Foo bar", size=14, bold=True)])], 11.0) is None


def test_unnumbered_bold_heading():
    h = _heading_from_group([_line([_span("Summary", size=15, bold=True)])], 11.0)
    assert isinstance(h, Heading) and h.starred and h.number is None
    # Plain body-sized text is never an unnumbered heading.
    assert _heading_from_group([_line([_span("Summary", size=11)])], 11.0) is None


def test_body_size_is_char_weighted_mode():
    pages = []  # build a fake page via models
    from pdf2tex.models import Block, Page
    big = _span("X", size=22)
    body = [_span("word " * 5, size=11) for _ in range(3)]
    pg = Page(number=1, width=595, height=842,
              blocks=[Block(lines=[_line([big])], bbox=(0, 0, 1, 1)),
                      Block(lines=[_line(body)], bbox=(0, 0, 1, 1))])
    assert _body_size([pg]) == 11.0


# --------------------------------------------------------------------------- #
# Integration test through a real PDF
# --------------------------------------------------------------------------- #
@pytest.fixture
def report_pdf(tmp_path):
    path = tmp_path / "report.pdf"
    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((72, 80), "Main Report Title", fontsize=22)            # title
    page.insert_text((72, 140), "1 Introduction", fontsize=14, fontname="hebo")
    page.insert_text((72, 175),
                     "This introduction has enough body words to set the body font size here.",
                     fontsize=11)
    page.insert_text((72, 220), "2.1 Methods", fontsize=12, fontname="hebo")
    page.insert_text((72, 255),
                     "The methods paragraph also contains a number of plain body words.",
                     fontsize=11)
    doc.save(path)
    doc.close()
    return str(path)


def _convert(path):
    pages = extract_document(path)
    pages, _ = filter_zones(pages)
    meta = DocMeta(source="report.pdf", date=datetime.date.today().isoformat(),
                   tool_version=__version__)
    doc = build_document(pages, meta)
    return doc, render(doc)


def test_report_headings_rendered(report_pdf):
    doc, tex = _convert(report_pdf)
    assert "\\title{Main Report Title}" in tex
    assert "\\section{Introduction}" in tex
    assert "\\subsection{Methods}" in tex
    # Source numbers are stripped; LaTeX renumbers.
    assert "1 Introduction" not in tex
    assert "set the body font size here" in tex


def test_report_compiles(report_pdf):
    import shutil
    import subprocess

    engine = shutil.which("pdflatex")
    if engine is None:
        pytest.skip("no LaTeX engine installed")
    _, tex = _convert(report_pdf)
    d = report_pdf.rsplit("/", 1)[0]
    with open(d + "/out.tex", "w", encoding="utf-8") as fh:
        fh.write(tex)
    proc = subprocess.run([engine, "-interaction=nonstopmode", "-halt-on-error", "out.tex"],
                          cwd=d, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout[-2000:]
