"""End-to-end and unit tests for the Phase 1 skeleton."""

from __future__ import annotations

import datetime

from pdf2tex import __version__
from pdf2tex.charmap import escape_text
from pdf2tex.emit import render
from pdf2tex.extract import _parse_pages_arg, extract_document, normalize_text
from pdf2tex.models import DocMeta
from pdf2tex.structure import build_document


def _convert(path, no_title=False):
    pages = extract_document(path)
    meta = DocMeta(
        source="simple.pdf",
        date=datetime.date.today().isoformat(),
        tool_version=__version__,
    )
    doc = build_document(pages, meta, detect_title_flag=not no_title)
    return doc, render(doc, no_title=no_title)


def test_end_to_end_produces_document(simple_pdf):
    doc, tex = _convert(simple_pdf)
    assert "\\documentclass[11pt,a4paper]{article}" in tex
    assert "\\begin{document}" in tex
    assert "\\end{document}" in tex
    assert "first paragraph of the body text" in tex


def test_title_detected_and_maketitle(simple_pdf):
    doc, tex = _convert(simple_pdf)
    assert doc.meta.title == "A Short Technical Report"
    assert "\\title{A Short Technical Report}" in tex
    assert "\\maketitle" in tex
    # The title line must not also appear as a body paragraph.
    assert tex.count("A Short Technical Report") == 1


def test_no_title_flag(simple_pdf):
    _, tex = _convert(simple_pdf, no_title=True)
    assert "\\maketitle" not in tex


def test_reserved_characters_escaped(simple_pdf):
    _, tex = _convert(simple_pdf)
    assert "5 \\% of total \\& include a value\\_of\\_x" in tex


def test_escape_text_unit():
    assert escape_text("a & b") == "a \\& b"
    assert escape_text("100%") == "100\\%"
    assert escape_text("a_b") == "a\\_b"
    assert escape_text("x^2") == "x\\textasciicircum{}2"
    assert escape_text("\\foo") == "\\textbackslash{}foo"


def test_normalize_ligatures_and_zero_width():
    assert normalize_text("ﬁle") == "file"
    assert normalize_text("a​b") == "ab"


def test_parse_pages_arg():
    assert _parse_pages_arg("1-3,5", 10) == [0, 1, 2, 4]
    assert _parse_pages_arg("2", 10) == [1]
    # Out-of-range pages are clamped away.
    assert _parse_pages_arg("8-20", 10) == [7, 8, 9]


def test_utf8_passthrough(tmp_path):
    import fitz

    path = tmp_path / "nordic.pdf"
    d = fitz.open()
    pg = d.new_page()
    pg.insert_text((72, 100), "Måling av størrelse på væske", fontsize=11)
    d.save(path)
    d.close()
    _, tex = _convert(str(path), no_title=True)
    assert "Måling av størrelse på væske" in tex
