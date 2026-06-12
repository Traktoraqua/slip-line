"""Tests for Phase 5 list detection (bullets, numbering, nesting)."""

from __future__ import annotations

import datetime
import shutil
import subprocess

import fitz
import pytest

from pdf2tex import __version__
from pdf2tex.emit import render
from pdf2tex.extract import extract_document
from pdf2tex.lists import build_list_block, match_marker, take_list
from pdf2tex.models import DocMeta, ListBlock, Line, Span
from pdf2tex.structure import _build_inlines, build_document
from pdf2tex.zones import filter_zones


def _line(text, x=72.0, y=100.0):
    sp = Span(text=text, font="Test", size=11.0, flags=0,
              bbox=(x, y - 8, x + 200, y), origin=(x, y))
    return Line(spans=[sp], bbox=(x, y - 8, x + 200, y))


# --------------------------------------------------------------------------- #
# Marker matching
# --------------------------------------------------------------------------- #
def test_match_marker_bullets():
    assert match_marker("• item")[0] == "itemize"
    assert match_marker("- item")[0] == "itemize"
    assert match_marker("* item")[0] == "itemize"


def test_match_marker_numbered():
    assert match_marker("1. item")[0] == "enumerate"
    assert match_marker("(a) item")[0] == "enumerate"
    assert match_marker("iv) item")[0] == "enumerate"


def test_match_marker_negatives():
    assert match_marker("1 Introduction") is None   # heading, no . or )
    assert match_marker("regular text") is None
    assert match_marker("e.g. example") is None      # no space after marker


# --------------------------------------------------------------------------- #
# Building blocks
# --------------------------------------------------------------------------- #
def _rows(lines):
    return [(ln, False) for ln in lines]


def test_flat_bullet_list():
    lines = [_line("• first", y=100), _line("• second", y=114), _line("• third", y=128)]
    rows = _rows(lines)
    end = take_list(rows, 0, leading=14)
    assert end == 3
    lb = build_list_block(rows, 0, end, _build_inlines)
    assert lb.kind == "itemize"
    assert len(lb.items) == 3
    assert lb.items[0][0].text == "first"


def test_continuation_appended():
    lines = [
        _line("• first item text", x=72, y=100),
        _line("continues here", x=86, y=114),   # indented continuation
        _line("• second", x=72, y=128),
    ]
    rows = _rows(lines)
    end = take_list(rows, 0, leading=14)
    assert end == 3
    lb = build_list_block(rows, 0, end, _build_inlines)
    assert len(lb.items) == 2
    assert lb.items[0][0].text == "first item text continues here"


def test_nested_list():
    lines = [
        _line("1. parent", x=72, y=100),
        _line("a) child one", x=92, y=114),   # deeper indent
        _line("b) child two", x=92, y=128),
        _line("2. parent two", x=72, y=142),
    ]
    rows = _rows(lines)
    end = take_list(rows, 0, leading=14)
    lb = build_list_block(rows, 0, end, _build_inlines)
    assert lb.kind == "enumerate"
    # items: parent1, nested ListBlock, parent2
    nested = [it for it in lb.items if isinstance(it, ListBlock)]
    assert len(nested) == 1
    assert len(nested[0].items) == 2
    assert nested[0].items[0][0].text == "child one"


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #
def test_render_nested_list():
    lines = [
        _line("• top", x=72, y=100),
        _line("– sub", x=90, y=114),
        _line("• top2", x=72, y=128),
    ]
    rows = _rows(lines)
    lb = build_list_block(rows, 0, take_list(rows, 0, 14), _build_inlines)
    from pdf2tex.emit import _render_list
    out = _render_list(lb)
    assert out.count("\\begin{itemize}") == 2
    assert out.count("\\end{itemize}") == 2
    assert "\\item top" in out
    assert "\\item sub" in out


# --------------------------------------------------------------------------- #
# Integration through a real PDF
# --------------------------------------------------------------------------- #
@pytest.fixture
def list_pdf(tmp_path):
    path = tmp_path / "list.pdf"
    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((72, 100), "Intro paragraph before the list here.", fontsize=11)
    page.insert_text((72, 140), "1. First numbered point", fontsize=11)
    page.insert_text((72, 158), "2. Second numbered point", fontsize=11)
    page.insert_text((72, 176), "3. Third numbered point", fontsize=11)
    page.insert_text((72, 210), "Closing paragraph after the list.", fontsize=11)
    doc.save(path)
    doc.close()
    return str(path)


def _convert(path):
    pages = extract_document(path)
    pages, _ = filter_zones(pages)
    meta = DocMeta(source="list.pdf", date=datetime.date.today().isoformat(),
                   tool_version=__version__)
    doc = build_document(pages, meta, detect_title_flag=False)
    return doc, render(doc, no_title=True)


def test_list_pdf_structure(list_pdf):
    doc, tex = _convert(list_pdf)
    assert sum(isinstance(e, ListBlock) for e in doc.elements) == 1
    assert "\\begin{enumerate}" in tex
    assert "\\item First numbered point" in tex
    assert "Intro paragraph before the list here." in tex
    assert "Closing paragraph after the list." in tex


def test_list_pdf_compiles(list_pdf, tmp_path):
    engine = shutil.which("pdflatex")
    if engine is None:
        pytest.skip("no LaTeX engine installed")
    _, tex = _convert(list_pdf)
    out = tmp_path / "out.tex"
    out.write_text(tex, encoding="utf-8")
    proc = subprocess.run(
        [engine, "-interaction=nonstopmode", "-halt-on-error", out.name],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout[-2000:]
