"""Tests for Phase 8: equation-image todos, warning summary and --debug."""

from __future__ import annotations

import datetime
import json
import os

import fitz
import pytest

from pdf2tex import __version__
from pdf2tex.cli import build_parser, convert
from pdf2tex.emit import render
from pdf2tex.extract import extract_document
from pdf2tex.models import DocMeta, TodoPlaceholder
from pdf2tex.structure import build_document
from pdf2tex.zones import filter_zones


def _img_stream(w, h, gray=160):
    pix = fitz.Pixmap(fitz.csRGB, fitz.IRect(0, 0, w, h))
    pix.clear_with(gray)
    return pix.tobytes("png")


@pytest.fixture
def equation_image_pdf(tmp_path):
    """A page with body text and a short isolated 'equation' image."""
    path = tmp_path / "eqimg.pdf"
    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((72, 100), "Text before the displayed equation image.", fontsize=11)
    page.insert_image(fitz.Rect(220, 150, 360, 168), stream=_img_stream(140, 18))
    page.insert_text((72, 210), "Text after the equation image.", fontsize=11)
    doc.save(path)
    doc.close()
    return str(path)


@pytest.fixture
def figure_image_pdf(tmp_path):
    """A page with a large figure image (should be dropped, not a todo)."""
    path = tmp_path / "figimg.pdf"
    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((72, 100), "Body text above a real figure.", fontsize=11)
    page.insert_image(fitz.Rect(120, 200, 400, 520), stream=_img_stream(280, 320))
    doc.save(path)
    doc.close()
    return str(path)


def _convert(path):
    pages = extract_document(path)
    pages, _ = filter_zones(pages)
    meta = DocMeta(source="x.pdf", date=datetime.date.today().isoformat(),
                   tool_version=__version__)
    doc = build_document(pages, meta, detect_title_flag=False)
    return doc, render(doc, no_title=True)


def test_equation_image_becomes_todo(equation_image_pdf):
    doc, tex = _convert(equation_image_pdf)
    todos = [e for e in doc.elements if isinstance(e, TodoPlaceholder)]
    assert len(todos) == 1
    assert "source p. 1" in todos[0].reason
    assert "\\todo[inline]{Equation (image), source p. 1}" in tex
    # Inserted in reading order: between the two paragraphs.
    kinds = [type(e).__name__ for e in doc.elements]
    assert kinds == ["Paragraph", "TodoPlaceholder", "Paragraph"]


def test_large_figure_not_todo(figure_image_pdf):
    doc, _ = _convert(figure_image_pdf)
    assert not any(isinstance(e, TodoPlaceholder) for e in doc.elements)


def _args(argv):
    return build_parser().parse_args(argv)


def test_summary_and_todo_via_cli(equation_image_pdf, capsys):
    convert(_args([equation_image_pdf, "--no-title"]))
    out = capsys.readouterr().out
    assert "Summary: 1 \\todo placeholder(s)" in out


def test_debug_dumps(equation_image_pdf, tmp_path):
    debug_dir = str(tmp_path / "dbg")
    convert(_args([equation_image_pdf, "--no-title", "--debug", debug_dir]))
    for name in ("extract.json", "zones.json", "document.json", "page-1.png"):
        assert os.path.exists(os.path.join(debug_dir, name)), name
    with open(os.path.join(debug_dir, "document.json"), encoding="utf-8") as fh:
        data = json.load(fh)
    types = [e["type"] for e in data["elements"]]
    assert "TodoPlaceholder" in types
    # Page geometry round-trips through JSON.
    with open(os.path.join(debug_dir, "extract.json"), encoding="utf-8") as fh:
        pages = json.load(fh)
    assert pages[0]["number"] == 1 and pages[0]["blocks"]
