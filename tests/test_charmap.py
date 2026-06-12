"""Tests for Phase 4 character handling (escaping + typographic + æøå)."""

from __future__ import annotations

import datetime
import shutil
import subprocess

import fitz
import pytest

from pdf2tex import __version__
from pdf2tex.charmap import escape_text
from pdf2tex.emit import render
from pdf2tex.extract import extract_document
from pdf2tex.models import DocMeta
from pdf2tex.structure import build_document
from pdf2tex.zones import filter_zones


def test_reserved_characters():
    assert escape_text("# $ % & _ { }") == r"\# \$ \% \& \_ \{ \}"
    assert escape_text("~") == r"\textasciitilde{}"
    assert escape_text("^") == r"\textasciicircum{}"
    assert escape_text("\\") == r"\textbackslash{}"


def test_smart_quotes():
    assert escape_text("“quoted”") == "``quoted''"
    assert escape_text("‘x’") == "`x'"
    assert escape_text("don’t") == "don't"


def test_dashes_and_ellipsis():
    assert escape_text("a–b") == "a--b"          # en dash
    assert escape_text("a—b") == "a---b"         # em dash
    assert escape_text("wait…") == r"wait\dots{}"
    assert escape_text("−5") == "-5"             # minus sign


def test_spaces():
    assert escape_text("a b") == "a~b"           # nbsp → tie
    assert escape_text("a b") == "a~b"           # narrow nbsp → tie
    assert escape_text("a b") == "a b"           # thin space → space
    assert escape_text("soft­hyphen") == "softhyphen"  # soft hyphen dropped


def test_nordic_passthrough():
    s = "Måling av størrelse på væske, é ü ñ"
    assert escape_text(s) == s


def test_no_double_escape():
    # The backslash we introduce for & must not be re-escaped.
    assert escape_text("&") == r"\&"
    assert "\\\\" not in escape_text("&%#_")


@pytest.fixture
def special_pdf(tmp_path):
    path = tmp_path / "special.pdf"
    doc = fitz.open()
    page = doc.new_page()
    # NB: PyMuPDF's base-14 insert_text cannot encode en-dash/curly-quote/
    # ellipsis glyphs, so those substitutions are exercised by the unit tests
    # above; here we round-trip æøå and reserved characters through a real PDF.
    page.insert_text((72, 100), "Måling av størrelse på væske.", fontsize=12)
    page.insert_text((72, 130), "Cost is 5 % of A&B with x_i terms.", fontsize=12)
    doc.save(path)
    doc.close()
    return str(path)


def _convert(path):
    pages = extract_document(path)
    pages, _ = filter_zones(pages)
    meta = DocMeta(source="special.pdf", date=datetime.date.today().isoformat(),
                   tool_version=__version__)
    return render(build_document(pages, meta, detect_title_flag=False))


def test_special_pdf_escaped(special_pdf):
    tex = _convert(special_pdf)
    assert "Måling av størrelse på væske." in tex
    assert "5 \\% of A\\&B with x\\_i terms." in tex


def test_special_pdf_compiles(special_pdf, tmp_path):
    engine = shutil.which("pdflatex")
    if engine is None:
        pytest.skip("no LaTeX engine installed")
    tex = _convert(special_pdf)
    out = tmp_path / "special.tex"
    out.write_text(tex, encoding="utf-8")
    proc = subprocess.run(
        [engine, "-interaction=nonstopmode", "-halt-on-error", out.name],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout[-2000:]
    assert (tmp_path / "special.pdf").exists()
