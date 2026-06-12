"""Tests for Phase 9 math-OCR fallback (with a fake backend)."""

from __future__ import annotations

import os

import fitz
import pytest

from pdf2tex.emit import _render_element
from pdf2tex.mathconv import ocr
from pdf2tex.mathconv.ocr import OcrResult, render_crop, validate
from pdf2tex.models import DisplayMath, Document, DocMeta, TodoPlaceholder


# --------------------------------------------------------------------------- #
# Validation (pure)
# --------------------------------------------------------------------------- #
BBOX = (100.0, 100.0, 240.0, 130.0)  # 140 x 30


def test_validate_accepts_reasonable():
    assert validate(r"E = mc^2", BBOX)
    assert validate(r"\frac{a}{b}", BBOX)
    assert validate(r"\left( x \right)", BBOX)


def test_validate_rejects_bad():
    assert not validate("", BBOX)
    assert not validate("   ", BBOX)
    assert not validate(r"\frac{a}{b", BBOX)          # unbalanced braces
    assert not validate(r"\left( x", BBOX)            # unbalanced \left/\right
    assert not validate("x" * 5000, BBOX)             # implausibly long


# --------------------------------------------------------------------------- #
# Crop rendering (no ML needed)
# --------------------------------------------------------------------------- #
@pytest.fixture
def image_pdf(tmp_path):
    path = tmp_path / "crop.pdf"
    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((110, 120), "E = m c^2", fontsize=12)
    doc.save(path)
    doc.close()
    return str(path)


def test_render_crop_returns_png(image_pdf):
    png = render_crop(image_pdf, 1, (100, 100, 240, 130), dpi=120)
    assert png[:8] == b"\x89PNG\r\n\x1a\n"


# --------------------------------------------------------------------------- #
# apply_ocr with a fake backend
# --------------------------------------------------------------------------- #
class FakeBackend:
    name = "fake"

    def __init__(self, latex):
        self._latex = latex

    def predict(self, png_bytes):
        return OcrResult(latex=self._latex)


def _doc_with_todo():
    return Document(
        meta=DocMeta(source="x", date="2026-01-01", tool_version="0.1.0"),
        elements=[
            TodoPlaceholder(reason="Equation (image), source p. 1", page=1,
                            bbox=(100, 100, 240, 130)),
        ],
    )


def test_apply_ocr_accepts_and_replaces(image_pdf):
    doc = _doc_with_todo()
    accepted, rejected = ocr.apply_ocr(doc, FakeBackend(r"E = mc^2"), image_pdf)
    assert (accepted, rejected) == (1, 0)
    el = doc.elements[0]
    assert isinstance(el, DisplayMath)
    assert el.latex == "E = mc^2"
    assert el.comment == "OCR (fake): verify, p. 1"
    assert "% OCR (fake): verify, p. 1" in _render_element(el)


def test_apply_ocr_rejects_keeps_todo(image_pdf):
    doc = _doc_with_todo()
    accepted, rejected = ocr.apply_ocr(doc, FakeBackend(r"\frac{a}{b"), image_pdf)
    assert (accepted, rejected) == (0, 1)
    assert isinstance(doc.elements[0], TodoPlaceholder)


def test_apply_ocr_debug_saves_pair(image_pdf, tmp_path):
    doc = _doc_with_todo()
    dbg = str(tmp_path / "d")
    os.makedirs(dbg)
    ocr.apply_ocr(doc, FakeBackend(r"E = mc^2"), image_pdf, debug_dir=dbg)
    assert os.path.exists(os.path.join(dbg, "ocr-p1-1.png"))
    assert os.path.exists(os.path.join(dbg, "ocr-p1-1.tex"))


# --------------------------------------------------------------------------- #
# Backend factory
# --------------------------------------------------------------------------- #
def test_get_backend_unknown():
    with pytest.raises(ValueError):
        ocr.get_backend("nope")


def test_get_backend_missing_dependency():
    # The ML libs are not installed in CI; loading must raise a clear hint.
    with pytest.raises(ImportError, match="pip install pdf2tex"):
        ocr.get_backend("pix2tex")
