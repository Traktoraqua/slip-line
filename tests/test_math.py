"""Tests for Phase 7 text-math: conversion, detection, inline/display."""

from __future__ import annotations

import shutil
import subprocess

import pytest

from pdf2tex.emit import _render_element, render
from pdf2tex.mathconv.detect import (
    build_display_math,
    is_math_span,
    line_ref,
    shift_class,
    spans_to_latex,
)
from pdf2tex.mathconv.unicode_map import convert_math
from pdf2tex.models import (
    DisplayMath,
    DocMeta,
    Document,
    InlineMath,
    Line,
    Page,
    Paragraph,
    Span,
    Text,
)
from pdf2tex.structure import _build_inlines


def _span(text, size=11.0, y=100.0, x=72.0, sup=False, font="Test"):
    flags = (1 << 0) if sup else 0
    return Span(text=text, font=font, size=size, flags=flags,
                bbox=(x, y - size, x + 10 * len(text), y), origin=(x, y))


def _line(spans):
    x0 = min(s.bbox[0] for s in spans)
    x1 = max(s.bbox[2] for s in spans)
    y0 = min(s.bbox[1] for s in spans)
    y1 = max(s.bbox[3] for s in spans)
    return Line(spans=spans, bbox=(x0, y0, x1, y1))


# --------------------------------------------------------------------------- #
# Unicode conversion
# --------------------------------------------------------------------------- #
def test_convert_greek_and_relations():
    assert convert_math("α").strip() == r"\alpha"
    assert r"\le" in convert_math("x ≤ y")
    assert r"\sum" in convert_math("∑")
    assert r"\partial" in convert_math("∂")


def test_convert_super_subscript():
    assert convert_math("x²") == "x^{2}"
    assert convert_math("H₂O") == "H_{2}O"
    assert convert_math("30°") == "30^{\\circ}"


def test_convert_functions():
    assert convert_math("sin") == r"\sin"
    assert convert_math("log") == r"\log"
    # A function name embedded in an identifier is left alone.
    assert convert_math("sinky") == "sinky"


# --------------------------------------------------------------------------- #
# Detection
# --------------------------------------------------------------------------- #
def test_is_math_span():
    assert is_math_span(_span("α"))
    assert is_math_span(_span("x", font="ABCDEF+CMMI10"))  # subset-prefixed font
    assert not is_math_span(_span("hello", font="Helvetica"))


def test_shift_class():
    spans = [_span("x", size=11, y=100), _span("2", size=7, y=96, sup=True)]
    ref_y, ref_size = line_ref(spans)
    assert shift_class(spans[0], ref_y, ref_size) == "normal"
    assert shift_class(spans[1], ref_y, ref_size) == "sup"
    sub = _span("i", size=7, y=103)
    assert shift_class(sub, ref_y, ref_size) == "sub"


def test_spans_to_latex_groups_scripts():
    x = _span("x", size=11, y=100)
    two = _span("2", size=7, y=96, sup=True)
    assert spans_to_latex([(x, "normal"), (two, "sup")]) == "x^{2}"


# --------------------------------------------------------------------------- #
# Inline math in paragraphs
# --------------------------------------------------------------------------- #
def test_inline_greek_in_text():
    line = _line([_span("the angle "), _span("θ", x=150)])
    inlines = _build_inlines([line])
    assert isinstance(inlines[-1], InlineMath)
    assert inlines[-1].latex == r"\theta"
    assert isinstance(inlines[0], Text)


def test_inline_subscript_steals_base():
    base = _span("the v", size=11, y=100)
    sub = _span("i", size=7, y=103, x=120)
    inlines = _build_inlines([_line([base, sub])])
    math = [i for i in inlines if isinstance(i, InlineMath)]
    assert math and math[0].latex == "v_{i}"
    text = "".join(i.text for i in inlines if isinstance(i, Text))
    assert text.strip() == "the"


# --------------------------------------------------------------------------- #
# Display math
# --------------------------------------------------------------------------- #
def _page():
    return Page(number=1, width=595.0, height=842.0, blocks=[])


def test_display_math_numbered():
    spans = [
        _span("E", x=260, y=100), _span("=", x=275, y=100),
        _span("m", x=290, y=100), _span("c", x=305, y=100),
        _span("2", x=315, y=96, size=7, sup=True),
        _span("(1)", x=520, y=100),
    ]
    dm = build_display_math([_line(spans)], _page(), 1, 0)
    assert isinstance(dm, DisplayMath)
    assert dm.latex == "E=mc^{2}"
    assert dm.number == "1"
    assert dm.label == "eq:p1-0"


def test_display_math_rejects_prose():
    spans = [_span(w + " ", x=72 + 40 * i)
             for i, w in enumerate(["This", "ordinary", "centered", "sentence"])]
    assert build_display_math([_line(spans)], _page(), 1, 0) is None


def test_display_math_renders_equation():
    dm = DisplayMath(latex="E=mc^{2}", number="1", label="eq:p1-0",
                     confidence=0.95, page=1)
    out = _render_element(dm)
    assert "\\begin{equation}" in out
    assert "\\label{eq:p1-0}" in out
    assert "% CHECK" not in out


def test_low_confidence_comment():
    dm = DisplayMath(latex="???", number=None, label=None, confidence=0.1, page=4)
    out = _render_element(dm)
    assert "\\[" in out
    assert "% CHECK: low-confidence math, p. 4" in out


# --------------------------------------------------------------------------- #
# Compile
# --------------------------------------------------------------------------- #
def test_math_document_compiles(tmp_path):
    engine = shutil.which("pdflatex")
    if engine is None:
        pytest.skip("no LaTeX engine installed")
    doc = Document(
        meta=DocMeta(source="m.pdf", date="2026-01-01", tool_version="0.1.0"),
        elements=[
            Paragraph(inlines=[Text("The angle "), InlineMath(latex=r"\theta"),
                               Text(" and value "), InlineMath(latex="v_{i}"),
                               Text(".")]),
            DisplayMath(latex="E=mc^{2}", number="1", label="eq:p1-0",
                        confidence=0.95, page=1),
            DisplayMath(latex=r"\sum_{i=1}^{n} x_{i}", confidence=0.9),
        ],
    )
    tex = render(doc, no_title=True)
    out = tmp_path / "out.tex"
    out.write_text(tex, encoding="utf-8")
    proc = subprocess.run(
        [engine, "-interaction=nonstopmode", "-halt-on-error", out.name],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout[-2000:]
