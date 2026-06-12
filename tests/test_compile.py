"""Compile check: emitted .tex must build under pdflatex when available.

Skipped automatically when no LaTeX engine is installed (e.g. dev machines
without TeX), so the rest of the suite stays runnable everywhere.
"""

from __future__ import annotations

import datetime
import shutil
import subprocess

import pytest

from pdf2tex import __version__
from pdf2tex.emit import render
from pdf2tex.extract import extract_document
from pdf2tex.models import DocMeta
from pdf2tex.structure import build_document

_ENGINE = shutil.which("pdflatex")

pytestmark = pytest.mark.skipif(_ENGINE is None, reason="no LaTeX engine installed")


def test_output_compiles(simple_pdf, tmp_path):
    pages = extract_document(simple_pdf)
    meta = DocMeta(source="simple.pdf", date=datetime.date.today().isoformat(),
                   tool_version=__version__)
    tex = render(build_document(pages, meta))

    tex_path = tmp_path / "out.tex"
    tex_path.write_text(tex, encoding="utf-8")

    proc = subprocess.run(
        [_ENGINE, "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout[-2000:]
    assert (tmp_path / "out.pdf").exists()
