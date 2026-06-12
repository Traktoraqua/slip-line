"""Shared pytest fixtures: synthetic PDFs generated with PyMuPDF."""

from __future__ import annotations

import fitz
import pytest


@pytest.fixture
def simple_pdf(tmp_path):
    """A two-paragraph PDF with a large title line on page 1."""
    path = tmp_path / "simple.pdf"
    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((72, 90), "A Short Technical Report", fontsize=22)
    page.insert_text(
        (72, 140),
        "This is the first paragraph of the body text. It has a few words.",
        fontsize=11,
    )
    page.insert_text(
        (72, 180),
        "Costs are 5 % of total & include a value_of_x term.",
        fontsize=11,
    )
    doc.save(path)
    doc.close()
    return str(path)
