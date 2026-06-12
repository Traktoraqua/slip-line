"""Tests for Phase 6 ruled-table detection and rendering."""

from __future__ import annotations

import datetime
import shutil
import subprocess

import fitz
import pytest

from pdf2tex import __version__
from pdf2tex.emit import _render_table, render
from pdf2tex.extract import extract_document
from pdf2tex.models import DocMeta, Table
from pdf2tex.structure import build_document
from pdf2tex.tables import detect_tables
from pdf2tex.zones import filter_zones


def _grid_pdf(path, cols_x, rows_y, cells, interior_top=True, caption=None):
    doc = fitz.open()
    page = doc.new_page()
    for y in rows_y:
        page.draw_line((cols_x[0], y), (cols_x[-1], y))
    # Outer verticals.
    page.draw_line((cols_x[0], rows_y[0]), (cols_x[0], rows_y[-1]))
    page.draw_line((cols_x[-1], rows_y[0]), (cols_x[-1], rows_y[-1]))
    # Interior verticals (optionally skipping the top row to merge the header).
    y_start = rows_y[0] if interior_top else rows_y[1]
    for x in cols_x[1:-1]:
        page.draw_line((x, y_start), (x, rows_y[-1]))
    for r, row in enumerate(cells):
        for c, val in enumerate(row):
            if val is not None:
                page.insert_text((cols_x[c] + 5, rows_y[r] + 14), str(val), fontsize=9)
    if caption:
        page.insert_text((cols_x[0], rows_y[-1] + 22), caption, fontsize=9)
    doc.save(path)
    doc.close()


@pytest.fixture
def simple_table_pdf(tmp_path):
    path = tmp_path / "t.pdf"
    _grid_pdf(
        path,
        cols_x=[72, 172, 272, 372],
        rows_y=[100, 120, 140, 160],
        cells=[["Name", "Mass", "Cost"], ["Rod", "12.5", "3"], ["Disk", "7.0", "42"]],
        caption="Table 1: Sample properties.",
    )
    return str(path)


@pytest.fixture
def merged_header_pdf(tmp_path):
    path = tmp_path / "m.pdf"
    _grid_pdf(
        path,
        cols_x=[72, 172, 272, 372],
        rows_y=[100, 120, 140, 160],
        cells=[["Properties", None, None], ["Rod", "12.5", "3"], ["Disk", "7.0", "42"]],
        interior_top=False,
    )
    return str(path)


def _convert(path, **kw):
    pages = extract_document(path)
    pages, _ = filter_zones(pages)
    meta = DocMeta(source="t.pdf", date=datetime.date.today().isoformat(),
                   tool_version=__version__)
    doc = build_document(pages, meta, detect_title_flag=False)
    return doc, render(doc, no_title=True, **kw)


def test_grid_dimensions_and_alignment(simple_table_pdf):
    page = extract_document(simple_table_pdf)[0]
    tds = detect_tables(page)
    assert len(tds) == 1
    grid = tds[0]["table"].grid
    assert (grid.nrows, grid.ncols) == (3, 3)
    # Column 0 is text (left); numeric columns 1, 2 are right-aligned.
    assert grid.col_aligns == ["l", "r", "r"]


def test_caption_label_stripped(simple_table_pdf):
    page = extract_document(simple_table_pdf)[0]
    cap = detect_tables(page)[0]["table"].caption
    assert cap == "Sample properties."  # "Table 1:" removed


def test_cell_contents(simple_table_pdf):
    page = extract_document(simple_table_pdf)[0]
    cells = {(c.r0, c.c0): c.lines for c in detect_tables(page)[0]["table"].grid.cells}
    assert cells[(0, 0)] == ["Name"]
    assert cells[(1, 0)] == ["Rod"]
    assert cells[(2, 2)] == ["42"]


def test_render_ruled(simple_table_pdf):
    _, tex = _convert(simple_table_pdf)
    assert "\\begin{tabular}{|l|r|r|}" in tex
    assert "Name & Mass & Cost \\\\" in tex
    assert tex.count("\\hline") == 4
    assert "\\caption{Sample properties.}" in tex


def test_render_booktabs(simple_table_pdf):
    _, tex = _convert(simple_table_pdf, booktabs=True)
    assert "\\begin{tabular}{lrr}" in tex
    assert "\\toprule" in tex and "\\midrule" in tex and "\\bottomrule" in tex
    assert "\\hline" not in tex


def test_merged_header_multicolumn(merged_header_pdf):
    page = extract_document(merged_header_pdf)[0]
    grid = detect_tables(page)[0]["table"].grid
    header = next(c for c in grid.cells if c.r0 == 0 and c.c0 == 0)
    assert header.colspan == 3
    assert header.lines == ["Properties"]
    out = _render_table(detect_tables(page)[0]["table"], booktabs=False)
    assert "\\multicolumn{3}{|c|}{Properties}" in out


def test_table_consumes_lines_no_stray_paragraphs(simple_table_pdf):
    doc, _ = _convert(simple_table_pdf)
    # Exactly one Table element and no paragraph echoing the cell text.
    assert sum(isinstance(e, Table) for e in doc.elements) == 1
    from pdf2tex.models import Paragraph
    paras = " ".join(
        i.text for e in doc.elements if isinstance(e, Paragraph) for i in e.inlines
    )
    assert "Rod" not in paras and "Disk" not in paras


@pytest.mark.parametrize("booktabs", [False, True])
def test_table_compiles(simple_table_pdf, tmp_path, booktabs):
    engine = shutil.which("pdflatex")
    if engine is None:
        pytest.skip("no LaTeX engine installed")
    _, tex = _convert(simple_table_pdf, booktabs=booktabs)
    out = tmp_path / "out.tex"
    out.write_text(tex, encoding="utf-8")
    proc = subprocess.run(
        [engine, "-interaction=nonstopmode", "-halt-on-error", out.name],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout[-2000:]
