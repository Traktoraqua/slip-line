# pdf2tex

A rule-based PDF → LaTeX converter (CLI) for digitally generated,
single-column technical reports and standards with a real text layer. Built on
[PyMuPDF](https://pymupdf.readthedocs.io/) — fully offline and deterministic,
no OCR and no LLM calls in the core pipeline.

See [`PLAN.md`](PLAN.md) for the full design and roadmap.

## Status — Phases 1–2

The end-to-end pipeline is in place and produces a complete, **compilable**
`.tex` document:

```
PDF ─▶ extract ─▶ zones ─▶ structure ─▶ emit ─▶ .tex
```

Implemented so far:

- **`extract.py`** — text spans (font, size, flags, bbox, origin), vector rule
  segments and image-block bounding boxes via PyMuPDF; Unicode NFC
  normalisation, ligature expansion and zero-width stripping.
- **`zones.py`** — header/footer removal via multi-page digit-masked repeat
  analysis, page-number stripping, and figure/caption removal (kept with
  `--keep-captions`); small inline images preserved as equation candidates.
- **`structure.py`** — body-font analysis; numbered headings
  (`\section`…`\paragraph` by depth) and unnumbered bold standalone headings;
  paragraph grouping with gap/indent breaks, dehyphenation, and bold/italic
  inline runs; title heuristic.
- **`charmap.py`** — escaping of LaTeX-reserved characters; UTF-8 Latin/Nordic
  letters pass through.
- **`emit.py`** — full preamble, title/`\maketitle`, paragraph rendering,
  line wrapping, custom-preamble templating.
- **`cli.py`** — argparse entry point (`pdf2tex input.pdf -o out.tex`).

Later phases add lists, tables and math — see `PLAN.md §5`.

## Install

```bash
pip install -e .          # core
pip install -e .[dev]     # + pytest
```

## Usage

```bash
pdf2tex input.pdf -o output.tex
pdf2tex input.pdf --pages 1-12,15 --lang norsk --no-title
```

Run `pdf2tex --help` for the full option list.

## Tests

```bash
pytest
```

The compile-check test builds the emitted `.tex` with `pdflatex` when a LaTeX
engine is installed, and is skipped otherwise.
