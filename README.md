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
PDF ─▶ extract ─▶ zones ─▶ structure ─▶ [lists | tables | math] ─▶ emit ─▶ .tex
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
- **`lists.py`** — bullet/numbered/lettered/roman marker detection, indent
  clustering into nesting levels (≤ 4) and continuation-line merging →
  nested `itemize`/`enumerate`.
- **`tables.py`** — ruled-table detection: rule clustering, grid snapping,
  span-to-cell assignment, merged-cell detection (`\multicolumn`/`\multirow`),
  numeric-column right-alignment, and caption attachment; `--booktabs` switches
  to `\toprule`/`\midrule`/`\bottomrule`.
- **`mathconv/`** — text-equation handling: Unicode→LaTeX map (`unicode_map`),
  math-span detection, sub/superscript reconstruction, inline `$…$` runs (with
  base-stealing) and isolated/numbered display equations (`equation`+`\label`
  or `\[ \]`), with `% CHECK` comments on low-confidence reconstructions
  (`detect`).  Image equations become `\todo[inline]{Equation (image),
  source p. N}` placeholders.
- **`debug.py`** — `--debug DIR` writes per-stage geometry JSON
  (`extract.json`, `zones.json`), the element tree (`document.json`) and
  annotated `page-N.png` images with colored boxes per element class.
- **`charmap.py`** — escaping of LaTeX-reserved characters plus typographic
  substitution (curly quotes, en/em dashes, ellipsis, non-breaking spaces,
  soft-hyphen removal); UTF-8 Latin/Nordic letters (æøå) pass through.
- **`emit.py`** — full preamble, title/`\maketitle`, paragraph rendering,
  line wrapping, custom-preamble templating.
- **`cli.py`** — argparse entry point (`pdf2tex input.pdf -o out.tex`).

A run prints a warning summary (count of `\todo` placeholders and
low-confidence math regions). The only remaining work is the optional math-OCR
fallback (Phase 9) — see `PLAN.md §5`.

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
