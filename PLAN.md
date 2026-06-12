# PLAN.md — pdf2tex: Rule-Based PDF → LaTeX Converter (CLI)

## 1. Goal and Scope

A Python CLI tool that converts **digitally generated, single-column technical
reports/standards** (real text layer, e.g. DNV-style documents) into a
**complete, compilable `.tex` document**.

**Converts:** body text, headings/sections, lists (bulleted, numbered, nested),
ruled tables, equations (inline and display, when present as selectable text),
and Latin/special characters (æ, ø, å, accents, LaTeX-reserved characters).

**Ignores:** figures (image blocks), page headers, page footers, page numbers.

**Equation fallback:** equations rendered as images cannot be recovered
rule-based → emit `\todo{Equation (image), source p. N}` with the source page
number, using the `todonotes` package.

**Engine:** pure rule-based, built on **PyMuPDF (`fitz`)**. No OCR, no LLM
calls, fully offline and deterministic.

-----

## 2. Architecture Overview

Pipeline of independent stages, each consuming/producing typed dataclasses
(`models.py`). Intermediate state can be dumped to JSON for debugging.

```
PDF ─▶ extract ─▶ zones ─▶ structure ─▶ [lists | tables | math] ─▶ charmap ─▶ emit ─▶ .tex
        (raw)    (filter)  (classify)        (specialize)          (escape)   (render)
```

### Module layout

```
pdf2tex/
├── cli.py            # argparse entry point
├── models.py         # dataclasses: Span, Line, Block, Element tree
├── extract.py        # PyMuPDF: spans, fonts, bboxes, drawings, image blocks
├── zones.py          # header/footer/page-number/figure removal
├── structure.py      # body-font analysis, headings, paragraph assembly
├── lists.py          # bullet/numbered list detection, nesting
├── tables.py         # ruled-table detection from line drawings → tabular
├── mathconv/
│   ├── detect.py     # math region detection (fonts, symbols, layout)
│   ├── unicode_map.py# Unicode math → LaTeX macro table
│   └── ocr.py        # OPTIONAL: math-OCR fallback for image equations
├── charmap.py        # text escaping, ligatures, æøå, reserved chars
└── emit.py           # preamble + element tree → .tex string
```

-----

## 3. Stage Specifications

### 3.1 `extract.py` — Raw extraction

- `page.get_text("dict")` → blocks → lines → spans with: text, font name,
  size, flags (bold/italic/superscript), color, bbox, origin (baseline).
- `page.get_drawings()` → vector line segments (for table rule detection).
- Image blocks recorded with bbox only (needed by `zones` and for the
  equation-image heuristic), never extracted.
- Normalize: Unicode NFC, ligature expansion (ﬁ→fi, ﬂ→fl), strip zero-width
  characters.
- Output: `Page` objects with `spans`, `rules` (h/v segments), `images`.

### 3.2 `zones.py` — Header/footer/figure filtering

- **Headers/footers:** cluster text lines whose y-position falls in the top or
  bottom ~8 % of the page AND whose (normalized) text repeats on ≥ 60 % of
  pages (allowing page-number variation via regex `\d+` masking). Drop them.
- **Page numbers:** isolated short numeric lines near top/bottom margin → drop.
- **Figures:** drop all image blocks. Drop caption lines matching
  `^(Figure|Fig\.?|Figur)\s+[\dA-Z][-.\d]*` that sit adjacent (within ~2 line
  heights) to a dropped image bbox. Flag `--keep-captions` disables caption
  removal.
- Exception: small inline images on a text line are *not* figures — they are
  candidate **equation images** (see 3.6).

### 3.3 `structure.py` — Document structure

- **Body font:** mode of (font, size) weighted by character count → body style.
- **Headings:** a line is a heading if (a) numbering regex matches
  `^\d+(\.\d+)*\s+\S` and (b) font is larger and/or bold relative to body.
  Map numbering depth → `\section`, `\subsection`, `\subsubsection`,
  `\paragraph` (depth ≥ 4). Unnumbered bold standalone lines of larger size →
  `\section*`-style (configurable).
- **Paragraphs:** merge consecutive lines with same style and normal line
  spacing; paragraph break on vertical gap > 1.5 × median leading or indent
  change. **Dehyphenation:** join `word-\n` + lowercase continuation.
- Inline style: bold spans → `\textbf{}`, italic → `\textit{}` (skip italics
  inside math regions).

### 3.4 `lists.py`

- **Bullets:** line starts with `•`, `◦`, `–`, `-`, `*` followed by space.
- **Numbered:** `^\(?(\d+|[a-z]|[ivx]+)[.)]\s`.
- **Nesting:** cluster start-x indents into levels (tolerance ~3 pt) → nested
  `itemize`/`enumerate` (max depth 4, LaTeX limit).
- Continuation lines (indented, no marker) append to the current `\item`.

### 3.5 `tables.py` — Ruled tables

- From `rules`: snap nearly-horizontal/vertical segments to a grid; a table
  region = connected set of intersecting h/v rules forming ≥ 2 × 2 cells.
- Build cell grid from rule intersections; assign spans to cells by bbox
  containment; detect merged cells (cell spanning multiple grid columns/rows)
  → `\multicolumn` (and note: `\multirow` via `multirow` package).
- Column alignment: left by default; numeric-majority columns → right.
- Emit `tabular` with `\hline`/`|` reproducing the ruled look (flag
  `--booktabs` switches to `\toprule/\midrule/\bottomrule`, no vertical rules).
- Multi-line cell text → `\makecell` or paragraph columns `p{width}` when cell
  text wraps; widths estimated from grid geometry as fractions of
  `\textwidth`.
- Table caption lines (`^Table\s+\d`) adjacent to the region → `\caption{}`
  inside a `table` float.

### 3.6 `mathconv/` — Equations

**Detection (`detect.py`)**, a span/line is math-flagged if any of:

- font name matches math families (`CMMI`, `CMSY`, `CMEX`, `Symbol`,
  `*Math*`, `MTMI`, …),
- it contains Unicode math symbols (Greek, ∑ ∫ √ ≤ ≥ ≈ ∂ ∞ ×, sub/superscript
  code points),
- vertical baseline offset + smaller size relative to neighbors
  (sub/superscripts, also flagged by PyMuPDF span flags).

**Display vs inline:**

- Display: an isolated line, horizontally centered (or consistently indented),
  majority math-flagged, optionally with right-aligned `(\d+[.\d]*)` equation
  number → `\begin{equation}` with `\label{eq:p<page>-<n>}`; without number →
  `\[ … \]`.
- Inline: math-flagged span run inside a normal paragraph → `$ … $`.

**Reconstruction:**

- `unicode_map.py`: table of Unicode → LaTeX (`α→\alpha`, `≤→\le`,
  `∑→\sum`, `√→\sqrt`, `°→^\circ`, …).
- Sub/superscript: baseline offset + size ratio → `_{…}` / `^{…}` (group
  consecutive shifted spans).
- Simple fraction heuristic: short centered span stack separated by a
  horizontal rule segment → `\frac{a}{b}` (best effort; bail out to flat text
  if geometry is ambiguous).
- **Image equations:** inline image whose bbox overlaps a text line, or a
  display-positioned isolated image → default `\todo[inline]{Equation (image), source p. N}`; with `--math-ocr` enabled, routed through `ocr.py` (see
  3.6b) instead.
- Any math region the reconstructor cannot parse confidently → emit the raw
  symbols converted via the Unicode map inside `$…$` and append
  `% CHECK: low-confidence math, p. N` comment.

### 3.6b `mathconv/ocr.py` — Optional math-OCR fallback (image equations)

Off by default; the core pipeline stays rule-based and deterministic.
Enabled with `--math-ocr {pix2tex,surya,pix2text}`.

- **Crop & render:** for each detected image-equation bbox, render at high
  resolution via `page.get_pixmap(clip=bbox, dpi=300)` with a small padding
  margin (~4 pt) → PNG in a temp dir.
- **Backend interface:** `MathOcrBackend.predict(png_bytes) -> OcrResult (latex, confidence|None)`. Three implementations behind lazy imports so the
  heavy ML dependencies are only loaded when the flag is used:
  - `pix2tex` (LaTeX-OCR) — block equations, lightweight, pip-installable.
  - `surya` (`surya_latex_ocr`) — handles inline math mixed with text;
    successor to Texify.
  - `pix2text` (MFR model) — strong accuracy on typeset formulas.
- **Sanity validation:** reject OCR output that is empty, exceeds a length
  ratio threshold vs. bbox area, or has unbalanced braces/`\left`/`\right`;
  attempt a quick parse with `pylatexenc` as a syntax gate.
- **Emission:** accepted output → `\begin{equation} … \end{equation}` (or
  `$…$` for inline-sized bboxes) annotated with `% OCR (<backend>): verify, p. N`; rejected output → fall back to the standard `\todo` placeholder.
- **Reporting:** OCR’d equation count and rejection count included in the
  end-of-run warning summary; `--debug` saves crop PNG + OCR string pairs
  for manual review.
- Backends are **optional extras** in packaging:
  `pip install pdf2tex[pix2tex]`, `[surya]`, `[pix2text]`.

### 3.7 `charmap.py` — Character handling

- Escape LaTeX-reserved: `# $ % & _ { } ~ ^ \` (text mode only — never inside
  emitted math).
- Quotes: `“ ” ‘ ’` → ``` '' ` '`; dashes `– —` → `-- ---`;
  non-breaking space → `~`; `… → \dots`.
- Norwegian/Latin letters (æøå, é, ü, …) pass through as UTF-8 (modern
  `pdflatex` handles UTF-8 natively; preamble still loads `[T1]{fontenc}`).

### 3.8 `emit.py` — Output

- Default preamble: `\documentclass[11pt,a4paper]{article}`, packages:
  `fontenc(T1)`, `babel` (option via `--lang`, default `english`), `geometry`,
  `amsmath, amssymb`, `booktabs`, `multirow`, `makecell`, `enumitem`,
  `todonotes`, `graphicx` (harmless, for future use).
- `--preamble FILE` substitutes a user template (body injected at
  `%%BODY%%`).
- Title heuristic: largest text on page 1 → `\title{}` + `\maketitle`
  (flag `--no-title` to disable).
- Wrap output at ~90 columns; blank line between paragraphs; comment header
  with source filename, date, tool version.

### 3.9 `cli.py`

```
pdf2tex input.pdf [-o output.tex]
  --pages 1-12,15       page selection
  --preamble FILE       custom preamble template
  --lang english|norsk  babel language
  --booktabs            booktabs-style tables
  --keep-captions       keep figure captions as plain text
  --no-title            skip title detection
  --math-ocr BACKEND    OCR image equations: pix2tex | surya | pix2text
                        (default: off → \todo placeholders)
  --debug DIR           dump per-stage JSON + annotated page PNGs
  -v / -q               logging verbosity
```

Exit code 0 with warnings summary (count of `\todo` placeholders and
low-confidence math regions) printed at the end.

-----

## 4. Data Model (models.py)

```python
Span(text, font, size, flags, bbox, origin, mathflag)
Line(spans, bbox)
Element  # union via subclasses:
  Heading(level, number, text)
  Paragraph(inlines)            # inline runs: Text/Bold/Italic/InlineMath
  ListBlock(kind, items, level)
  Table(grid, caption)
  DisplayMath(latex, number, label, confidence)
  TodoPlaceholder(reason, page)
Document(meta, elements)
```

-----

## 5. Implementation Phases

Each phase ends with a compilable `.tex` from the test corpus.

1. **Skeleton:** models, extract, naive paragraphs, emit, CLI. Plain text in,
   compilable doc out.
1. **Zones:** header/footer/page-number/figure removal (needs multi-page
   repeat analysis).
1. **Structure:** body-font analysis, numbered headings, dehyphenation,
   bold/italic inline styles.
1. **Charmap:** escaping, quotes/dashes, æøå round-trip test.
1. **Lists:** bullets, numbering, nesting.
1. **Tables:** rule snapping, grid, multicolumn, captions, `--booktabs`.
1. **Math (text):** detection, Unicode map, sub/superscripts, display vs
   inline, equation numbers.
1. **Math (images) + polish:** `\todo` placeholders, confidence comments,
   `--debug` visual dumps (page PNGs with colored bboxes per element class),
   warning summary.
1. **Math-OCR fallback (optional):** `ocr.py` backend interface, bbox
   crop/render, pix2tex backend first (lightest install), then surya and
   pix2text; sanity validation, `% OCR: verify` annotations, optional-extras
   packaging. Evaluate backends on the synthetic round-trip corpus (known
   LaTeX → PDF with rasterized equations → OCR → compare).

-----

## 6. Testing

- **Corpus:** 3–5 representative report/standard PDFs + tiny synthetic PDFs
  generated from known LaTeX (compile LaTeX → PDF → convert back → compare),
  giving ground truth for tables, lists, math.
- **Golden files:** per-stage JSON snapshots and final `.tex` diffs (pytest).
- **Compile check:** run `latexmk -pdf` (or `tectonic`) on every output in CI;
  failure = test failure.
- **Unit tests:** charmap escaping, unicode_map completeness, list-nesting
  edge cases, rule-snapping tolerance.

-----

## 7. Known Risks / Mitigations

|Risk                                                        |Mitigation                                                                                                             |
|------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------|
|Sub/superscript misclassification (footnote markers vs math)|require math context or math font; footnote pattern `^\d$` at line end → `\footnotemark` comment                       |
|Reading order anomalies (text boxes, margin notes)          |sort blocks by (y, x); flag overlapping blocks in `--debug`                                                            |
|Borderless cells inside ruled tables (partial rules)        |grid inferred from outer frame + column x-clusters of span positions                                                   |
|Math fonts embedded with subset names (`ABCDEF+CMMI10`)     |strip subset prefix before font-family matching                                                                        |
|False heading hits (bold lead-ins)                          |require numbering OR (size delta AND standalone line)                                                                  |
|Fraction/radical geometry too complex                       |bail to flat Unicode-mapped math + `% CHECK` comment — never silently wrong                                            |
|Math-OCR hallucination (neural backends)                    |off by default; syntax gate + brace balancing; every OCR result annotated `% OCR: verify`; rejects fall back to `\todo`|

-----

## 8. Dependencies

- Python ≥ 3.10, `pymupdf`, `pytest` (dev), `latexmk`/`tectonic` (CI compile
  check only). Core: no network, no OCR, no API calls — fully deterministic.
- Optional extras (only with `--math-ocr`): `pix2tex` / `surya-ocr` /
  `pix2text`, plus `pylatexenc` for OCR-output syntax validation. Installed
  via `pip install pdf2tex[pix2tex]` etc.; lazy-imported so the core tool
  runs without them.