"""Math detection and reconstruction for text equations (Phase 7).

* :func:`is_math_span` flags spans by math font family or Unicode symbol.
* :func:`shift_class` classifies a span as normal / superscript / subscript
  from its baseline offset and size relative to the line.
* :func:`spans_to_latex` reconstructs ``_{…}`` / ``^{…}`` from classified
  spans.
* :func:`build_display_math` recognises isolated, centred / numbered math
  lines and returns a :class:`~pdf2tex.models.DisplayMath`.
"""

from __future__ import annotations

import re

from ..models import DisplayMath, Line, Page, Span
from .unicode_map import MATH_CHARS, convert_math

# Math font families (subset prefix like "ABCDEF+" stripped before matching).
_MATH_FONTS = (
    "cmmi", "cmsy", "cmex", "msam", "msbm", "mtmi", "mtsy", "symbol",
    "mathjax", "stix", "xits", "euclid", "mathematical",
)
# Operators that signal math even when set in a normal font.
_OPERATOR_CHARS = set("=+−<>≤≥≠≈±×÷∑∫√≡∂∇·⋅/^_")
# Equation number at the right margin, e.g. "(3.1)" or "3".
_EQNUM_RE = re.compile(r"^\(?\d+(?:\.\d+)*\)?$")
_SAFE_MATH = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
                 "0123456789+-*/=<>()[].,|^_{} \\")

LOW_CONFIDENCE = 0.6


def _font_base(font: str) -> str:
    name = font.lower()
    if "+" in name:  # subset prefix, e.g. "abcdef+cmmi10"
        name = name.split("+", 1)[1]
    return name


def is_math_span(span: Span) -> bool:
    name = _font_base(span.font)
    if any(fam in name for fam in _MATH_FONTS):
        return True
    return any(ch in MATH_CHARS for ch in span.text)


def is_equation_number(text: str) -> bool:
    return bool(_EQNUM_RE.match(text.strip()))


def has_math_indicator(spans: list[Span]) -> bool:
    text = "".join(s.text for s in spans)
    return any(is_math_span(s) for s in spans) or any(ch in _OPERATOR_CHARS for ch in text)


def line_ref(spans: list[Span]) -> tuple[float, float]:
    """Baseline (origin y) and size of the line's dominant (largest) span."""
    ref = max(spans, key=lambda s: s.size)
    return ref.origin[1], ref.size


def shift_class(span: Span, ref_y: float, ref_size: float) -> str:
    """Return 'normal', 'sup' or 'sub' for *span* relative to the baseline."""
    if span.is_superscript:
        return "sup"
    dy = ref_y - span.origin[1]  # > 0 ⇒ raised above the baseline
    if span.size < 0.85 * ref_size:
        if dy > 0.15 * ref_size:
            return "sup"
        if dy < -0.10 * ref_size:
            return "sub"
    return "normal"


def spans_to_latex(items: list[tuple[Span, str]]) -> str:
    """Reconstruct LaTeX from (span, shift_class) pairs."""
    parts: list[str] = []
    i = 0
    while i < len(items):
        _, cls = items[i]
        if cls in ("sup", "sub"):
            j = i
            chunk = []
            while j < len(items) and items[j][1] == cls:
                chunk.append(items[j][0].text)
                j += 1
            inner = convert_math("".join(chunk)).strip()
            parts.append(("^" if cls == "sup" else "_") + "{" + inner + "}")
            i = j
        else:
            parts.append(convert_math(items[i][0].text))
            i += 1
    return "".join(parts)


def confidence(latex: str) -> float:
    """Fraction of characters that look like valid math (vs. garbage)."""
    if not latex:
        return 0.0
    good = sum(ch in _SAFE_MATH for ch in latex)
    return good / len(latex)


def _looks_like_prose(text: str) -> bool:
    words = [w for w in re.split(r"[^A-Za-z]+", text) if len(w) >= 4]
    return len(words) > 3


def build_display_math(
    lines: list[Line],
    page: Page,
    page_number: int,
    eq_index: int,
) -> DisplayMath | None:
    """Return a DisplayMath if *lines* form an isolated display equation."""
    spans: list[Span] = []
    for ln in lines:
        spans.extend(sorted((s for s in ln.spans if s.text.strip()),
                            key=lambda s: s.bbox[0]))
    if not spans:
        return None

    # Detect a right-aligned equation number on the (first) line.
    eqnum = None
    if _EQNUM_RE.match(spans[-1].text.strip()) and spans[-1].bbox[0] > page.width * 0.7:
        eqnum = spans[-1]
    content = [s for s in spans if s is not eqnum]
    if not content:
        return None

    text = "".join(s.text for s in content)
    has_indicator = (
        any(is_math_span(s) for s in content)
        or any(ch in _OPERATOR_CHARS for ch in text)
    )
    if not has_indicator or _looks_like_prose(text):
        return None

    first = lines[0]
    center = (first.bbox[0] + first.bbox[2]) / 2.0
    centered = abs(center - page.width / 2.0) <= 0.18 * page.width
    if not (centered or eqnum is not None):
        return None

    ref_y, ref_size = line_ref(content)
    items = [(s, shift_class(s, ref_y, ref_size)) for s in content]
    latex = spans_to_latex(items).strip()
    if not latex:
        return None

    number = label = None
    if eqnum is not None:
        number = eqnum.text.strip().strip("()")
        label = f"eq:p{page_number}-{eq_index}"

    return DisplayMath(
        latex=latex,
        number=number,
        label=label,
        confidence=confidence(latex),
        page=page_number,
    )
