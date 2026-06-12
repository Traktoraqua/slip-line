"""Text-mode character handling for LaTeX output.

Two concerns, applied in a single pass:

* **Reserved characters** — ``# $ % & _ { } ~ ^ \\`` are escaped so they
  render literally (text mode only — never call this on emitted math).
* **Typographic substitution** — curly quotes, en/em dashes, the ellipsis and
  the various non-breaking / exotic spaces map to their LaTeX equivalents; the
  soft hyphen is dropped.

Latin/Norwegian letters (æ, ø, å, é, ü, …) pass through unchanged as UTF-8;
modern ``pdflatex`` handles them natively (the preamble loads ``[T1]{fontenc}``
and ``[utf8]{inputenc}``).
"""

from __future__ import annotations

# Reserved characters. Backslash and braces get dedicated escapes so the
# backslashes we introduce are never themselves re-escaped (single pass).
_RESERVED = {
    "\\": r"\textbackslash{}",
    "{": r"\{",
    "}": r"\}",
    "#": r"\#",
    "$": r"\$",
    "%": r"\%",
    "&": r"\&",
    "_": r"\_",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}

# Typographic punctuation and spacing (keys as \u escapes to stay unambiguous).
_TYPOGRAPHIC = {
    # Curly quotes.
    "‘": "`",        # ' left single
    "’": "'",        # ' right single / apostrophe
    "“": "``",       # " left double
    "”": "''",       # " right double
    "‚": "`",        # ‚ single low-9
    "„": "``",       # „ double low-9
    # Dashes / hyphen-likes.
    "–": "--",       # – en dash
    "—": "---",      # — em dash
    "−": "-",        # − minus sign (text mode)
    "‑": "-",        # ‑ non-breaking hyphen
    # Ellipsis.
    "…": r"\dots{}",  # …
    # Non-breaking spaces → tie.
    " ": "~",        # no-break space
    " ": "~",        # narrow no-break space
    " ": "~",        # figure space
    # Other Unicode spaces → ordinary space.
    " ": " ", " ": " ", " ": " ", " ": " ",
    " ": " ", " ": " ", " ": " ", " ": " ",
    "　": " ",        # ideographic space
    # Soft hyphen → removed.
    "­": "",
}

# Single combined table; the two key sets are disjoint.
_MAP = {**_TYPOGRAPHIC, **_RESERVED}


def escape_text(text: str) -> str:
    """Escape reserved characters and apply typographic substitutions."""
    return "".join(_MAP.get(ch, ch) for ch in text)
