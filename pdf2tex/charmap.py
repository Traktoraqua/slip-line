"""Text-mode character handling for LaTeX output.

Phase 1 ships the essential escaping needed for a compilable document:
LaTeX-reserved characters in text mode.  Quotes, dashes and richer
substitutions arrive in the dedicated charmap phase; Latin/Norwegian letters
(æ, ø, å, accents) pass through as UTF-8, which modern ``pdflatex`` handles.
"""

from __future__ import annotations

# Order matters: backslash and braces are handled via dedicated escapes so we
# never double-escape the backslashes we introduce.
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


def escape_text(text: str) -> str:
    """Escape LaTeX-reserved characters for use in text mode."""
    return "".join(_RESERVED.get(ch, ch) for ch in text)
