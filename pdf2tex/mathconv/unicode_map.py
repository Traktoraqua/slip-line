"""Unicode → LaTeX conversion for math content.

:data:`MAP` covers Greek letters, operators and relations; :data:`SUPERSCRIPTS`
and :data:`SUBSCRIPTS` cover the Unicode super/subscript code points.
:func:`convert_math` turns a raw span string into a LaTeX math-mode fragment.
"""

from __future__ import annotations

import re

# Greek lowercase / uppercase, operators, relations, arrows, misc symbols.
MAP: dict[str, str] = {
    # Greek lowercase
    "α": r"\alpha", "β": r"\beta", "γ": r"\gamma", "δ": r"\delta",
    "ε": r"\epsilon", "ζ": r"\zeta", "η": r"\eta", "θ": r"\theta",
    "ι": r"\iota", "κ": r"\kappa", "λ": r"\lambda", "μ": r"\mu",
    "ν": r"\nu", "ξ": r"\xi", "π": r"\pi", "ρ": r"\rho",
    "σ": r"\sigma", "τ": r"\tau", "υ": r"\upsilon", "φ": r"\phi",
    "ϕ": r"\phi", "χ": r"\chi", "ψ": r"\psi", "ω": r"\omega",
    "ϵ": r"\epsilon", "ϑ": r"\vartheta", "ϱ": r"\varrho", "ς": r"\varsigma",
    # Greek uppercase (only the non-Latin-looking ones)
    "Γ": r"\Gamma", "Δ": r"\Delta", "Θ": r"\Theta", "Λ": r"\Lambda",
    "Ξ": r"\Xi", "Π": r"\Pi", "Σ": r"\Sigma", "Φ": r"\Phi",
    "Ψ": r"\Psi", "Ω": r"\Omega", "Υ": r"\Upsilon",
    # Micro sign (distinct code point from Greek mu)
    "µ": r"\mu",
    # Big operators
    "∑": r"\sum", "∏": r"\prod", "∫": r"\int", "∬": r"\iint",
    "∮": r"\oint", "√": r"\sqrt", "∂": r"\partial", "∇": r"\nabla",
    "∞": r"\infty",
    # Relations
    "≤": r"\le", "≥": r"\ge", "≠": r"\neq", "≈": r"\approx",
    "≡": r"\equiv", "≅": r"\cong", "∼": r"\sim", "∝": r"\propto",
    "≪": r"\ll", "≫": r"\gg", "≃": r"\simeq",
    # Binary operators
    "±": r"\pm", "∓": r"\mp", "×": r"\times", "÷": r"\div",
    "⋅": r"\cdot", "·": r"\cdot", "∗": "*", "∘": r"\circ",
    # Set / logic
    "∈": r"\in", "∉": r"\notin", "⊂": r"\subset", "⊆": r"\subseteq",
    "⊃": r"\supset", "⊇": r"\supseteq", "∪": r"\cup", "∩": r"\cap",
    "∅": r"\emptyset", "∀": r"\forall", "∃": r"\exists", "¬": r"\neg",
    "∧": r"\wedge", "∨": r"\vee",
    # Arrows
    "→": r"\rightarrow", "←": r"\leftarrow", "↔": r"\leftrightarrow",
    "⇒": r"\Rightarrow", "⇐": r"\Leftarrow", "⇔": r"\Leftrightarrow",
    "↦": r"\mapsto",
    # Geometry / misc
    "°": r"^{\circ}", "⊥": r"\perp", "∥": r"\parallel", "∠": r"\angle",
    "′": "'", "″": "''", "‰": r"\permil",
    # Minus / fraction slash
    "−": "-", "⁄": "/",
}

# Convenient fraction glyphs.
MAP.update({
    "½": r"\frac{1}{2}", "¼": r"\frac{1}{4}", "¾": r"\frac{3}{4}",
    "⅓": r"\frac{1}{3}", "⅔": r"\frac{2}{3}",
})

SUPERSCRIPTS: dict[str, str] = {
    "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5",
    "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁺": "+", "⁻": "-",
    "⁼": "=", "⁽": "(", "⁾": ")", "ⁿ": "n", "ⁱ": "i",
}

SUBSCRIPTS: dict[str, str] = {
    "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4", "₅": "5",
    "₆": "6", "₇": "7", "₈": "8", "₉": "9", "₊": "+", "₋": "-",
    "₌": "=", "₍": "(", "₎": ")", "ₐ": "a", "ₑ": "e", "ₒ": "o",
    "ₓ": "x", "ₙ": "n", "ᵢ": "i", "ⱼ": "j",
}

# Function names that should become control words.
_FUNCTIONS = [
    "sinh", "cosh", "tanh", "sin", "cos", "tan", "cot", "sec", "csc",
    "arcsin", "arccos", "arctan", "log", "ln", "exp", "lim", "max", "min",
    "det", "gcd", "arg", "dim", "deg", "ker",
]
_FUNCTION_RE = re.compile(
    r"(?<![\\A-Za-z])(" + "|".join(_FUNCTIONS) + r")(?![A-Za-z])"
)

# The set of code points treated as math symbols (for detection).
MATH_CHARS = set(MAP) | set(SUPERSCRIPTS) | set(SUBSCRIPTS)


def convert_math(text: str) -> str:
    """Convert a raw string to a LaTeX math-mode fragment (best effort)."""
    out: list[str] = []
    i = 0
    n = len(text)
    while i < n:
        ch = text[i]
        if ch in SUPERSCRIPTS:
            j = i
            run = []
            while j < n and text[j] in SUPERSCRIPTS:
                run.append(SUPERSCRIPTS[text[j]])
                j += 1
            out.append("^{" + "".join(run) + "}")
            i = j
            continue
        if ch in SUBSCRIPTS:
            j = i
            run = []
            while j < n and text[j] in SUBSCRIPTS:
                run.append(SUBSCRIPTS[text[j]])
                j += 1
            out.append("_{" + "".join(run) + "}")
            i = j
            continue
        if ch in MAP:
            val = MAP[ch]
            out.append(val)
            if val[-1].isalpha():
                out.append(" ")  # avoid gluing a macro to a following letter
            i += 1
            continue
        if ch in "%#&":  # would otherwise be special even in math mode
            out.append("\\" + ch)
            i += 1
            continue
        out.append(ch)
        i += 1
    result = "".join(out)
    result = _FUNCTION_RE.sub(lambda m: "\\" + m.group(1), result)
    return result
