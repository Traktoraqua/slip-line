"""Typed data model shared across pipeline stages.

The geometry types (:class:`Span`, :class:`Line`, :class:`Block`,
:class:`Page`) describe what was extracted from the PDF.  The element types
(:class:`Heading`, :class:`Paragraph`, …) describe the logical document that
:mod:`pdf2tex.emit` renders to LaTeX.

Only a subset is populated in Phase 1; the remaining element types are
defined here so later phases can build on a stable model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Union

BBox = tuple[float, float, float, float]
Point = tuple[float, float]


# PyMuPDF span ``flags`` bit positions (see ``TextPage`` docs).
FLAG_SUPERSCRIPT = 1 << 0
FLAG_ITALIC = 1 << 1
FLAG_SERIF = 1 << 2
FLAG_MONOSPACE = 1 << 3
FLAG_BOLD = 1 << 4


# --------------------------------------------------------------------------- #
# Low-level geometry
# --------------------------------------------------------------------------- #
@dataclass
class Span:
    """A run of text sharing one font, size and style, with its geometry."""

    text: str
    font: str
    size: float
    flags: int
    bbox: BBox
    origin: Point
    color: int = 0
    mathflag: bool = False

    @property
    def is_superscript(self) -> bool:
        return bool(self.flags & FLAG_SUPERSCRIPT)

    @property
    def is_italic(self) -> bool:
        # Honour both the PyMuPDF flag and an italic/oblique font name.
        name = self.font.lower()
        return bool(self.flags & FLAG_ITALIC) or "italic" in name or "oblique" in name

    @property
    def is_bold(self) -> bool:
        name = self.font.lower()
        return bool(self.flags & FLAG_BOLD) or "bold" in name or "black" in name


@dataclass
class Line:
    spans: list[Span]
    bbox: BBox

    @property
    def text(self) -> str:
        return "".join(span.text for span in self.spans)


@dataclass
class Block:
    """A text block as reported by PyMuPDF (a contiguous group of lines)."""

    lines: list[Line]
    bbox: BBox

    @property
    def text(self) -> str:
        return "\n".join(line.text for line in self.lines)


@dataclass
class Rule:
    """A near-horizontal or near-vertical vector segment (table rules)."""

    x0: float
    y0: float
    x1: float
    y1: float

    @property
    def horizontal(self) -> bool:
        return abs(self.y1 - self.y0) <= abs(self.x1 - self.x0)

    @property
    def vertical(self) -> bool:
        return not self.horizontal


@dataclass
class ImageBlock:
    """An image region — recorded by bbox only, never extracted."""

    bbox: BBox


@dataclass
class Page:
    number: int  # 1-based source page number
    width: float
    height: float
    blocks: list[Block] = field(default_factory=list)
    rules: list[Rule] = field(default_factory=list)
    images: list[ImageBlock] = field(default_factory=list)

    @property
    def lines(self) -> list[Line]:
        return [line for block in self.blocks for line in block.lines]

    @property
    def spans(self) -> list[Span]:
        return [span for line in self.lines for span in line.spans]


# --------------------------------------------------------------------------- #
# Inline runs (paragraph contents)
# --------------------------------------------------------------------------- #
@dataclass
class Text:
    text: str


@dataclass
class Bold:
    text: str


@dataclass
class Italic:
    text: str


@dataclass
class InlineMath:
    latex: str


Inline = Union[Text, Bold, Italic, InlineMath]


# --------------------------------------------------------------------------- #
# Logical elements
# --------------------------------------------------------------------------- #
class Element:
    """Base class for renderable document elements."""


@dataclass
class Heading(Element):
    level: int
    number: Optional[str]
    text: str


@dataclass
class Paragraph(Element):
    inlines: list[Inline]


@dataclass
class ListBlock(Element):
    kind: str  # "itemize" | "enumerate"
    items: list  # list of (list[Inline] | nested ListBlock)
    level: int = 0


@dataclass
class Table(Element):
    grid: list  # rows of cells
    caption: Optional[str] = None


@dataclass
class DisplayMath(Element):
    latex: str
    number: Optional[str] = None
    label: Optional[str] = None
    confidence: Optional[float] = None


@dataclass
class TodoPlaceholder(Element):
    reason: str
    page: int


# --------------------------------------------------------------------------- #
# Document
# --------------------------------------------------------------------------- #
@dataclass
class DocMeta:
    source: str
    date: str
    tool_version: str
    title: Optional[str] = None
    lang: str = "english"


@dataclass
class Document:
    meta: DocMeta
    elements: list[Element] = field(default_factory=list)
