"""Optional math-OCR fallback for image equations (Phase 9).

Off by default — the core pipeline stays rule-based and deterministic.  When
``--math-ocr {pix2tex,surya,pix2text}`` is given, each detected equation-image
placeholder is rendered to a high-resolution crop, passed to the chosen
backend, sanity-checked, and (on success) replaced by an annotated
``DisplayMath``; rejected output keeps the ``\\todo`` placeholder.

Backends lazy-import their heavy ML dependencies, so importing this module is
always cheap and the core tool runs without them installed.
"""

from __future__ import annotations

import io
import logging
import os
from dataclasses import dataclass
from typing import Optional, Protocol

import fitz

from ..models import DisplayMath, Document, TodoPlaceholder

log = logging.getLogger(__name__)

RENDER_DPI = 300
CROP_PADDING = 4.0  # pt around the bbox


@dataclass
class OcrResult:
    latex: str
    confidence: Optional[float] = None


class MathOcrBackend(Protocol):
    name: str

    def predict(self, png_bytes: bytes) -> OcrResult:
        ...


# --------------------------------------------------------------------------- #
# Crop / render
# --------------------------------------------------------------------------- #
def render_crop(pdf_path: str, page_number: int, bbox, dpi: int = RENDER_DPI,
                padding: float = CROP_PADDING) -> bytes:
    """Render the padded *bbox* on *page_number* (1-based) to PNG bytes."""
    doc = fitz.open(pdf_path)
    try:
        page = doc[page_number - 1]
        clip = fitz.Rect(
            bbox[0] - padding, bbox[1] - padding,
            bbox[2] + padding, bbox[3] + padding,
        )
        pix = page.get_pixmap(clip=clip, dpi=dpi)
        return pix.tobytes("png")
    finally:
        doc.close()


# --------------------------------------------------------------------------- #
# Sanity validation
# --------------------------------------------------------------------------- #
def _parse_ok(latex: str) -> bool:
    try:
        from pylatexenc.latexwalker import LatexWalker
    except ImportError:
        return True  # syntax gate is optional
    try:
        LatexWalker(latex).get_latex_nodes()
        return True
    except Exception:
        return False


def validate(latex: str, bbox) -> bool:
    """Reject empty, unbalanced or implausibly long OCR output."""
    s = (latex or "").strip()
    if not s:
        return False
    if s.count("{") != s.count("}"):
        return False
    if s.count(r"\left") != s.count(r"\right"):
        return False
    w, h = bbox[2] - bbox[0], bbox[3] - bbox[1]
    if len(s) > 4 * (w + h) + 40:  # length ratio vs. bbox extent
        return False
    return _parse_ok(s)


# --------------------------------------------------------------------------- #
# Backends (lazy imports)
# --------------------------------------------------------------------------- #
class Pix2TexBackend:
    """LaTeX-OCR (pix2tex): lightweight block-equation backend."""

    name = "pix2tex"

    def __init__(self) -> None:
        try:
            from pix2tex.cli import LatexOCR
        except ImportError as exc:  # pragma: no cover - exercised without the lib
            raise ImportError(
                "pix2tex is not installed; install with 'pip install pdf2tex[pix2tex]'"
            ) from exc
        self._model = LatexOCR()

    def predict(self, png_bytes: bytes) -> OcrResult:  # pragma: no cover
        from PIL import Image

        img = Image.open(io.BytesIO(png_bytes))
        return OcrResult(latex=self._model(img))


class SuryaBackend:
    """Surya LaTeX-OCR: handles inline math mixed with text."""

    name = "surya"

    def __init__(self) -> None:
        try:
            from surya.recognition import RecognitionPredictor  # noqa: F401
        except ImportError as exc:  # pragma: no cover
            raise ImportError(
                "surya is not installed; install with 'pip install pdf2tex[surya]'"
            ) from exc
        from surya.recognition import RecognitionPredictor

        self._predictor = RecognitionPredictor()

    def predict(self, png_bytes: bytes) -> OcrResult:  # pragma: no cover
        from PIL import Image

        img = Image.open(io.BytesIO(png_bytes))
        preds = self._predictor([img])
        text = preds[0].text_lines[0].text if preds and preds[0].text_lines else ""
        return OcrResult(latex=text)


class Pix2TextBackend:
    """pix2text MFR model: strong on typeset formulas."""

    name = "pix2text"

    def __init__(self) -> None:
        try:
            from pix2text import Pix2Text  # noqa: F401
        except ImportError as exc:  # pragma: no cover
            raise ImportError(
                "pix2text is not installed; install with 'pip install pdf2tex[pix2text]'"
            ) from exc
        from pix2text import Pix2Text

        self._model = Pix2Text()

    def predict(self, png_bytes: bytes) -> OcrResult:  # pragma: no cover
        from PIL import Image

        img = Image.open(io.BytesIO(png_bytes))
        return OcrResult(latex=self._model.recognize_formula(img))


_BACKENDS = {
    "pix2tex": Pix2TexBackend,
    "surya": SuryaBackend,
    "pix2text": Pix2TextBackend,
}


def get_backend(name: str) -> MathOcrBackend:
    try:
        factory = _BACKENDS[name]
    except KeyError:
        raise ValueError(f"unknown math-OCR backend: {name!r}")
    return factory()


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def _is_equation_image(element) -> bool:
    return (
        isinstance(element, TodoPlaceholder)
        and element.bbox is not None
        and element.reason.startswith("Equation (image)")
    )


def _save_pair(debug_dir: str, page: int, idx: int, png: bytes, latex: str) -> None:
    base = os.path.join(debug_dir, f"ocr-p{page}-{idx}")
    with open(base + ".png", "wb") as fh:
        fh.write(png)
    with open(base + ".tex", "w", encoding="utf-8") as fh:
        fh.write(latex or "")


def apply_ocr(
    document: Document,
    backend: MathOcrBackend,
    pdf_path: str,
    debug_dir: Optional[str] = None,
) -> tuple[int, int]:
    """Replace accepted equation-image todos with annotated DisplayMath.

    Returns ``(accepted, rejected)`` counts.
    """
    accepted = rejected = 0
    idx = 0
    for i, element in enumerate(document.elements):
        if not _is_equation_image(element):
            continue
        idx += 1
        try:
            png = render_crop(pdf_path, element.page, element.bbox)
            result = backend.predict(png)
        except Exception as exc:  # backend / render failure → keep placeholder
            log.warning("OCR failed on p.%s: %s", element.page, exc)
            rejected += 1
            continue
        if debug_dir:
            _save_pair(debug_dir, element.page, idx, png, result.latex)
        if validate(result.latex, element.bbox):
            document.elements[i] = DisplayMath(
                latex=result.latex.strip(),
                comment=f"OCR ({backend.name}): verify, p. {element.page}",
            )
            accepted += 1
        else:
            rejected += 1  # keep the original \todo placeholder
    log.info("OCR: %d accepted, %d rejected", accepted, rejected)
    return accepted, rejected
