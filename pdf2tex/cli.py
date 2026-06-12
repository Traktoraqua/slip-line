"""Command-line entry point for pdf2tex."""

from __future__ import annotations

import argparse
import datetime
import logging
import os
import sys

from . import __version__
from .emit import render
from .extract import extract_document
from .models import DocMeta, TodoPlaceholder
from .structure import build_document
from .zones import filter_zones

log = logging.getLogger("pdf2tex")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pdf2tex",
        description="Rule-based PDF → LaTeX converter.",
    )
    p.add_argument("input", help="input PDF file")
    p.add_argument("-o", "--output", help="output .tex file (default: alongside input)")
    p.add_argument("--pages", help='page selection, e.g. "1-12,15"')
    p.add_argument("--preamble", help="custom preamble template (body injected at %%BODY%%)")
    p.add_argument(
        "--lang", default="english", help="babel language (default: english)"
    )
    # Flags reserved for later phases — accepted now so the CLI stays stable.
    p.add_argument("--booktabs", action="store_true", help="booktabs-style tables")
    p.add_argument(
        "--keep-captions", action="store_true", help="keep figure captions as plain text"
    )
    p.add_argument("--no-title", action="store_true", help="skip title detection")
    p.add_argument(
        "--math-ocr",
        choices=["pix2tex", "surya", "pix2text"],
        help="OCR image equations (default: off → \\todo placeholders)",
    )
    p.add_argument("--debug", metavar="DIR", help="dump per-stage debug artifacts")
    p.add_argument("--version", action="version", version=f"pdf2tex {__version__}")
    g = p.add_mutually_exclusive_group()
    g.add_argument("-v", "--verbose", action="count", default=0, help="increase verbosity")
    g.add_argument("-q", "--quiet", action="store_true", help="suppress non-error output")
    return p


def _configure_logging(verbose: int, quiet: bool) -> None:
    if quiet:
        level = logging.ERROR
    elif verbose >= 2:
        level = logging.DEBUG
    elif verbose == 1:
        level = logging.INFO
    else:
        level = logging.WARNING
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def _default_output(input_path: str) -> str:
    base, _ = os.path.splitext(input_path)
    return base + ".tex"


def convert(args: argparse.Namespace) -> str:
    pages = extract_document(args.input, pages=args.pages)
    pages, _ = filter_zones(pages, keep_captions=args.keep_captions)
    meta = DocMeta(
        source=os.path.basename(args.input),
        date=datetime.date.today().isoformat(),
        tool_version=__version__,
        lang=args.lang,
    )
    document = build_document(pages, meta, detect_title_flag=not args.no_title)

    preamble_text = None
    if args.preamble:
        with open(args.preamble, "r", encoding="utf-8") as fh:
            preamble_text = fh.read()

    tex = render(document, no_title=args.no_title, preamble=preamble_text)

    todo_count = sum(isinstance(e, TodoPlaceholder) for e in document.elements)
    if todo_count:
        log.warning("%d \\todo placeholder(s) emitted for manual review", todo_count)
    return tex


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    _configure_logging(args.verbose, args.quiet)

    if not os.path.isfile(args.input):
        log.error("input file not found: %s", args.input)
        return 2

    try:
        tex = convert(args)
    except Exception as exc:  # surface a clean error, no traceback by default
        log.error("conversion failed: %s", exc)
        if args.verbose >= 2:
            raise
        return 1

    output = args.output or _default_output(args.input)
    with open(output, "w", encoding="utf-8") as fh:
        fh.write(tex)
    log.info("Wrote %s", output)
    if not args.quiet:
        print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
