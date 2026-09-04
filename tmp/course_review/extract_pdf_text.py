#!/usr/bin/env python3
"""Extract page-delimited text and basic metadata from course PDFs."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from pypdf import PdfReader


def safe_name(path: Path) -> str:
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", path.stem).strip("_")
    return stem or "document"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("pdfs", nargs="+", type=Path)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    manifest = []

    for pdf in args.pdfs:
        record = {"source": str(pdf)}
        try:
            reader = PdfReader(str(pdf))
            output = args.out / f"{safe_name(pdf)}.txt"
            chunks = []
            chars = 0
            nonempty_pages = 0
            for page_number, page in enumerate(reader.pages, start=1):
                text = page.extract_text() or ""
                text = text.replace("\x00", "").strip()
                if text:
                    nonempty_pages += 1
                chars += len(text)
                chunks.append(f"\n===== PAGE {page_number} =====\n{text}\n")
            output.write_text("".join(chunks), encoding="utf-8")
            record.update(
                {
                    "text": str(output),
                    "pages": len(reader.pages),
                    "nonempty_pages": nonempty_pages,
                    "characters": chars,
                    "title": (reader.metadata.title if reader.metadata else None),
                    "author": (reader.metadata.author if reader.metadata else None),
                    "status": "ok",
                }
            )
        except Exception as exc:
            record.update({"status": "error", "error": repr(exc)})
        manifest.append(record)

    (args.out / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
