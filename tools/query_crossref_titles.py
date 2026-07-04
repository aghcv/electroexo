#!/usr/bin/env python3
"""Query Crossref titles for DOI sanity checks."""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request


def main() -> int:
    for title in sys.argv[1:]:
        url = "https://api.crossref.org/works?rows=5&query.title=" + urllib.parse.quote(title)
        data = json.load(urllib.request.urlopen(url, timeout=20))
        print("QUERY", title)
        for item in data["message"]["items"][:5]:
            print(" ", item.get("DOI"), "|", (item.get("title") or [""])[0], "|", (item.get("container-title") or [""])[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
