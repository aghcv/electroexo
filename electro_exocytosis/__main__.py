from __future__ import annotations

import sys

from electro_exocytosis.cli import app


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        sys.argv.pop(1)
    app()
