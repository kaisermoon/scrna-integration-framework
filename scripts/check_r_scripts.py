#!/usr/bin/env python3
"""Parse repository R entry points and preserve the Rscript exit status."""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="*", type=Path)
    parser.add_argument("--rscript", default=shutil.which("Rscript") or "Rscript")
    args = parser.parse_args()
    files = args.files or sorted((ROOT / "scripts").glob("*.R"))
    if not files:
        parser.error("no R scripts found")

    result = subprocess.run(
        [
            args.rscript,
            "--vanilla",
            "-e",
            "for (f in commandArgs(trailingOnly=TRUE)) parse(file=f)",
            *map(str, files),
        ],
        check=False,
    )
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
