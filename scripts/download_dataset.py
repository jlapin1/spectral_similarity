#!/usr/bin/env python3
"""Download a processed dataset from ProteomeXchange / PRIDE.

This supports extending the *observed* branch to a new instrument. It lists and
downloads files for a given accession using the ``ppx`` client.

Important constraints in this project's environment:

* Vendor raw files (Thermo ``.raw``, Bruker ``.d``) are large (tens to hundreds
  of GB per project) and still need conversion (``ThermoRawFileParser`` /
  ``msconvert``) and a database search (MaxQuant / FragPipe / Sage) before they
  fit the pipeline. Those tools are not bundled here.
* Prefer projects that publish already-processed peak lists (``.mgf`` / ``.mzML``)
  and search results, so they can be conformed to the schemas in ``CLAUDE.md``
  section 6 without re-searching.

Requires: ``pip install ppx`` and network access.

Examples
--------
    # list the files available for an accession
    python scripts/download_dataset.py PXD004732 --list

    # download only MGF and search-result files into temp_data/PXD004732/
    python scripts/download_dataset.py PXD004732 \
        --pattern ".mgf" --pattern "msms.txt" --pattern ".mzML" \
        --dest temp_data/PXD004732
"""

import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accession", help="ProteomeXchange/PRIDE accession, e.g. PXD004732.")
    parser.add_argument("--list", action="store_true", help="List remote files and exit (no download).")
    parser.add_argument(
        "--pattern",
        action="append",
        default=[],
        help="Substring filter for filenames; repeatable. If omitted, all files match.",
    )
    parser.add_argument(
        "--dest",
        default=None,
        help="Download directory (default: temp_data/<accession>).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        import ppx
    except ImportError:
        sys.exit("ppx is not installed. Run: pip install ppx")

    dest = args.dest or os.path.join(REPO_ROOT, "temp_data", args.accession)
    os.makedirs(dest, exist_ok=True)

    project = ppx.find_project(args.accession, local=dest)
    remote_files = project.remote_files()

    patterns = args.pattern
    if patterns:
        selected = [f for f in remote_files if any(p.lower() in f.lower() for p in patterns)]
    else:
        selected = remote_files

    print(f"{args.accession}: {len(remote_files)} remote files, {len(selected)} match.")
    for name in selected:
        print(f"  {name}")

    if args.list:
        return

    if not selected:
        sys.exit("No files matched the given --pattern filters; nothing to download.")

    print(f"Downloading {len(selected)} files into {dest} ...")
    project.download(selected)
    print("Download complete. Next: conform the files to the schemas in CLAUDE.md section 6.")


if __name__ == "__main__":
    main()
