#!/usr/bin/env python3
"""Download and convert timsTOF and Orbitrap Astral benchmark data for the I/L study.

Fetches the DDA arms of the multi-species LFQ benchmarks that match the existing
Orbitrap pipeline, then optionally converts the vendor raw files to MGF so they can
feed the observed-spectrum figures (figure2-6) and, for timsTOF, the observed-CCS arm.

Datasets (see CLAUDE.md section 13):
  * astral_dda  : Orbitrap Astral DDA, Generation Beta (PXD070049 / PXD071205).
                  HCD NCE 30, 15 min. Files like
                  LFQ_Astral_DDA_15min_50ng_Condition_[A/B]_REP[1-3].raw
  * timstof_dda : Bruker timsTOF Pro DDA-PASEF, Van Puyvelde Generation Alpha
                  (PXD028735). .d folders; also carries ion mobility for observed CCS.

Requires ppx (pip install ppx) and network for downloads. Conversion is external:
  * .raw -> MGF : ThermoRawFileParser (https://github.com/compomics/ThermoRawFileParser)
  * .d  -> mzML : ProteoWizard msconvert or tdf2mzml (preserve ion mobility for CCS)
These converters are not bundled. The script detects them and, if missing, prints the
exact commands to run rather than failing silently.

Identification results: use ProteoBench's curated outputs rather than re-searching. The
per-module results repositories live under https://github.com/Proteobench . NOTE: the
figure3-6 spectrum-annotation path needs per-PSM scan info (Raw file, MS/MS scan
number), present in MaxQuant evidence/msms and FragPipe psm.tsv but not in DIA-NN
precursor reports; conform whichever table you use to the columns in CLAUDE.md section 6.

Examples
--------
    python scripts/download_benchmark_data.py astral_dda --list
    python scripts/download_benchmark_data.py astral_dda --convert --identifications
    python scripts/download_benchmark_data.py timstof_dda --list --max-files 6
"""

import argparse
import os
import shutil
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Dataset specifications. See CLAUDE.md section 13.
DATASETS = {
    "astral_dda": {
        "description": "Orbitrap Astral DDA, Generation Beta (HCD NCE 30, 15 min).",
        "accessions": ["PXD070049", "PXD071205"],
        # tokens that must appear, plus the raw extension (last item).
        "patterns": ["lfq_astral_dda", "astral", ".raw"],
        "raw_ext": ".raw",
        "proteobench_results": "https://github.com/Proteobench/Results_quant_ion_DDA_Astral",
    },
    "timstof_dda": {
        "description": "Bruker timsTOF Pro DDA-PASEF, Van Puyvelde Generation Alpha (with ion mobility).",
        "accessions": ["PXD028735"],
        "patterns": ["timstof", "tims", ".d"],
        "raw_ext": ".d",
        "proteobench_results": "https://github.com/Proteobench/Results_quant_ion_DIA_diaPASEF",
    },
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("dataset", choices=sorted(DATASETS), help="Which benchmark dataset to fetch.")
    parser.add_argument("--list", action="store_true", help="List matching remote files and exit.")
    parser.add_argument("--dest", default=None, help="Destination dir (default temp_data/<dataset>).")
    parser.add_argument("--convert", action="store_true", help="Convert downloaded raw files to MGF if a converter is available.")
    parser.add_argument("--identifications", action="store_true", help="Also clone the ProteoBench results repo for identifications.")
    parser.add_argument("--max-files", type=int, default=None, help="Cap the number of files downloaded (quick test).")
    return parser.parse_args()


def select_files(remote_files, patterns):
    """Keep files that carry the raw extension and at least one dataset token."""
    tokens = [p.lower() for p in patterns[:-1]]
    ext = patterns[-1].lower()
    selected = []
    for f in remote_files:
        fl = f.lower()
        if ext in fl and (not tokens or any(t in fl for t in tokens)):
            selected.append(f)
    return selected


def download(dataset, dest, list_only, max_files):
    try:
        import ppx
    except ImportError:
        sys.exit("ppx is not installed. Run: pip install ppx")

    spec = DATASETS[dataset]
    for accession in spec["accessions"]:
        acc_dest = os.path.join(dest, accession)
        os.makedirs(acc_dest, exist_ok=True)
        project = ppx.find_project(accession, local=acc_dest)
        remote = project.remote_files()
        selected = select_files(remote, spec["patterns"])
        print(f"{accession}: {len(remote)} remote files, {len(selected)} match {spec['patterns']}.")
        for f in selected:
            print(f"  {f}")
        if list_only or not selected:
            continue
        todo = selected if max_files is None else selected[:max_files]
        print(f"Downloading {len(todo)} files into {acc_dest} ...")
        project.download(todo)


def find_raw_converter():
    for name in ("ThermoRawFileParser", "ThermoRawFileParser.sh", "thermorawfileparser"):
        path = shutil.which(name)
        if path:
            return [path]
    dll = os.environ.get("THERMORAWFILEPARSER_DLL")
    if dll and shutil.which("dotnet"):
        return ["dotnet", dll]
    return None


def convert(dest, raw_ext):
    if raw_ext == ".raw":
        converter = find_raw_converter()
        if not converter:
            print(
                "\nNo ThermoRawFileParser found. To convert .raw -> MGF install it from\n"
                "https://github.com/compomics/ThermoRawFileParser and run, per file:\n"
                "  ThermoRawFileParser -i <file>.raw -o <outdir> -f 0   # 0 = MGF\n"
                "or set THERMORAWFILEPARSER_DLL and put dotnet on PATH."
            )
            return
        raws = [
            os.path.join(root, fn)
            for root, _dirs, files in os.walk(dest)
            for fn in files
            if fn.lower().endswith(".raw")
        ]
        for raw in raws:
            print(f"Converting {raw} -> MGF ...")
            subprocess.run(converter + ["-i", raw, "-o", os.path.dirname(raw), "-f", "0"], check=True)
        print(f"Converted {len(raws)} file(s).")
    elif raw_ext == ".d":
        print(
            "\ntimsTOF .d folders need ProteoWizard msconvert or tdf2mzml (not auto-run).\n"
            "Convert to mzML preserving ion mobility (needed for observed CCS):\n"
            "  msconvert <file>.d --mzML --filter \"peakPicking true 1-\"\n"
            "or tdf2mzml (https://github.com/mafreitas/tdf2mzml). The ambiguity_search\n"
            "module already reads mzML via pyteomics, so MGF is optional for timsTOF."
        )


def clone_identifications(dataset, dest):
    url = DATASETS[dataset]["proteobench_results"]
    if not shutil.which("git"):
        print(f"\ngit not found. Manually clone identifications from: {url}")
        return
    target = os.path.join(dest, "proteobench_results")
    if os.path.exists(target):
        print(f"Identifications already present at {target}")
        return
    print(f"Cloning ProteoBench results: {url}")
    subprocess.run(["git", "clone", "--depth", "1", url, target], check=True)
    print(
        "NOTE: ProteoBench results are ion/precursor-level quant tables. The figure3-6\n"
        "spectrum-annotation path needs per-PSM scan info (Raw file, MS/MS scan number),\n"
        "present in MaxQuant evidence/msms and FragPipe psm.tsv but not DIA-NN reports.\n"
        "Conform the chosen table to the columns in CLAUDE.md section 6 before the figures."
    )


def main():
    args = parse_args()
    dest = args.dest or os.path.join(REPO_ROOT, "temp_data", args.dataset)
    os.makedirs(dest, exist_ok=True)

    print(f"Dataset: {args.dataset}  ({DATASETS[args.dataset]['description']})")
    download(args.dataset, dest, args.list, args.max_files)

    if args.list:
        return
    if args.convert:
        convert(dest, DATASETS[args.dataset]["raw_ext"])
    if args.identifications:
        clone_identifications(args.dataset, dest)

    print("\nDone. Next: conform identifications and converted spectra to the schemas in CLAUDE.md, then run the figures.")


if __name__ == "__main__":
    main()
