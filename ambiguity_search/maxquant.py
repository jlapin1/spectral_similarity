import argparse
from collections import defaultdict
from pathlib import Path
from typing import ClassVar, Dict, List, Set

import pandas as pd
import pyteomics.mzml
from tqdm import tqdm



class MaxQuantAmbiguitySearch:
    """
    Search peptides with I/L substitutions in MaxQuant results.
    """

    AMBIGUOUS_AA: ClassVar[Set[str]] = {"I", "L"}
    """
    Ambiguities of interes
    """

    AMBIGUOUS_AA_REPL: ClassVar[Dict[str, str]] = {
        "I": "L",
        "L": "I"
    }
    """
    Replacements for ambiguous amino acids.
    """

    FIRST_N: ClassVar[int] = 1
    """
    Number of ambiguities to consider for each peptide. 1 means only the first ambiguity is considered.
    """

    def __init__(self, maxquant_folders: List[Path], mzml_folders: Path):
        """
        Initialize the MaxQuantAmbiguitySearch with a list of MaxQuant results.
        
        Parameters
        ----------
        maxquant_folders : List[Path]
            A list of paths to MaxQuant results containing the 'msms.txt' files.
        """
        self.maxquant_folders = maxquant_folders
        self.mzml_folders = mzml_folders
        self.mzmls = {}

    def search(self) -> pd.DataFrame:
        """
        Search for peptides with I/L substitutions
        """
        peptide_index = defaultdict(set)

        for folder in tqdm(self.maxquant_folders, desc="Builidng peptide index"):
            if not folder.is_dir():
                tqdm.write(f"Skipping non-directory: {folder}")
                continue
            msms_file = folder.joinpath("msms.txt")
            if not msms_file.exists():
                tqdm.write(f"Skipping folder without msms.txt: {folder}")
                continue

            msms_df = pd.read_csv(
                msms_file,
                sep="\t",
                usecols=["Sequence", "Raw file", "Score", "Scan number"]
            )
            msms_df.sort_values(by=["Score"], inplace=True, ascending=False)

            for row in msms_df[["Sequence", "Raw file", "Scan number"]].itertuples():
                for aa in row[1]:
                    if aa in self.AMBIGUOUS_AA:
                        peptide_index[row[1]].add(f"{row[2]}:{row[3]}")
                        break

        ambiguities = []

        for seq, raw_files in tqdm(peptide_index.items(), desc="Searching for ambiguities and loading spectra"):
            ambiguity_matches = 0

            for idx, aa in enumerate(seq):
                if ambiguity_matches == self.FIRST_N:
                    break

                if aa not in self.AMBIGUOUS_AA:
                    continue

                ambiguous_seq = list(seq)
                ambiguous_seq[idx] = self.AMBIGUOUS_AA_REPL[aa]
                ambiguous_seq = "".join(ambiguous_seq) # type: ignore

                ambiguous_raw_files = peptide_index.get(ambiguous_seq)

                if ambiguous_raw_files is not None:
                    raw_files = list(raw_files) # type: ignore
                    ambiguous_raw_files = list(ambiguous_raw_files) # type: ignore
                    seq_spectrum = self.get_spectrum(raw_files[0])
                    ambiguous_seq_spectrum = self.get_spectrum(ambiguous_raw_files[0])  # type: ignore
                    ambiguities.append(
                        [
                            seq,
                            ambiguous_seq,
                            ", ".join(raw_files),
                            ", ".join(ambiguous_raw_files),
                            seq_spectrum[0],
                            seq_spectrum[1],
                            ambiguous_seq_spectrum[0],
                            ambiguous_seq_spectrum[1]
                        ]
                    )
                    ambiguity_matches += 1

        ambiguous_df = pd.DataFrame(
            ambiguities,
            columns=[
                "sequence",
                "ambigous_sequence",
                "sequence_raw_files",
                "ambiguous_sequence_raw_files",
                "sequence_mz",
                "sequence_intensity",
                "ambiguous_sequence_mz",
                "ambiguous_sequence_intensity"
            ]
        )
        return ambiguous_df

    def get_spectrum(self, raw_file: str):
        """
        Get the spectrum for a given raw file and scan number.
        
        Parameters
        ----------
        raw_file : str
            The basename of the raw file and the scan number,
            separated by a colon (e.g., "file.raw:123").
        
        Returns
        -------
        Tupel[List[float], List[float]]
            A tuple containing two lists: m/z values and intensity values of the spectrum.s
        """
        raw_file, scan_number = raw_file.split(":")
        mzml_file = raw_file + ".mzML"

        if mzml_file not in self.mzmls:
            mzml_path = self.mzml_folders.joinpath(mzml_file)
            if not mzml_path.exists():
                tqdm.write(f"Skipping missing mzML file: {mzml_path}")
                return ([], [])

            self.mzmls[mzml_file] = pyteomics.mzml.MzML(str(mzml_path), use_index=True)

        spectrum = self.mzmls[mzml_file].get_by_id(
            f"controllerType=0 controllerNumber=1 scan={scan_number}"
        )
        return spectrum["m/z array"], spectrum["intensity array"]

def get_cli():
    """
    Command line interface for MaxQuantAmbiguitySearch
    """

    parser = argparse.ArgumentParser(
        description="Search for peptides with I/L substitutions in MaxQuant results."
    )
    parser.add_argument(
        "outfile",
        type=Path,
        help="Path the outfile. (.tsv == tab-separated values, .parquet == Parquet format, unknown: .tsv)"
    )
    parser.add_argument(
        "mzml_folder",
        type=Path,
        help="Path to the folder containing mzML files."
    )
    parser.add_argument(
        "maxquant_folders",
        nargs="+",
        type=Path,
        help="Paths to MaxQuant result folders containing 'msms.txt' files."
    )


    return parser



def main():
    """
    Main function for testing
    """
    cli = get_cli()
    args = cli.parse_args()

    searcher = MaxQuantAmbiguitySearch(args.maxquant_folders, args.mzml_folder)
    ambiguous_peptides = searcher.search()

    match args.outfile.suffix:
        case ".parquet":
            ambiguous_peptides.to_parquet(args.outfile, index=False)
        case _:
            ambiguous_peptides.to_csv(args.outfile, index=False, sep="\t")
        
if __name__ == "__main__":
    main()
