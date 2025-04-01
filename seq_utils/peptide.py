import re


def remove_non_il(peptides: list) -> list:
    """
    Function to remove non-I/L characters from a peptide sequence.
    """
    return [p for p in peptides if "I" in p or "L" in p]


def remove_ux_containing(peptides: list) -> list:
    """
    Function to remove non-I/L characters from a peptide sequence.
    """
    return [p for p in peptides if "U" not in p and "X" not in p]


def switch_first_il(peptide):
    return re.sub(r"[IL]", lambda x: "L" if x.group() == "I" else "I", peptide, count=1)
