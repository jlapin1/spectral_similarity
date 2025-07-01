import re
import random


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


import re
import random


def switch_random_il(peptide):
    """
    Randomly swap an occurrence of I or L, ignoring any I or L inside square brackets.
    If only a single occurrence is found outside brackets, it will be switched.
    """
    # Find all spans corresponding to brackets
    bracket_spans = [m.span() for m in re.finditer(r"\[[^\]]*\]", peptide)]

    def is_inside_brackets(position):
        for start, end in bracket_spans:
            if start <= position < end:
                return True
        return False

    # Find positions of I and L that are NOT inside brackets
    positions = [
        m.start()
        for m in re.finditer(r"[IL]", peptide)
        if not is_inside_brackets(m.start())
    ]

    if not positions:
        # No I or L found outside brackets
        return peptide

    if len(positions) == 1:
        # Only one occurrence; switch it
        pos = positions[0]
    else:
        # Randomly select one occurrence excluding the first
        pos = random.choice(positions[1:])

    # Perform the swap
    swapped_char = "L" if peptide[pos] == "I" else "I"

    # Reconstruct the peptide
    peptide = peptide[:pos] + swapped_char + peptide[pos + 1 :]
    return peptide
