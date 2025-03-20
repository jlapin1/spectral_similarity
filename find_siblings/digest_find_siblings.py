"""
A simple script to digest a FASTA file and find sibling peptides
(i.e. peptides that are similar except for I/L exchange, e.g. NLFLSK and NIFISK
"""

from pyteomics import fasta, parser
import argparse

minlength = 6
maxlength = 60

def digest_fasta_keep_with_leucines(file_path):
    peps_by_length = {}

    for _, sequence in fasta.read(file_path):
        peptides = parser.cleave(sequence, parser.expasy_rules['trypsin'])

        for peptide in peptides:
            peplen = len(peptide)
            if peplen >= minlength and peplen <= maxlength and ("L" in peptide or "I" in peptide):
                if peplen not in peps_by_length.keys():
                    peps_by_length[peplen] = {}

                pepgroup = peptide.replace('I', 'J').replace('L', 'J')
                if pepgroup not in peps_by_length[peplen].keys():
                    peps_by_length[peplen][pepgroup] = set()
                
                peps_by_length[peplen][pepgroup].add(peptide)
    
    return peps_by_length


if __name__ == '__main__':
    arparser = argparse.ArgumentParser(description='Digest FASTA file and find sibling peptides (similar except for I/L exchange).')
    arparser.add_argument('fasta', type=str, help='Path to the input FASTA file')
    arparser.add_argument('out', type=str, help='Path to the output file to save sibling peptides')
    args = arparser.parse_args()

    peps_by_len = digest_fasta_keep_with_leucines(args.fasta)

    count_siblings = 0
    for length, peps in peps_by_len.items():
        for group_pep in list(peps.keys()):
            if len(peps[group_pep]) == 1:
                del peps[group_pep]
            else:
                count_siblings += len(peps[group_pep])
    
    with open(args.out, 'w') as f:
        f.write(f"Total number of sibling peptides: {count_siblings}\n")
        print(f"Total number of sibling peptides: {count_siblings}")

        for length in sorted(peps_by_len.keys()):
            for group_pep in peps_by_len[length]:
                f.write(f"{peps_by_len[length][group_pep]}\n")

        # for length, peps in peps_by_len.items():
        #     for group_pep in list(peps.keys()):
        #         f.write(f"{peps[group_pep]}\n")
