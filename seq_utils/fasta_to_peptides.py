from Bio import SeqIO
from pyteomics import parser

def tryptic_digest(sequence):
    # Define trypsin digestion rules using pyteomics parser
    # The rule 'K' or 'R', except when followed by 'P'
    peptides = list(parser.cleave(sequence, parser.expasy_rules['trypsin']))
    
    return peptides

# Function to read FASTA file and create peptides
def create_tryptic_peptides(fasta_file):
    peptides = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_sequence = str(record.seq)  # Get the protein sequence as string
        protein_peptides = tryptic_digest(protein_sequence)  # Get the tryptic peptides
        peptides.extend(protein_peptides)  # Add them to the overall list
    return peptides

if __name__ == "__main__":
    # Specify the path to your FASTA file
    fasta_file = "fasta/UP000005640_9606.fasta"

    # Generate tryptic peptides
    peptides = create_tryptic_peptides(fasta_file)

    # Display the first 10 peptides (if you have a large number of peptides)
    print("First 10 Tryptic Peptides:")
    for peptide in peptides[:10]:
        print(peptide)

