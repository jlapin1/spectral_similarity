from Bio import SeqIO

# Function to simulate tryptic digestion
def tryptic_digest(sequence):
    # Define the cleavage rules: cleave after K or R, unless followed by P
    peptides = []
    current_peptide = []
    for i in range(len(sequence)):
        aa = sequence[i]
        current_peptide.append(aa)
        
        # Check for cleavage sites (K or R, not followed by P)
        if aa in ['K', 'R'] and (i == len(sequence) - 1 or sequence[i + 1] != 'P'):
            peptides.append(''.join(current_peptide))
            current_peptide = []
    
    # Add the last peptide if it exists
    if current_peptide:
        peptides.append(''.join(current_peptide))
    
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

