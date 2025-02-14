import sys
from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner
import matplotlib.pyplot as plt

def read_fasta(file_path):
    """Read a sequence from a FASTA file."""
    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            return str(record.seq)
    raise ValueError("No sequence found in the FASTA file.")

def align_sequences(seq1, seq2):
    """Align two sequences using global alignment."""
    aligner = PairwiseAligner()
    aligner.mode = 'global'  # Use global alignment
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")  # Use BLOSUM62 matrix
    alignments = aligner.align(seq1, seq2)
    return alignments[0]  # Return the best alignment

def trim_unpaired_ends(alignment):
    """Trim unpaired sequences from the ends of the alignment."""
    aligned_seq1 = alignment.aligned[0]
    aligned_seq2 = alignment.aligned[1]

    # Find the start and end of the aligned region
    start = max(aligned_seq1[0][0], aligned_seq2[0][0])
    end = min(aligned_seq1[-1][1], aligned_seq2[-1][1])

    # Trim the sequences
    trimmed_seq1 = alignment.target[start:end]
    trimmed_seq2 = alignment.query[start:end]

    return trimmed_seq1, trimmed_seq2

def plot_alignment(seq1, seq2, alignment):
    """Visualize the alignment using matplotlib."""
    aligned_seq1 = alignment.target
    aligned_seq2 = alignment.query

    fig, ax = plt.subplots(figsize=(10, 2))
    ax.text(0.5, 0.8, f"Sequence 1: {aligned_seq1}", fontsize=12, ha='center')
    ax.text(0.5, 0.5, f"Sequence 2: {aligned_seq2}", fontsize=12, ha='center')
    ax.text(0.5, 0.2, f"Score: {alignment.score}", fontsize=12, ha='center')
    ax.axis('off')
    plt.show()

def main():
    if len(sys.argv) != 3:
        print("Usage: python sequence_alignment.py <fasta_file1> <fasta_file2>")
        sys.exit(1)

    # Read sequences from FASTA files
    seq1 = read_fasta(sys.argv[1])
    seq2 = read_fasta(sys.argv[2])

    # Align sequences
    alignment = align_sequences(seq1, seq2)

    # Trim unpaired ends
    trimmed_seq1, trimmed_seq2 = trim_unpaired_ends(alignment)

    # Display results
    print("Original Alignment:")
    print(alignment)
    print("\nTrimmed Sequences:")
    print(f"Sequence 1: {trimmed_seq1}")
    print(f"Sequence 2: {trimmed_seq2}")

    # Plot alignment
    plot_alignment(trimmed_seq1, trimmed_seq2, alignment)

if __name__ == "__main__":
    main()
