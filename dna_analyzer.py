# DNA Analyzer Script
# Author: Your Name
# Date: 2025-08-15
# Description: Reads DNA sequences from FASTA, performs analysis, and prints results

# ----------------------------
# 1. Function to read FASTA file
# ----------------------------
def read_fasta(file_path):
    sequence = ""
    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):   # Ignore header lines
                sequence += line.strip().upper()  # Remove newline & convert to uppercase
    return sequence

# ----------------------------
# 2. Function to count nucleotides
# ----------------------------
def count_nucleotides(seq):
    counts = {
        "A": seq.count("A"),
        "T": seq.count("T"),
        "G": seq.count("G"),
        "C": seq.count("C")
    }
    return counts

# ----------------------------
# 3. Function to calculate GC content
# ----------------------------
def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return ((g + c) / len(seq)) * 100

# ----------------------------
# 4. Function to generate reverse complement
# ----------------------------
def reverse_complement(seq):
    complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
    return "".join([complement[base] for base in reversed(seq)])

# ----------------------------
# 5. Function to find motif positions
# ----------------------------
def find_motif(seq, motif):
    positions = []
    motif = motif.upper()
    for i in range(len(seq) - len(motif) + 1):
        if seq[i:i+len(motif)] == motif:
            positions.append(i+1)  # 1-based indexing
    return positions

# ----------------------------
# 6. Main Program
# ----------------------------
def main():
    fasta_file = input("Enter path to FASTA file: ")  # Example: "sequence.fasta"
    dna_seq = read_fasta(fasta_file)

    print("\n--- DNA Sequence Analysis ---")
    print("Sequence:", dna_seq)
    
    # Nucleotide count
    nuc_count = count_nucleotides(dna_seq)
    print("Nucleotide Count:", nuc_count)
    
    # GC content
    gc = gc_content(dna_seq)
    print(f"GC Content: {gc:.2f}%")
    
    # Reverse complement
    rev_comp = reverse_complement(dna_seq)
    print("Reverse Complement:", rev_comp)
    
    # Motif search
    motif = input("Enter motif to search (e.g., ATG): ")
    motif_positions = find_motif(dna_seq, motif)
    print(f"Motif {motif} found at positions:", motif_positions)

# ----------------------------
# Run the program
# ----------------------------
if __name__ == "__main__":
    main()
