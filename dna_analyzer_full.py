# DNA Analyzer with Disease Association
# Author: SANA
# Date: 2025-08-15
# Description: Reads DNA sequences from FASTA, performs analysis, motif search, and disease association

# ----------------------------
# 1. Function to read FASTA file
# ----------------------------
def read_fasta(file_path):
    sequence = ""
    with open(file_path, "r") as f:
        for line in f:
            if not line.startswith(">"):  # Ignore header lines
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
# 6. Function to check disease-associated mutations
# ----------------------------
def find_disease(seq, mutation_db):
    found = {}
    for motif, disease in mutation_db.items():
        if motif in seq:
            found[motif] = disease
    return found

# ----------------------------
# 7. Main Program
# ----------------------------
def main():
    # Input FASTA file
    fasta_file = input("Enter path to FASTA file: ")
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
    if motif_positions:
        print(f"Motif {motif} found at positions:", motif_positions)
    else:
        print(f"Motif {motif} not found in sequence.")
    
    # Disease-associated mutation search
    # Example database (replace with real mutations for actual use)
    mutation_db = {
        "ATGCGT": "Disease A - Example mutation in gene X",
        "CGTTAG": "Disease B - Example mutation in gene Y",
        "GGCCTA": "Disease C - Example mutation in gene Z"
    }
    disease_hits = find_disease(dna_seq, mutation_db)
    if disease_hits:
        print("\n--- Potential Disease Associations ---")
        for motif, disease in disease_hits.items():
            print(f"Motif: {motif} → {disease}")
    else:
        print("\nNo known disease-associated mutations found.")

# ----------------------------
# Run the program
# ----------------------------
if __name__ == "__main__":
    main()
