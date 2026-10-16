import sys
from collections import Counter

# ---------- Genetic code (RNA codons) ----------
BASES = "UCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: AMINO[16 * i + 4 * j + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)
}

AMINO_NAMES = {
    "F": "Phenylalanine", "L": "Leucine", "S": "Serine", "Y": "Tyrosine",
    "C": "Cysteine", "W": "Tryptophan", "P": "Proline", "H": "Histidine",
    "Q": "Glutamine", "R": "Arginine", "I": "Isoleucine", "M": "Methionine (START)",
    "T": "Threonine", "N": "Asparagine", "K": "Lysine", "V": "Valine",
    "A": "Alanine", "D": "Aspartic acid", "E": "Glutamic acid",
    "G": "Glycine", "*": "STOP",
}


# ---------- 1. Read the FASTA file ----------
def read_fasta(path):
    """Return the genome as one uppercase string (header lines are skipped)."""
    sequence = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sequence:      # a second header means a second record: stop
                    break
                continue
            sequence.append(line.upper())
    return "".join(sequence)


# ---------- 2. Sliding window of size 3 ----------
def count_codons(seq):
    counts = Counter()
    valid = set("ACGT")
    for i in range(len(seq) - 2):          # window moves 1 position at a time
        codon = seq[i:i + 3]
        if set(codon) <= valid:            # ignore windows with N, R, Y, etc.
            counts[codon] += 1
    return counts


# ---------- 3. DNA -> RNA -> amino acid ----------
def dna_to_rna(codon):
    return codon.replace("T", "U")


def translate(rna_codon):
    return CODON_TABLE[rna_codon]


# ---------- 4. Main ----------
def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "genome.fasta"
    seq = read_fasta(path)
    print(f"File: {path}")
    print(f"Genome length: {len(seq):,} bases")

    counts = count_codons(seq)
    total = sum(counts.values())
    print(f"Total codons counted: {total:,}\n")

    print(f"{'Rank':<5}{'DNA':<6}{'RNA':<6}{'Count':<10}{'Freq %':<9}{'AA':<4}Amino acid")
    print("-" * 60)
    for rank, (codon, n) in enumerate(counts.most_common(10), start=1):
        rna = dna_to_rna(codon)
        aa = translate(rna)
        print(f"{rank:<5}{codon:<6}{rna:<6}{n:<10}{n / total * 100:<9.3f}{aa:<4}{AMINO_NAMES[aa]}")


if __name__ == "__main__":
    main()