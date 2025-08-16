# DNA Sequence Analyzer

# Genetic code dictionary for translation
genetic_code = {
    "AUG":"M", "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
    "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"
}

def transcribe_dna_to_rna(dna_seq):
    """Convert DNA sequence to mRNA sequence"""
    return dna_seq.replace("T", "U")

def translate_rna_to_protein(rna_seq):
    """Convert RNA sequence to protein sequence"""
    protein = []
    # Read codons in triplets
    for i in range(0, len(rna_seq)-2, 3):
        codon = rna_seq[i:i+3]
        amino_acid = genetic_code.get(codon, "?")
        if amino_acid == "*":  # stop codon
            break
        protein.append(amino_acid)
    return "".join(protein)

def gc_content(dna_seq):
    """Calculate GC percentage in DNA"""
    g = dna_seq.count("G")
    c = dna_seq.count("C")
    return (g + c) / len(dna_seq) * 100

# Example usage
if __name__ == "__main__":
    dna = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    print("DNA:", dna)

    rna = transcribe_dna_to_rna(dna)
    print("RNA:", rna)

    protein = translate_rna_to_protein(rna)
    print("Protein:", protein)

    gc = gc_content(dna)
    print("GC Content:", f"{gc:.2f}%")
