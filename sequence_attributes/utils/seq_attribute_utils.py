from typing import Union
import pandas as pd
from collections import Counter, namedtuple


def gc_content(sequence: str, round_to: int = 2, percentage: bool = True) -> float:
    """Return GC content as % or fraction."""
    seq = sequence.upper()
    assert len(seq) > 0, "Sequence must not be empty."
    valid = set("ACGTN")
    assert set(seq).issubset(valid), "Sequence contains invalid characters."
    gc = seq.count("G") + seq.count("C")
    frac = gc / len(seq)
    return round(frac * 100, round_to) if percentage else round(frac, round_to)


def lookup_by_ccds(ccds_id: str, chrom: str, df: pd.DataFrame) -> pd.DataFrame:
    """Filter rows matching CCDS ID AND chromosome."""
    return df[(df["ccds_id"] == ccds_id) & (df["chrom"] == chrom)]


def get_tm_from_dna_sequence(sequence: str, round_to: int = 2) -> float:
    """Calculate melting temperature using simple rule."""
    seq = sequence.upper()
    assert len(seq) > 0, "Sequence must not be empty."
    valid = set("ACGTN")
    assert set(seq).issubset(valid), "Invalid characters in sequence."
    g = seq.count("G")
    c = seq.count("C")
    a = seq.count("A")
    t = seq.count("T")
    tm = 4 * (g + c) + 2 * (a + t)
    return round(tm, round_to)


def get_sequence_composition(sequence: Union[str, list, tuple]) -> dict:
    """Return count of each element."""
    if isinstance(sequence, str):
        seq = list(sequence)
    else:
        seq = list(sequence)
    return dict(Counter(seq))


def calculate_amino_acid_content(amino_acid: str,
                                 sequence: str,
                                 round_to: int = 2,
                                 percentage: bool = True) -> float:
    """Return % or fraction of a given amino acid in sequence."""
    seq = sequence.upper()
    aa = amino_acid.upper()
    assert len(aa) == 1, "Amino acid must be one character."
    assert len(seq) > 0, "Sequence must not be empty."
    fraction = seq.count(aa) / len(seq)
    return round(fraction * 100, round_to) if percentage else round(fraction, round_to)


def extract_kmers(sequence: str, kmer_size: int = 3) -> list:
    """Extract all kmers of length kmer_size."""
    seq = sequence.upper()
    assert kmer_size > 0
    if len(seq) < kmer_size:
        return []
    return [seq[i:i + kmer_size] for i in range(len(seq) - kmer_size + 1)]


def return_standard_genetic_code() -> dict:
    """Return the standard genetic code."""
    return {
        "UUU": "F", "UUC": "F",
        "UUA": "L", "UUG": "L",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "UAU": "Y", "UAC": "Y",
        "UAA": None, "UAG": None,
        "UGU": "C", "UGC": "C",
        "UGA": None,
        "UGG": "W",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAU": "H", "CAC": "H",
        "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AUU": "I", "AUC": "I", "AUA": "I",
        "AUG": "M",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAU": "N", "AAC": "N",
        "AAA": "K", "AAG": "K",
        "AGU": "S", "AGC": "S",
        "AGA": "R", "AGG": "R",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAU": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }


def protein_translation(sequence: str, genetic_code: dict) -> str:
    """Translate DNA → protein until stop codon."""
    rna = sequence.upper().replace("T", "U")
    protein = []
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        if len(codon) != 3:
            break
        aa = genetic_code.get(codon)
        if aa is None:
            break
        protein.append(aa)
    return "".join(protein)


# Namedtuple required by assignment
FastaTuple = namedtuple(
    "FastaTuple",
    [
        "header",
        "dna_sequence",
        "protein_sequence",
        "protein_sequence_len",
        "tm",
        "gc",
        "nucleotide_comp",
        "amino_acid_comp",
        "kmer_comp",
        "proline_comp",
        "ccds",
        "refseq_gene_id",
        "biotype",
        "ensembl_gene_id",
        "description",
    ]
)


def get_additional_sequence_attributes(header: str,
                                       dna_sequence: str,
                                       genetic_code: dict,
                                       attribute_df: pd.DataFrame) -> FastaTuple:
    """Generate all required attributes for one FASTA entry."""

    protein = protein_translation(dna_sequence, genetic_code)
    tm_val = get_tm_from_dna_sequence(dna_sequence)
    gc_val = gc_content(dna_sequence)
    nt_comp = get_sequence_composition(dna_sequence)
    aa_comp = get_sequence_composition(protein)
    kmers = extract_kmers(dna_sequence, 3)
    kmer_comp = get_sequence_composition(kmers)
    proline = calculate_amino_acid_content("P", protein)
    row = attribute_df.loc[attribute_df["ccds_id"] == header]
    refseq = row["refseq_gene_id"].iloc[0] if not row.empty else None
    biotype = row["biotype"].iloc[0] if not row.empty else None
    ensg = row["ensembl_gene_id"].iloc[0] if not row.empty else None
    desc = row["description"].iloc[0] if not row.empty else None

    return FastaTuple(
        header,
        dna_sequence,
        protein,
        len(protein),
        tm_val,
        gc_val,
        nt_comp,
        aa_comp,
        kmer_comp,
        proline,
        header,
        refseq,
        biotype,
        ensg,
        desc,
    )

