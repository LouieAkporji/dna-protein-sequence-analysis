from sequence_attributes import (
    gc_content,
    get_sequence_composition,
    calculate_amino_acid_content,
    extract_kmers,
    return_standard_genetic_code,
    protein_translation,
    get_tm_from_dna_sequence,
    lookup_by_ccds,
    get_additional_sequence_attributes,
)


def test_gc_content_basic():
    assert gc_content("AGCTCGGTA") == 55.56

