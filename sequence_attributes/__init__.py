"""
Public exports for sequence_attributes package.
"""

from .sequence_formats.fasta_format import get_fasta_lists, _verify_lists
from .utils.io_utils import FileHandler
from .utils.seq_attribute_utils import (
    gc_content,
    lookup_by_ccds,
    get_tm_from_dna_sequence,
    get_sequence_composition,
    calculate_amino_acid_content,
    extract_kmers,
    return_standard_genetic_code,
    protein_translation,
    get_additional_sequence_attributes,
)

__all__ = [
    "get_fasta_lists",
    "_verify_lists",
    "FileHandler",
    "gc_content",
    "lookup_by_ccds",
    "get_tm_from_dna_sequence",
    "get_sequence_composition",
    "calculate_amino_acid_content",
    "extract_kmers",
    "return_standard_genetic_code",
    "protein_translation",
    "get_additional_sequence_attributes",
]

