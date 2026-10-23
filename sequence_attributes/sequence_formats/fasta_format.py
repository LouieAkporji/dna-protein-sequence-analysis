"""
FASTA parsing utilities.
"""

from typing import Tuple
from sequence_attributes.utils.io_utils import FileHandler


def get_fasta_lists(infile: str = None) -> Tuple[list, list]:
    """
    Parse a FASTA file and return a list of headers and list of sequences.
    """
    headers = []
    sequences = []
    seq = []

    with FileHandler(infile) as infile_handle:
        for line in infile_handle:
            line = line.strip()

            if line.startswith(">"):
                # remove leading >
                header = line[1:].strip()

                if seq:
                    sequences.append("".join(seq))
                    seq = []

                headers.append(header)

            else:
                seq.append(line)

        # append final sequence
        if seq:
            sequences.append("".join(seq))

    _verify_lists(headers, sequences)
    return headers, sequences


def _verify_lists(header_list: list = None, seq_list: list = None) -> bool:
    """
    Ensure header list and sequence list are the same length.
    """
    if len(header_list) != len(seq_list):
        print("Header and sequence list lengths do not match.")
        raise SystemExit(1)
    return True

