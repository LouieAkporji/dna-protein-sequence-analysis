import os
from sequence_attributes import get_fasta_lists, _verify_lists

def test_get_fasta_lists():
    # Create a small temporary FASTA file
    fasta_content = ">seq1\nATGC\n>seq2\nGGGCCC\n"
    tmp_file = "tmp_test.fasta"
    with open(tmp_file, "w") as f:
        f.write(fasta_content)

    headers, seqs = get_fasta_lists(tmp_file)

    assert headers == ["seq1", "seq2"]
    assert seqs == ["ATGC", "GGGCCC"]

    os.remove(tmp_file)


def test_verify_lists():
    assert _verify_lists(["a", "b"], ["AAA", "TTT"]) is True

