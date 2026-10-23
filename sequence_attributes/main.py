import argparse
import sys
import pandas as pd

from sequence_attributes import (
    get_fasta_lists,
    get_additional_sequence_attributes,
    return_standard_genetic_code,
)


def get_cli_args():
    """Command-line interface for Assignment 5."""
    parser = argparse.ArgumentParser(description="Assignment 5 CLI")

    parser.add_argument(
        "--infile_ccds_fasta",
        required=True,
        help="Input CCDS FASTA file",
    )

    parser.add_argument(
        "--infile_ccds_attributes",
        required=True,
        help="Input CCDS attributes file",
    )

    parser.add_argument(
        "--infile_ensembl_gene",
        required=True,
        help="Input Ensembl gene TSV file",
    )

    parser.add_argument(
        "--excel_outfile",
        required=True,
        help="Output Excel file",
    )

    return parser.parse_args()


def main():
    args = get_cli_args()

    # Load CCDS attributes
    ccds_df = pd.read_csv(args.infile_ccds_attributes, sep="\t")
    ccds_df = ccds_df.rename(
        columns={"gene_id": "refseq_gene_id", "#chromosome": "chrom"}
    )

    # Load Ensembl data
    ens_df = pd.read_csv(args.infile_ensembl_gene, sep="\t")
    ens_df = ens_df.rename(
        columns={
            "ID": "ensembl_gene_id",
            "Canonical Transcript": "ensembl_canonical_transcript_id",
            "Gene": "gene",
            "Description": "description",
            "Biotype": "biotype",
        }
    )

    # Merge CCDS + Ensembl on gene
    merged_df = pd.merge(
        ccds_df, ens_df, how="left", on="gene"
    )

    # Parse FASTA
    headers, sequences = get_fasta_lists(args.infile_ccds_fasta)

    genetic_code = return_standard_genetic_code()
    results = []

    # Process sequences
    for i, (header, dna_sequence) in enumerate(zip(headers, sequences)):
        attributes = get_additional_sequence_attributes(
            header=header,
            dna_sequence=dna_sequence,
            genetic_code=genetic_code,
            attribute_df=merged_df,
        )
        results.append(attributes)

        if i == 50:  # TEMP LIMIT FOR TESTING — CHANGE LATER
            break

    # Convert results → DataFrame
    df = pd.DataFrame(results)

    # Sort
    df = df.sort_values(
        by=["proline_comp", "protein_sequence_len"],
        ascending=[False, True]
    )

    # Full Excel output
    df.to_excel(args.excel_outfile, index=False)

    # Top/bottom 10 TSV
    tsv_file = args.excel_outfile.replace(".xlsx", ".tsv")
    top_bottom = pd.concat([df.head(10), df.tail(10)])
    cols = [
        "ccds", "refseq_gene_id", "biotype",
        "ensembl_gene_id", "description",
        "protein_sequence_len", "proline_comp",
    ]
    top_bottom[cols].to_csv(tsv_file, sep="\t", index=False)


if __name__ == "__main__":
    main()

