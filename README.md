# DNA and Protein Sequence Analysis (Python)

Academic bioinformatics project from Northeastern University. Processes CCDS FASTA nucleotide sequences, joins CCDS attributes with Ensembl gene annotations, and calculates sequence features including GC content, estimated DNA melting temperature, translated protein length, amino-acid composition (including proline), and k-mers.

## Tools
Python, pandas, openpyxl, pytest.

## Structure
- `sequence_attributes/main.py`: command-line pipeline and Excel/TSV exports
- `sequence_attributes/utils/`: sequence metrics and file utilities
- `sequence_attributes/sequence_formats/`: FASTA parsing
- `sequence_attributes/tests/`: unit tests

## Setup
```bash
python -m pip install -r requirements.txt
```

## Run
From the repository root, supply downloaded CCDS FASTA, CCDS attributes, and Ensembl gene TSV inputs:
```bash
python -m sequence_attributes.main \
  --infile_ccds_fasta path/to/CCDS_nucleotide.current.fna \
  --infile_ccds_attributes path/to/CCDS.current.txt \
  --infile_ensembl_gene path/to/ensembl_gene_data.tsv \
  --excel_outfile output.xlsx
```

The program produces an Excel workbook and a TSV summary. **Note:** the original submitted `main.py` contains a temporary `if i == 50: break` limit, so it processes only the first 51 FASTA records. Remove this limit after validating full-dataset performance if you want to analyze every record.

## Tests
```bash
python -m pytest sequence_attributes/tests
```

## Data
The original input datasets are intentionally excluded from this repository package because they are large. Obtain data from their respective CCDS/Ensembl sources, and review applicable licensing and university course policies before public sharing. The project is provided as coursework and is not presented as a validated production bioinformatics pipeline.
