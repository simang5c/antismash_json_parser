# AntiSMASH JSON Parser for BGCs and Genes

## Description

This Python program extracts Biosynthetic Gene Clusters (BGCs) and their associated genes from antiSMASH JSON output files. It handles multi-exon and single-exon CDS features, retrieves strand information, and captures all relevant gene annotations, including gene IDs, names, transcripts, products, and translations. The output is formatted into two CSV files for easy downstream analysis.

## Requirements

- Python 3.8 or higher
 
  - Pandas
  - json
  - re
  - sys

## Usage


Place your antiSMASH JSON file (e.g., antismash_output.json) in the same folder as the script.

## Run the parser:
```
# Change the json file with your file
python process_json_antismash.py Fungi_genome_antismash_downloaded.json
```

The program will generate two CSV files:
- antismash_areas.csv – contains BGC coordinates, products, and protocluster types
- antismash_area_genes.csv – contains gene coordinates, strand, annotations, and translation sequences
