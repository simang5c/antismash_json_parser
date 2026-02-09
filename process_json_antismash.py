import json
import re
import pandas as pd
import sys

def parse_location(loc_str):
    """
    Parse feature location string.
    """
    strand_match = re.search(r'\((\+|-)\)$', loc_str)
    strand = strand_match.group(1) if strand_match else '+'

    # extract all the coordinate
    coords = [list(map(int, c.split(':'))) for c in re.findall(r'\[(\d+:\d+)\]', loc_str)]
    starts = [s for s, e in coords]
    ends = [e for s, e in coords]
    return min(starts), max(ends), strand

# Load your json file
json_file = sys.argv[1]
with open(json_file) as f:
    data = json.load(f)

region_rows = []
gene_rows = []

for record in data.get("records", []):
    genome_id = record.get("id")
    features = record.get("features", [])
    areas = record.get("areas", [])

    if not areas:
        print(f"No BGC areas identified {genome_id}")
        continue

    for i, area in enumerate(areas, start=1):
        # BGC area info
        r_start = area["start"]
        r_end = area["end"]
        products = ", ".join(area.get("products", []))
        protoclusters = ", ".join([p.get("category", "") for p in area.get("protoclusters", {}).values()])

        region_rows.append({
            "genome": genome_id,
            "area_number": i,
            "start": r_start,
            "end": r_end,
            "products": products,
            "protoclusters": protoclusters
        })

        # Genes that are overlapping CDS regions
        # When uploading your own fasta file, you must also upload a gff file that should contain CDS regions
        for feature in features:
            if feature.get("type") != "CDS":
                continue

            loc_str = feature["location"]
            g_start, g_end, strand = parse_location(loc_str)

            if g_start >= r_start and g_end <= r_end:
                q = feature.get("qualifiers", {})
                gene_rows.append({
                    "genome": genome_id, # your genome id
                    "area_number": i, 
                    "locus_tag": q.get("gene", [""])[0],
                    "gene_id": q.get("gene_id", [""])[0],
                    "gene_name": q.get("gene_name", [""])[0],
                    "transcript_id": q.get("transcript_id", [""])[0],
                    "parent_transcript": q.get("Parent", [""])[0],
                    "source": q.get("source", [""])[0],
                    "score": q.get("score", [""])[0],
                    "start": g_start,
                    "end": g_end,
                    "strand": strand,
                    "product": q.get("product", [""])[0],
                    "translation": q.get("translation", [""])[0]
                })

# save all results in a pddataframe
pd.DataFrame(region_rows).to_csv("antismash_areas.csv", index=False)
pd.DataFrame(gene_rows).to_csv("antismash_area_genes.csv", index=False)

print(f"Saved {len(region_rows)} BGC areas and {len(gene_rows)} genes.")
