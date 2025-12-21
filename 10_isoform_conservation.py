import csv
from pathlib import Path
from collections import defaultdict

def isoform_coverage(inp):
    # key by window sequence; count unique transcripts
    seq_to_tx = defaultdict(set)
    rows = []
    with Path(inp).open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            seq_to_tx[row["window_seq"]].add(row["transcript_id"])
            rows.append(row)
    coverage = {s: len(txs) for s, txs in seq_to_tx.items()}
    return rows, coverage

def annotate(inp, outp):
    rows, cov = isoform_coverage(inp)
    cols = list(rows[0].keys()) + ["isoform_conservation"]
    with Path(outp).open("w") as f_out:
        f_out.write("\t".join(cols) + "\n")
        for row in rows:
            row["isoform_conservation"] = str(cov.get(row["window_seq"], 1))
            f_out.write("\t".join([row[c] for c in cols]) + "\n")

annotate("outputs/design/ASO_candidates_with_access.tsv", "outputs/design/ASO_candidates_with_conservation.tsv")
annotate("outputs/design/Cas13_guides_with_access.tsv", "outputs/design/Cas13_guides_with_conservation.tsv")
print("Wrote outputs/design/*_with_conservation.tsv")
