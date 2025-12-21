import csv
from pathlib import Path

def accessibility(gc, homopolymer_flag):
    # Heuristic: start from 1.0, penalize GC>70 and homopolymers
    score = 1.0
    try:
        gc_val = float(gc)
    except:
        gc_val = None
    if gc_val is not None:
        if gc_val > 70: score -= 0.3
        if gc_val > 80: score -= 0.2
    if homopolymer_flag == "yes":
        score -= 0.3
    return max(0.0, round(score, 2))

def add_access(inp, outp):
    with Path(inp).open() as f_in, Path(outp).open("w") as f_out:
        r = csv.DictReader(f_in, delimiter="\t")
        cols = r.fieldnames + ["accessibility_score"]
        f_out.write("\t".join(cols) + "\n")
        for row in r:
            gc = row.get("gc")
            homo = row.get("has_homopolymer_5+", "no")
            row["accessibility_score"] = str(accessibility(gc, homo))
            f_out.write("\t".join([row[c] for c in cols]) + "\n")

add_access("outputs/design/ASO_candidates_with_repeats.tsv", "outputs/design/ASO_candidates_with_access.tsv")
add_access("outputs/design/Cas13_guides_with_repeats.tsv", "outputs/design/Cas13_guides_with_access.tsv")
print("Wrote outputs/design/*_with_access.tsv")
