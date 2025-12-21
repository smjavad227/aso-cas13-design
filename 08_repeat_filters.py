import re, csv
from pathlib import Path

def has_homopolymer(seq, k=5):
    return any(n*k in seq for n in "ATGC")

def repeat_score(seq):
    # Count simple di-/tri-nucleotide repeats (rough heuristic)
    di = sum(seq.count(d*3) for d in ["AT","TA","CG","GC","AG","GA","TC","CT","AC","CA","TG","GT"])
    tri = sum(seq.count(t*2) for t in ["ATG","CAG","CTG","GAA","GTT","TGG","CCC","GGG","AAA","TTT"])
    return di + tri

def annotate(inp, outp):
    with Path(inp).open() as f_in, Path(outp).open("w") as f_out:
        r = csv.DictReader(f_in, delimiter="\t")
        cols = r.fieldnames + ["has_homopolymer_5+", "simple_repeat_score"]
        f_out.write("\t".join(cols) + "\n")
        for row in r:
            s = row["window_seq"]
            row["has_homopolymer_5+"] = "yes" if has_homopolymer(s, 5) else "no"
            row["simple_repeat_score"] = str(repeat_score(s))
            f_out.write("\t".join([row[c] for c in cols]) + "\n")

annotate("outputs/design/ASO_candidates_with_dG.tsv", "outputs/design/ASO_candidates_with_repeats.tsv")
annotate("outputs/design/Cas13_guides_with_dG.tsv", "outputs/design/Cas13_guides_with_repeats.tsv")
print("Wrote outputs/design/*_with_repeats.tsv")
