import csv
from pathlib import Path

def safe_float(x, default=0.0):
    try: return float(x)
    except: return default

def integrate(inp, outp, typ):
    with Path(inp).open() as f_in, Path(outp).open("w") as f_out:
        r = csv.DictReader(f_in, delimiter="\t")
        cols = r.fieldnames + ["final_integrated_score"]
        f_out.write("\t".join(cols) + "\n")
        for row in r:
            base = safe_float(row.get("score_final", "0"))
            dG = safe_float(row.get("rnafold_dG", "0"))
            homo = row.get("has_homopolymer_5+", "no") == "yes"
            rep = safe_float(row.get("simple_repeat_score", "0"))
            acc = safe_float(row.get("accessibility_score", "1"))
            iso = int(row.get("isoform_conservation", "1"))

            score = base
            # Homopolymer penalty
            if homo: score -= 3.0
            # Repeat penalty
            score -= 0.5 * rep
            # Accessibility multiplier
            if acc < 0.7: score *= 0.8
            # RNAfold penalty if extremely stable
            if dG < -20: score -= 2.0
            # Bonus for isoform conservation
            if iso >= 2: score += 1.0
            if iso >= 3: score += 0.5

            row["final_integrated_score"] = f"{score:.2f}"
            f_out.write("\t".join([row[c] for c in cols]) + "\n")

integrate("outputs/design/ASO_candidates_with_conservation.tsv", "outputs/design/ASO_candidates_final_integrated.tsv", "ASO")
integrate("outputs/design/Cas13_guides_with_conservation.tsv", "outputs/design/Cas13_guides_final_integrated.tsv", "Cas13")
print("Wrote outputs/design/ASO_candidates_final_integrated.tsv and outputs/design/Cas13_guides_final_integrated.tsv")
