import pandas as pd
from pathlib import Path

aso = pd.read_csv("outputs/design/ASO_candidates_final_integrated.tsv", sep="\t")
cas = pd.read_csv("outputs/design/Cas13_guides_final_integrated.tsv", sep="\t")

aso_strict = aso[(aso["offtargets_genome"] == 0) & (aso["has_homopolymer_5+"] == "no")].sort_values("final_integrated_score", ascending=False).head(20)
cas_strict = cas[(cas["offtargets_genome"] <= 1) & (cas["has_homopolymer_5+"] == "no")].sort_values("final_integrated_score", ascending=False).head(20)

Path("outputs/design").mkdir(parents=True, exist_ok=True)
aso_strict.to_csv("outputs/design/ASO_candidates_final_stringent.tsv", sep="\t", index=False)
cas_strict.to_csv("outputs/design/Cas13_guides_final_stringent.tsv", sep="\t", index=False)

print("Wrote stringent final sets.")
