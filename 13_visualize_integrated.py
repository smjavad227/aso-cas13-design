import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

fig_dir = Path("outputs/figures")
fig_dir.mkdir(parents=True, exist_ok=True)

aso = pd.read_csv("outputs/design/ASO_candidates_final_integrated.tsv", sep="\t")
cas = pd.read_csv("outputs/design/Cas13_guides_final_integrated.tsv", sep="\t")

for df in (aso, cas):
    for col in ["gc","tm","score_base","offtargets_genome","score_final","rnafold_dG","simple_repeat_score","accessibility_score","isoform_conservation","final_integrated_score"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

# 1) Integrated score vs off-targets
plt.figure(figsize=(9,5))
sns.scatterplot(x="offtargets_genome", y="final_integrated_score", data=aso, color="#1f77b4", s=40, label="ASO")
sns.scatterplot(x="offtargets_genome", y="final_integrated_score", data=cas, color="#d62728", s=40, label="Cas13")
plt.xlabel("Genome-wide off-target count")
plt.ylabel("Final integrated score")
plt.title("Integrated score vs genome-wide off-targets")
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir / "integrated_score_vs_offtargets.png", dpi=300)
plt.close()

# 2) ΔG distributions
plt.figure(figsize=(8,5))
sns.histplot(aso["rnafold_dG"], color="#1f77b4", bins=20, alpha=0.5, label="ASO")
sns.histplot(cas["rnafold_dG"], color="#d62728", bins=20, alpha=0.5, label="Cas13")
plt.xlabel("RNAfold ΔG (kcal/mol)")
plt.ylabel("Count")
plt.title("ΔG distribution (ASO vs Cas13)")
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir / "dG_distribution.png", dpi=300)
plt.close()

# 3) Accessibility distributions
plt.figure(figsize=(8,5))
sns.histplot(aso["accessibility_score"], color="#1f77b4", bins=10, alpha=0.5, label="ASO")
sns.histplot(cas["accessibility_score"], color="#d62728", bins=10, alpha=0.5, label="Cas13")
plt.xlabel("Accessibility score")
plt.ylabel("Count")
plt.title("Accessibility score distribution (ASO vs Cas13)")
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir / "accessibility_distribution.png", dpi=300)
plt.close()

# 4) Top 10 by integrated score (grouped bar)
aso_top10 = aso.sort_values("final_integrated_score", ascending=False).head(10).copy()
cas_top10 = cas.sort_values("final_integrated_score", ascending=False).head(10).copy()
aso_top10["label"] = aso_top10["type"] + "_" + aso_top10["transcript_id"] + "_" + aso_top10["win_start"].astype(str) + "_" + aso_top10["win_end"].astype(str)
cas_top10["label"] = cas_top10["type"] + "_" + cas_top10["transcript_id"] + "_" + cas_top10["win_start"].astype(str) + "_" + cas_top10["win_end"].astype(str)

combo = pd.concat([
    aso_top10[["label","final_integrated_score"]].assign(group="ASO"),
    cas_top10[["label","final_integrated_score"]].assign(group="Cas13")
], ignore_index=True)

plt.figure(figsize=(12,6))
sns.barplot(data=combo, x="label", y="final_integrated_score", hue="group", palette=["#1f77b4", "#d62728"])
plt.xticks(rotation=60, ha="right")
plt.xlabel("Candidate (type_transcript_win_start_win_end)")
plt.ylabel("Final integrated score")
plt.title("Top 10 candidates by integrated score (ASO vs Cas13)")
plt.tight_layout()
plt.savefig(fig_dir / "top_candidates_integrated_barplot.png", dpi=300)
plt.close()

print("Saved figures: integrated_score_vs_offtargets.png, dG_distribution.png, accessibility_distribution.png, top_candidates_integrated_barplot.png")
