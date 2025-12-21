import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

aso_path = Path("outputs/design/ASO_candidates_final_grch38.tsv")
cas_path = Path("outputs/design/Cas13_guides_final_grch38.tsv")
fig_dir = Path("outputs/figures")
fig_dir.mkdir(parents=True, exist_ok=True)

# Load data
aso = pd.read_csv(aso_path, sep="\t")
cas = pd.read_csv(cas_path, sep="\t")

# Ensure numeric
for df in (aso, cas):
    df["gc"] = pd.to_numeric(df["gc"], errors="coerce")
    df["tm"] = pd.to_numeric(df["tm"], errors="coerce")
    df["score_base"] = pd.to_numeric(df["score_base"], errors="coerce")
    df["offtargets_genome"] = pd.to_numeric(df["offtargets_genome"], errors="coerce")
    df["score_final"] = pd.to_numeric(df["score_final"], errors="coerce")

# 1) GC% distribution (overlaid histograms)
plt.figure(figsize=(8,5))
sns.histplot(aso["gc"], color="#1f77b4", bins=20, alpha=0.5, kde=False, label="ASO")
sns.histplot(cas["gc"], color="#d62728", bins=20, alpha=0.5, kde=False, label="Cas13")
plt.xlabel("GC%")
plt.ylabel("Count")
plt.title("GC% distribution of ASO and Cas13 candidates")
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir / "gc_distribution.png", dpi=300)
plt.close()

# 2) Score vs off-targets scatter
plt.figure(figsize=(9,5))
sns.scatterplot(x="offtargets_genome", y="score_final", data=aso, color="#1f77b4", s=40, label="ASO")
sns.scatterplot(x="offtargets_genome", y="score_final", data=cas, color="#d62728", s=40, label="Cas13")
plt.xlabel("Genome-wide off-target count")
plt.ylabel("Final score")
plt.title("Final score vs genome-wide off-targets")
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir / "score_vs_offtargets.png", dpi=300)
plt.close()

# 3) Top 10 by final score (ASO and Cas13)
aso_top10 = aso.sort_values("score_final", ascending=False).head(10).copy()
cas_top10 = cas.sort_values("score_final", ascending=False).head(10).copy()

aso_top10["label"] = aso_top10["type"] + "_" + aso_top10["transcript_id"] + "_" + aso_top10["win_start"].astype(str) + "_" + aso_top10["win_end"].astype(str)
cas_top10["label"] = cas_top10["type"] + "_" + cas_top10["transcript_id"] + "_" + cas_top10["win_start"].astype(str) + "_" + cas_top10["win_end"].astype(str)

# Combine for a grouped bar chart
combo = pd.concat([
    aso_top10[["label","score_final"]].assign(group="ASO"),
    cas_top10[["label","score_final"]].assign(group="Cas13")
], ignore_index=True)

plt.figure(figsize=(12,6))
sns.barplot(data=combo, x="label", y="score_final", hue="group", palette=["#1f77b4", "#d62728"])
plt.xticks(rotation=60, ha="right")
plt.xlabel("Candidate (type_transcript_win_start_win_end)")
plt.ylabel("Final score")
plt.title("Top 10 candidates by final score (ASO vs Cas13)")
plt.tight_layout()
plt.savefig(fig_dir / "top_candidates_barplot.png", dpi=300)
plt.close()

print("Saved figures: outputs/figures/gc_distribution.png, score_vs_offtargets.png, top_candidates_barplot.png")
