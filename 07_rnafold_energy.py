import csv
from pathlib import Path

def parse_rnafold(txt_path):
    energies = {}
    with txt_path.open() as f:
        lines = [l.rstrip() for l in f if l.strip()]
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            qid = lines[i][1:]
            # Expect sequence at i+1 and structure+energy at i+2
            if i + 2 < len(lines):
                struct = lines[i+2]
                try:
                    energy_str = struct.split("(")[-1].split(")")[0].strip()
                    energies[qid] = float(energy_str)
                except:
                    energies[qid] = None
            i += 3
        else:
            i += 1
    return energies

def merge_energy(tsv_in, tsv_out, energies, typ):
    inp = Path(tsv_in); outp = Path(tsv_out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with inp.open() as f_in, outp.open("w") as f_out:
        r = csv.DictReader(f_in, delimiter="\t")
        cols = r.fieldnames + ["rnafold_dG"]
        f_out.write("\t".join(cols) + "\n")
        for row in r:
            qid = f"{typ}_{row['transcript_id']}_{row['win_start']}_{row['win_end']}"
            dG = energies.get(qid, None)
            row["rnafold_dG"] = "" if dG is None else f"{dG:.2f}"
            f_out.write("\t".join([row[c] for c in cols]) + "\n")

# Inputs/outputs
aso_tsv = "outputs/design/ASO_candidates_final_grch38.tsv"
cas_tsv = "outputs/design/Cas13_guides_final_grch38.tsv"
aso_txt = Path("outputs/results/ASO_rnafold.txt")
cas_txt = Path("outputs/results/Cas13_rnafold.txt")

# Parse and merge
aso_E = parse_rnafold(aso_txt)
cas_E = parse_rnafold(cas_txt)

merge_energy(aso_tsv, "outputs/design/ASO_candidates_with_dG.tsv", aso_E, "ASO")
merge_energy(cas_tsv, "outputs/design/Cas13_guides_with_dG.tsv", cas_E, "Cas13")

print("Wrote outputs/design/ASO_candidates_with_dG.tsv and outputs/design/Cas13_guides_with_dG.tsv")
