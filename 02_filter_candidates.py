import csv, math
from pathlib import Path

inp_csv = Path("outputs/results/01_target_windows.csv")
aso_out = Path("outputs/design/ASO_candidates_top.tsv")
cas_out = Path("outputs/design/Cas13_guides_top.tsv")
aso_out.parent.mkdir(parents=True, exist_ok=True)

def gc_content(seq):
    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)

def tm_proxy(seq):
    # Wallace rule: 2*(A+T) + 4*(G+C)
    return 2*(seq.count("A")+seq.count("T")) + 4*(seq.count("G")+seq.count("C"))

rows = []
with inp_csv.open() as f:
    reader = csv.DictReader(f)
    for r in reader:
        seq = r["window_seq"]
        gc = gc_content(seq)
        tm = tm_proxy(seq)
        score = 0.4*gc + 0.3*tm - 0.1*len(seq)  # placeholder scoring
        r["gc"] = f"{gc:.1f}"
        r["tm"] = f"{tm:.1f}"
        r["score"] = f"{score:.2f}"
        rows.append(r)

# Split ASO vs Cas13
aso_rows = [r for r in rows if r["type"]=="ASO"]
cas_rows = [r for r in rows if r["type"]=="Cas13"]

# Sort by score descending
aso_rows.sort(key=lambda r: float(r["score"]), reverse=True)
cas_rows.sort(key=lambda r: float(r["score"]), reverse=True)

# Write top 50 each
with aso_out.open("w") as f:
    f.write("type\ttranscript_id\tchr\tstrand\twin_start\twin_end\tgc\ttm\tscore\twindow_seq\n")
    for r in aso_rows[:50]:
        f.write("\t".join([r["type"],r["transcript_id"],r["chr"],r["strand"],
                           r["win_start"],r["win_end"],r["gc"],r["tm"],r["score"],r["window_seq"]])+"\n")

with cas_out.open("w") as f:
    f.write("type\ttranscript_id\tchr\tstrand\twin_start\twin_end\tgc\ttm\tscore\twindow_seq\n")
    for r in cas_rows[:50]:
        f.write("\t".join([r["type"],r["transcript_id"],r["chr"],r["strand"],
                           r["win_start"],r["win_end"],r["gc"],r["tm"],r["score"],r["window_seq"]])+"\n")

print(f"Wrote {aso_out} and {cas_out} with top candidates.")
