import csv
from pathlib import Path
from collections import defaultdict

aso_top = Path("outputs/design/ASO_candidates_top.tsv")
cas_top = Path("outputs/design/Cas13_guides_top.tsv")
aso_blast = Path("outputs/results/ASO_offtargets_grch38.txt")
cas_blast = Path("outputs/results/Cas13_offtargets_grch38.txt")

aso_final = Path("outputs/design/ASO_candidates_final_grch38.tsv")
cas_final = Path("outputs/design/Cas13_guides_final_grch38.tsv")

def load_top(path, typ):
    items = []
    with path.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            qid = f"{typ}_{row['transcript_id']}_{row['win_start']}_{row['win_end']}"
            row['qid'] = qid
            items.append(row)
    return items

def count_hits(blast_path):
    hits = defaultdict(int)
    with blast_path.open() as f:
        for line in f:
            if not line.strip(): continue
            qid = line.split("\t", 1)[0]
            hits[qid] += 1
    return hits

def finalize(items, hits, penalty=3.0):
    out = []
    for row in items:
        base = float(row['score'])
        off = hits.get(row['qid'], 0)
        final = base - penalty*off
        out.append({
            "type": row["type"],
            "transcript_id": row["transcript_id"],
            "chr": row["chr"],
            "strand": row["strand"],
            "win_start": row["win_start"],
            "win_end": row["win_end"],
            "gc": row["gc"],
            "tm": row["tm"],
            "score_base": row["score"],
            "offtargets_genome": str(off),
            "score_final": f"{final:.2f}",
            "window_seq": row["window_seq"]
        })
    out.sort(key=lambda r: float(r["score_final"]), reverse=True)
    return out

aso_items = load_top(aso_top, "ASO")
cas_items = load_top(cas_top, "Cas13")
aso_hits = count_hits(aso_blast)
cas_hits = count_hits(cas_blast)

aso_final_rows = finalize(aso_items, aso_hits, penalty=3.0)
cas_final_rows = finalize(cas_items, cas_hits, penalty=3.0)

def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("type\ttranscript_id\tchr\tstrand\twin_start\twin_end\tgc\ttm\tscore_base\tofftargets_genome\tscore_final\twindow_seq\n")
        for r in rows[:50]:
            f.write("\t".join([r[k] for k in ["type","transcript_id","chr","strand","win_start","win_end","gc","tm","score_base","offtargets_genome","score_final","window_seq"]])+"\n")

write(aso_final, aso_final_rows)
write(cas_final, cas_final_rows)
print(f"Wrote {aso_final} and {cas_final} with genome-wide off-target penalties.")
