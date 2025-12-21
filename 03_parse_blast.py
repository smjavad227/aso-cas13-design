import csv
from pathlib import Path
from collections import defaultdict

# Inputs
aso_top = Path("outputs/design/ASO_candidates_top.tsv")
cas_top = Path("outputs/design/Cas13_guides_top.tsv")
aso_blast = Path("outputs/results/ASO_offtargets.txt")
cas_blast = Path("outputs/results/Cas13_offtargets.txt")

# Outputs
aso_final = Path("outputs/design/ASO_candidates_final.tsv")
cas_final = Path("outputs/design/Cas13_guides_final.tsv")

# Load top candidates into dicts keyed by query id used in FASTA
def load_top(path, typ):
    items = []
    with path.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            # Query ID format used in FASTA: >type_transcript_id_win_start_win_end
            qid = f"{typ}_{row['transcript_id']}_{row['win_start']}_{row['win_end']}"
            row['qid'] = qid
            items.append(row)
    return items

aso_items = load_top(aso_top, "ASO")
cas_items = load_top(cas_top, "Cas13")

# Count BLAST hits per qid
def count_hits(blast_path):
    counts = defaultdict(int)
    with blast_path.open() as f:
        for line in f:
            if not line.strip(): continue
            cols = line.strip().split("\t")
            qid = cols[0]
            counts[qid] += 1
    return counts

aso_hits = count_hits(aso_blast)
cas_hits = count_hits(cas_blast)

# Apply off-target penalty and sort
def finalize(items, hits, base_weight=1.0, penalty=2.0):
    out = []
    for row in items:
        base_score = float(row['score'])
        off = hits.get(row['qid'], 0)
        final_score = base_weight*base_score - penalty*off
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
            "offtargets_chr19": str(off),
            "score_final": f"{final_score:.2f}",
            "window_seq": row["window_seq"]
        })
    out.sort(key=lambda r: float(r["score_final"]), reverse=True)
    return out

aso_final_rows = finalize(aso_items, aso_hits, base_weight=1.0, penalty=2.0)
cas_final_rows = finalize(cas_items, cas_hits, base_weight=1.0, penalty=2.0)

# Write outputs
def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("type\ttranscript_id\tchr\tstrand\twin_start\twin_end\tgc\ttm\tscore_base\tofftargets_chr19\tscore_final\twindow_seq\n")
        for r in rows[:50]:  # keep top 50 for presentation
            f.write("\t".join([r[k] for k in ["type","transcript_id","chr","strand","win_start","win_end","gc","tm","score_base","offtargets_chr19","score_final","window_seq"]])+"\n")

write_tsv(aso_final, aso_final_rows)
write_tsv(cas_final, cas_final_rows)
print(f"Wrote {aso_final} and {cas_final} with off-target penalties applied.")
