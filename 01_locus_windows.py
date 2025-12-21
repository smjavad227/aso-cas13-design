import csv
from pathlib import Path

ref_fa = Path("inputs/reference/chr19.fa")
gtf = Path("inputs/reference/annotation.gtf")
bed = Path("inputs/reference/DMPK_CTGrepeat_GRCh38.bed")
out_csv = Path("outputs/results/01_target_windows.csv")
out_csv.parent.mkdir(parents=True, exist_ok=True)

def norm_chr(ch):
    ch = ch.strip()
    if ch.startswith("chr"):
        ch = ch[3:]
    return ch  # e.g., "chr19" -> "19"

# Load chr19 sequence (Ensembl uses "19")
seq_lines = []
with ref_fa.open() as f:
    for line in f:
        if line.startswith(">"): continue
        seq_lines.append(line.strip())
chr19_seq = "".join(seq_lines).upper()

def rc(s):
    comp = str.maketrans("ACGTN", "TGCAN")
    return s.translate(comp)[::-1]

# Parse GTF for DMPK exons/UTRs
regions = []  # (chr, start, end, strand, transcript_id)
with gtf.open() as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9: continue
        chrom, src, feature, start, end, score, strand, frame, attrs = parts
        if feature not in ("three_prime_utr", "exon"):
            continue

        # parse attributes
        at = {}
        for kv in attrs.split(";"):
            kv = kv.strip()
            if not kv: continue
            if " " in kv:
                k, v = kv.split(" ", 1)
                at[k] = v.strip().strip('"')
        gene_name = at.get("gene_name") or at.get("gene")
        tid = at.get("transcript_id")
        if gene_name == "DMPK" and tid:
            regions.append((norm_chr(chrom), int(start), int(end), strand, tid))

# Load BED locus
with bed.open() as f:
    b = f.readline().strip().split("\t")
    bed_chr, bed_start, bed_end = norm_chr(b[0]), int(b[1]), int(b[2])

# Intersect BED with regions
targets = []
for chrom, start, end, strand, tid in regions:
    if chrom != bed_chr:
        continue
    ov_start = max(start, bed_start)
    ov_end = min(end, bed_end)
    if ov_start < ov_end:
        s = chr19_seq[ov_start-1:ov_end]  # 1-based -> 0-based
        if strand == "-":
            s = rc(s)
        targets.append({
            "transcript_id": tid, "chr": chrom, "strand": strand,
            "start": ov_start, "end": ov_end, "sequence": s
        })

# Windowing
ASO_W = 18
CAS_W = 28
rows = []
for t in targets:
    seq = t["sequence"]
    for i in range(0, max(0, len(seq) - ASO_W + 1)):
        wseq = seq[i:i+ASO_W]
        rows.append([
            "ASO", t["transcript_id"], t["chr"], t["strand"],
            t["start"], t["end"], t["start"]+i, t["start"]+i+ASO_W-1,
            t["sequence"], wseq
        ])
    for i in range(0, max(0, len(seq) - CAS_W + 1)):
        wseq = seq[i:i+CAS_W]
        rows.append([
            "Cas13", t["transcript_id"], t["chr"], t["strand"],
            t["start"], t["end"], t["start"]+i, t["start"]+i+CAS_W-1,
            t["sequence"], wseq
        ])

with out_csv.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["type","transcript_id","chr","strand","start","end","win_start","win_end","sequence","window_seq"])
    w.writerows(rows)

print(f"Wrote {out_csv} with {len(rows)} window rows.")
