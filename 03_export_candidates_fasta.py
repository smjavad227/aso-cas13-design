import argparse, csv
p = argparse.ArgumentParser()
p.add_argument("--input", required=True)
p.add_argument("--seqcol", default="window_seq")
p.add_argument("--idfmt", default="{type}|{transcript_id}|{win_start}|{win_end}")
p.add_argument("--out", required=True)
a = p.parse_args()
with open(a.input) as f, open(a.out, "w") as o:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        seq = row[a.seqcol].strip().upper()
        _id = a.idfmt.format(**row)
        o.write(f">{_id}\n{seq}\n")
