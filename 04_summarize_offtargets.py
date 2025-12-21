import sys, csv
from collections import defaultdict
thr_ident = 90.0
thr_len = 18
infile, outfile = sys.argv[1], sys.argv[2]
counts = defaultdict(int)
with open(infile) as f:
    r = csv.reader(f, delimiter='\t')
    for row in r:
        qid, pident, length = row[0], float(row[2]), int(row[3])
        if pident >= thr_ident and length >= thr_len:
            counts[qid] += 1
with open(outfile, 'w') as o:
    w = csv.writer(o, delimiter='\t')
    w.writerow(['candidate_id','offtargets_genome'])
    for k,v in counts.items():
        w.writerow([k,v])
