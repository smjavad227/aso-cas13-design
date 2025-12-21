import sys, csv

inp, offt, outp = sys.argv[1], sys.argv[2], sys.argv[3]

# Load off-target counts
offt_map = {}
with open(offt) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        offt_map[row['candidate_id']] = row['offtargets_genome']

def make_id(row):
    return f"{row['type']}|{row['transcript_id']}|{row['win_start']}|{row['win_end']}"

# Read input file explicitly
with open(inp) as f:
    header = f.readline().strip().split('\t')   # read header line
    # ✅ Add the new column if not present
    if 'offtargets_genome' not in header:
        header.append('offtargets_genome')

    rows = []
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == len(header) - 1:  # minus the new column
            row = dict(zip(header[:-1], parts))
            rows.append(row)

# Write output with updated counts
with open(outp, 'w', newline='') as o:
    writer = csv.DictWriter(o, fieldnames=header, delimiter='\t')
    writer.writeheader()
    for row in rows:
        key = make_id(row)
        if key in offt_map:
            row['offtargets_genome'] = offt_map[key]
        else:
            row['offtargets_genome'] = "0"  # default if not found
        writer.writerow(row)
