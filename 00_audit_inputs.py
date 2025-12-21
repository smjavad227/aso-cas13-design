import os, sys, csv, json, re
from pathlib import Path

ROOT = Path(".")
TARGETS = [
    # Reference
    "outputs/reference/DMPK_ref_GRCh38.fasta",
    # Pre-integration design
    "outputs/design/ASO_candidates_final_grch38.tsv",
    "outputs/design/Cas13_guides_final_grch38.tsv",
    # Merged with conservation
    "outputs/design/ASO_candidates_with_conservation.tsv",
    "outputs/design/Cas13_guides_with_conservation.tsv",
    # Final integrated
    "outputs/design/ASO_candidates_final_integrated.tsv",
    "outputs/design/Cas13_guides_final_integrated.tsv",
    # Stringent sets
    "outputs/design/ASO_candidates_final_stringent.tsv",
    "outputs/design/Cas13_guides_final_stringent.tsv",
    # RNAfold raw outputs
    "outputs/results/ASO_rnafold.txt",
    "outputs/results/Cas13_rnafold.txt",
    # Figures
    "outputs/figures/integrated_score_vs_offtargets.png",
    "outputs/figures/dG_distribution.png",
    "outputs/figures/accessibility_distribution.png",
    "outputs/figures/top_candidates_integrated_barplot.png",
]

REQ_COLS_DESIGN = [
    "type","transcript_id","chr","strand","win_start","win_end",
    "gc","tm","score_base","offtargets_genome","score_final","window_seq"
]
REQ_COLS_MERGED_EXTRAS = [
    "rnafold_dG","has_homopolymer_5+","simple_repeat_score",
    "accessibility_score","isoform_conservation"
]
REQ_COLS_FINAL = REQ_COLS_DESIGN + REQ_COLS_MERGED_EXTRAS + ["final_integrated_score"]

def size_lines(p: Path):
    if not p.exists():
        return {"exists": False, "bytes": 0, "lines": 0}
    try:
        b = p.stat().st_size
    except Exception:
        b = 0
    lines = 0
    try:
        with p.open("rb") as f:
            for _ in f:
                lines += 1
    except Exception:
        lines = 0
    return {"exists": True, "bytes": b, "lines": lines}

def audit_fasta(p: Path):
    res = {"type": "FASTA", "headers": 0, "non_iupac": 0}
    if not p.exists() or p.stat().st_size == 0:
        return res
    headers = 0
    non_iupac = 0
    with p.open() as f:
        for line in f:
            if line.startswith(">"):
                headers += 1
            else:
                non_iupac += len(re.sub(r"[ACGTNacgtn]", "", line.strip()))
    res["headers"] = headers
    res["non_iupac"] = non_iupac
    try:
        with p.open() as f:
            res["preview"] = "".join([next(f) for _ in range(5)])
    except Exception:
        res["preview"] = ""
    return res

def audit_tsv(p: Path, expected_cols=None, numeric_checks=None):
    res = {"type": "TSV", "columns": [], "missing_cols": [], "numeric_invalid": {}}
    if not p.exists() or p.stat().st_size == 0:
        return res
    with p.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        res["columns"] = r.fieldnames or []
        if expected_cols:
            missing = [c for c in expected_cols if c not in res["columns"]]
            res["missing_cols"] = missing
        # Preview first 3 records
        rows = []
        for i, row in enumerate(r):
            rows.append(row)
            if i >= 2: break
        res["preview_rows"] = rows
    # Numeric checks
    if numeric_checks:
        import pandas as pd
        try:
            df = pd.read_csv(p, sep="\t")
            invalid = {}
            for col, kind in numeric_checks.items():
                if col in df.columns:
                    ser = pd.to_numeric(df[col], errors="coerce")
                    invalid[col] = int(ser.isna().sum())
            res["numeric_invalid"] = invalid
        except Exception as e:
            res["numeric_invalid_error"] = str(e)
    return res

def audit_text(p: Path, grep_keywords=None):
    res = {"type": "TEXT", "contains": {}}
    if not p.exists() or p.stat().st_size == 0:
        return res
    try:
        text = p.read_text(errors="ignore")
    except Exception:
        text = ""
    if grep_keywords:
        for k in grep_keywords:
            res["contains"][k] = (k in text)
    # short preview
    try:
        with p.open() as f:
            res["preview"] = "".join([next(f) for _ in range(5)])
    except Exception:
        res["preview"] = ""
    return res

def audit_image(p: Path):
    # Simple size check; optional PNG signature
    res = {"type": "IMAGE", "png_sig_ok": None}
    if not p.exists() or p.stat().st_size == 0:
        return res
    try:
        with p.open("rb") as f:
            sig = f.read(8)
        res["png_sig_ok"] = (sig == b"\x89PNG\r\n\x1a\n")
    except Exception:
        res["png_sig_ok"] = None
    return res

def main():
    report = {}
    for t in TARGETS:
        p = Path(t)
        meta = size_lines(p)
        entry = {"path": str(p), **meta}

        if t.endswith(".fasta"):
            entry["detail"] = audit_fasta(p)

        elif t.endswith(".tsv"):
            expected = None
            numeric = None
            if "final_grch38.tsv" in t:
                expected = REQ_COLS_DESIGN
                numeric = {
                    "gc":"float","tm":"float","score_base":"float","score_final":"float",
                    "win_start":"int","win_end":"int","offtargets_genome":"int"
                }
            if "with_conservation.tsv" in t:
                expected = REQ_COLS_DESIGN + REQ_COLS_MERGED_EXTRAS
                numeric = {
                    "gc":"float","tm":"float","score_base":"float","score_final":"float",
                    "win_start":"int","win_end":"int","offtargets_genome":"int",
                    "rnafold_dG":"float","simple_repeat_score":"float","accessibility_score":"float",
                    "isoform_conservation":"int"
                }
            if "final_integrated.tsv" in t:
                expected = REQ_COLS_FINAL
                numeric = {
                    "gc":"float","tm":"float","score_base":"float","score_final":"float",
                    "win_start":"int","win_end":"int","offtargets_genome":"int",
                    "rnafold_dG":"float","simple_repeat_score":"float","accessibility_score":"float",
                    "isoform_conservation":"int","final_integrated_score":"float"
                }
            if "stringent.tsv" in t:
                expected = REQ_COLS_FINAL
                numeric = {
                    "offtargets_genome":"int","final_integrated_score":"float"
                }
            entry["detail"] = audit_tsv(p, expected, numeric)

        elif t.endswith(".txt"):
            entry["detail"] = audit_text(p, grep_keywords=["MFE","dG","ΔG"])

        elif t.endswith(".png"):
            entry["detail"] = audit_image(p)

        report[t] = entry

    # Also list any unexpected empties in key directories
    extra = {}
    for d in ["outputs/reference","outputs/design","outputs/results","outputs/figures"]:
        dp = Path(d)
        if dp.exists():
            empties = []
            for f in dp.iterdir():
                try:
                    if f.is_file() and f.stat().st_size == 0:
                        empties.append(str(f))
                except Exception:
                    pass
            extra[d] = {"empty_files": empties}
        else:
            extra[d] = {"missing_dir": True}
    report["_directory_empty_scan"] = extra

    print(json.dumps(report, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
