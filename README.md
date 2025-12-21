# ASO-Cas13 Design Pipeline

This repository contains a reproducible Python pipeline for designing and evaluating antisense oligonucleotides (ASOs) and Cas13 candidates.  
The workflow is modular, step-by-step, and intended for open science and reproducibility.

## Quick Start
Clone the repository and install dependencies:

```bash
git clone https://github.com/smjavad227/aso-cas13-design.git
cd aso-cas13-design
pip install -r requirements.txt
 ```
## External tools required (install separately):
- BLAST+
- BEDTools
- ViennaRNA (RNAfold)
## How to Run
Run the scripts in order:
python scripts/00_audit_inputs.py
python scripts/01_locus_windows.py
python scripts/02_filter_candidates.py
python scripts/03_export_candidates_fasta.py
python scripts/03_parse_blast.py
python scripts/04_parse_blast_genome.py
python scripts/04_summarize_offtargets.py
python scripts/05_merge_offtargets.py
python scripts/05_visualize_candidates.py
python scripts/07_rnafold_energy.py
python scripts/08_repeat_filters.py
python scripts/09_accessibility_proxy.py
python scripts/10_isoform_conservation.py
python scripts/11_final_integration.py
python scripts/12_stringent_final_sets.py
python scripts/13_visualize_integrated.py
## Utility scripts:
python scripts/check_data_files.py
python scripts/check_inputs_reference.py
## Outputs
Candidate tables (ASO and Cas13 designs)
Summary statistics
Visualizations (plots, graphs)
Supplementary tables for manuscript use
## License
This project is released under the MIT License for open scientific use.
## Citation
If you use this pipeline in your research, please cite:
Hashemi, S.M.J. ASO-Cas13 Design Pipeline (2025). GitHub repository: https://github.com/smjavad227/aso-cas13-design
