# ASO-Cas13 Design Pipeline

This repository contains a reproducible Python pipeline for designing and evaluating antisense oligonucleotides (ASOs) and Cas13 candidates.  
The workflow is modular, step-by-step, and intended for open science and reproducibility.

## Repository Structure
requirements.txt # Python dependencies 
 README.md # Project documentati
 
## Requirements
- Python 3.9+
- Install dependencies:
  ```bash
  pip install -r requirements.txt
External tools required (install separately):

BLAST+

BEDTools

ViennaRNA (RNAfold)

How to Run
Run the scripts in order:
python scripts/00_prepare_inputs.py
python scripts/01_design_windows.py
python scripts/02_filter_candidates.py
...
python scripts/13_visualize_integrated.py
Utility scripts:
python scripts/check_data_files.py
python scripts/check_inputs_reference.py
Outputs
Candidate tables (ASO and Cas13 designs)

Summary statistics

Visualizations (plots, graphs)

Supplementary tables for manuscript use

License
This project is released under the MIT License for open scientific use.

Citation
If you use this pipeline in your research, please cite:

Hashemi, S.M.J. ASO-Cas13 Design Pipeline (2025). GitHub repository:
https://github.com/smjavad227/aso-cas13-design