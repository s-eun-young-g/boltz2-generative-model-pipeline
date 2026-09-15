#!/usr/bin/env python3
"""
sasa_vs_rl.py
-------------
Solvent-accessible surface area of the predicted structure against ribosomal
load. This is the negative control for the energy result: if relaxed energy
correlates with translation because the structures carry real physical
information, a cruder geometric summary of the same structures should not.

It does not. See the README.

Reads `rl_sasa_pairs.csv` (pdb, total_sasa, rl) and prints the correlation.
To regenerate that file from structures, compute per-PDB total SASA with
freesasa and join to the ribosomal load column of the Sample et al. table.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr

PAIRS = Path(__file__).resolve().parent / "rl_sasa_pairs.csv"


def main() -> None:
    d = pd.read_csv(PAIRS).dropna(subset=["total_sasa", "rl"])
    r, p = pearsonr(d.total_sasa, d.rl)
    rho, p_rho = spearmanr(d.total_sasa, d.rl)
    print(f"n = {len(d):,}")
    print(f"Pearson  r   = {r:+.4f}  (p = {p:.3f})")
    print(f"Spearman rho = {rho:+.4f}  (p = {p_rho:.3f})")
    print(f"r^2          = {r ** 2:.5f}")


if __name__ == "__main__":
    main()
