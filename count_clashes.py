#!/usr/bin/env python3
"""
count_clashes.py
----------------
Count steric clashes in predicted vs relaxed RNA structures, swept over a range
of distance thresholds.

A "clash" here is a pair of NON-BONDED heavy atoms closer than the threshold.
Bonded pairs are excluded from the pair list before distances are computed, so a
covalent bond never counts as a clash.

Run it once over the native (Boltz-2 output) PDBs and once over the relaxed PDBs,
then join on model id to get the native-vs-relaxed table in `clash_summary.csv`.

Usage:
    python count_clashes.py --pdb-dir <dir> --out native_clashes.csv
    python count_clashes.py --pdb-dir <dir> --out relaxed_clashes.csv

Note on the threshold sweep: below about 1.5 Å the count is a clash count. Above
that it starts picking up ordinary 1-3 covalent geometry, which is why the top of
the sweep inverts and is excluded from the headline number. See the README.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd

# Thresholds in Angstrom. mdtraj works in nm, so these are divided by 10 below.
THRESHOLDS_A = np.round(np.arange(0.025, 1.601, 0.025), 3)


def clash_counts(pdb_file: Path, thresholds_a: np.ndarray) -> dict[float, int]:
    """Non-bonded atom pairs closer than each threshold, for one structure."""
    traj = md.load(str(pdb_file))

    bonded = {
        (min(b[0].index, b[1].index), max(b[0].index, b[1].index))
        for b in traj.topology.bonds
    }
    pairs = np.array(
        [p for p in traj.topology.select_pairs("all", "all")
         if (min(p), max(p)) not in bonded]
    )
    if pairs.size == 0:
        return {t: 0 for t in thresholds_a}

    dists_nm = md.compute_distances(traj, pairs)[0]
    return {t: int((dists_nm < t / 10.0).sum()) for t in thresholds_a}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    rows, failed = [], []
    for pdb in sorted(args.pdb_dir.glob("*.pdb")):
        try:
            counts = clash_counts(pdb, THRESHOLDS_A)
            rows.append({"pdb": pdb.name, **{f"{t:.3f}": c for t, c in counts.items()}})
        except Exception as exc:                      # a malformed PDB is logged, never silently dropped
            print(f"FAILED {pdb.name}: {exc}")
            failed.append(pdb.name)

    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(f"wrote {len(rows):,} structures to {args.out}")
    if failed:
        Path(str(args.out) + ".failed.txt").write_text("\n".join(failed))
        print(f"{len(failed)} structures failed, names written beside the output")


if __name__ == "__main__":
    main()
