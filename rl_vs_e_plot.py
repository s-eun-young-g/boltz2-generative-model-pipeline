# !/usr/bin/env python3
"""
rl vs E_relaxed (ATG/AUG-excluded) with IQR outlier removal
-----------------------------------------------------------
1) Map each PDB → (model, canonical sequence) from native PDBs.
2) EXCLUDE sequences containing AUG/ATG (after U->T canonicalisation).
3) Join to UTR CSV by sequence to bring in 'rl'.
4) Join to energy CSV by model to bring in 'E_relaxed'.
5) Compute mean, Q1, Q3, IQR for E_relaxed and remove outliers using the
   Tukey rule: values < Q1-1.5*IQR or > Q3+1.5*IQR.
6) Print & save the excluded model IDs.
7) Plot rl vs E_relaxed (linear), with Pearson r in the title.

Outputs:
- RESULTS_DIR/merged_no_ATG.csv                    (pre-outlier table)
- RESULTS_DIR/excluded_outliers_by_IQR.csv         (model, E_relaxed)
- RESULTS_DIR/merged_no_ATG_iqr_filtered.csv       (post-outlier table)
- RESULTS_DIR/rl_vs_E_relaxed_linear_no_ATG.png    (plot)
"""

from pathlib import Path
from multiprocessing import Pool, cpu_count
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from scipy.stats import pearsonr

# ---------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------
ROOT             = Path(__file__).resolve().parent
UTR_CSV_PATH     = ROOT / "../Therascript_master/data/GSE114002/GSM3130435_egfp_unmod_1.csv"
ENERGY_CSV_PATH  = ROOT / "PATH_TO_ENERGY_CSV"
PDBS_PATH        = ROOT / "PATH_TO_NATIVE_PDBS"
RESULTS_DIR      = ROOT / "DESIRED_PLOT_LOCATION"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
# CSV columns
UTR_SEQ_COL        = "utr"
UTR_RL_COL         = "rl"
ENERGY_KEY_COL     = "model"
ENERGY_RELAXED_COL = "E_relaxed"

# Optional exclusions (by model id)
EXCLUDE_MODELS = None

# Optional energy filters (applied BEFORE IQR; set None to disable)
RELAXED_MAX      = None   # keep only x < RELAXED_MAX
RELAXED_MIN_ABS  = None   # keep only |x| > RELAXED_MIN_ABS

# 3-letter → 1-letter mapping (nucleotides)
_3_TO_1 = {
    "A": "A", "ADE": "A", "DA": "A",
    "C": "C", "CYT": "C", "DC": "C",
    "G": "G", "GUA": "G", "DG": "G",
    "U": "U", "URA": "U", "DT": "U",
}

# ---------------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------------
def canonical(seq: str) -> str:
    """RNA → DNA (U→T) and uppercase."""
    return str(seq).upper().replace("U", "T")

def pdb_to_seq(pdb_path: Path) -> tuple[str, str]:
    """Return (model_id, canonical DNA sequence) for one PDB."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(pdb_path.stem, pdb_path)
    letters = []
    for model in struct:        # first model only
        for chain in model:
            for res in chain:
                r = res.get_resname().strip()
                if r in _3_TO_1:
                    letters.append(_3_TO_1[r])
        break
    return pdb_path.stem, canonical("".join(letters))

def extract_sequences(pdb_dir: Path, threads=None) -> pd.DataFrame:
    files = sorted(pdb_dir.glob("*.pdb"))
    if not files:
        sys.exit(f"No PDBs found in {pdb_dir}")
    threads = threads or max(1, cpu_count() - 1)
    with Pool(threads) as pool:
        data = pool.map(pdb_to_seq, files)
    df = pd.DataFrame(data, columns=[ENERGY_KEY_COL, "sequence"])
    print(f" extracted sequences for {len(df):,} PDBs")
    return df

def linear_scatter(df: pd.DataFrame, x: str, y: str, out_png: Path):
    """Scatter with Pearson r shown in the title."""
    df = df.dropna(subset=[x, y]).copy()
    df[x] = pd.to_numeric(df[x], errors="coerce")
    df[y] = pd.to_numeric(df[y], errors="coerce")
    df = df.dropna(subset=[x, y])
    df = df[np.isfinite(df[x]) & np.isfinite(df[y])]
    if df.empty:
        print(f"[WARN] No data to plot for {out_png.name}")
        return
    r, p = pearsonr(df[x], df[y])
    title = f"{y} vs {x} (linear)\nPearson r = {r:.3f} (p = {p:.1e})"
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df[x], df[y], s=16, alpha=0.4)
    ax.set_xlabel(x, fontsize=12); ax.set_ylabel(y, fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.tick_params(axis='both', labelsize=10)
    ax.grid(True, which="both", lw=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close()
    print(f" Saved: {out_png} ({len(df):,} points)")

# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------
def main():
    # 1) PDB → (model, sequence)
    seq_df = extract_sequences(PDBS_PATH)

    # 2) EXCLUDE any sequence containing AUG/ATG (after canonicalisation, AUG→ATG)
    before = len(seq_df)
    seq_df = seq_df[~seq_df["sequence"].str.contains("ATG")]
    print(f" kept {len(seq_df):,} / {before:,} PDBs WITHOUT ATG (AUG excluded via U→T)")

    # 3) Load UTR table and canonicalize sequences the same way
    utr_df = pd.read_csv(UTR_CSV_PATH).rename(columns=str.strip)
    utr_df[UTR_SEQ_COL] = utr_df[UTR_SEQ_COL].map(canonical)

    # 4) Join by sequence to bring in 'rl'
    merged = seq_df.merge(
        utr_df[[UTR_SEQ_COL, UTR_RL_COL]],
        left_on="sequence", right_on=UTR_SEQ_COL,
        how="left", validate="many_to_one"
    )

    # 5) Load energy table (minimized energies) and coerce numeric
    energy_df = pd.read_csv(ENERGY_CSV_PATH).rename(columns=str.strip)

    # Normalize key column
    energy_df[ENERGY_KEY_COL] = energy_df[ENERGY_KEY_COL].astype(str).str.strip()
    energy_df[ENERGY_RELAXED_COL] = pd.to_numeric(energy_df[ENERGY_RELAXED_COL], errors="coerce")

    # (Optional) quick diagnostics
    dups = energy_df[ENERGY_KEY_COL].duplicated(keep=False).sum()
    if dups:
        print(f"  energy CSV has {dups} duplicate model rows; collapsing to one per model.")

    # Strategy: keep the last non-null record per model (change to 'first' or use an aggregate if you prefer)
    energy_df = (
        energy_df
        .sort_values([ENERGY_KEY_COL])        # stable order if file is chronological
        .drop_duplicates(subset=[ENERGY_KEY_COL], keep="last")
    )

    # Ensure left table also has unique models (defensive)
    merged = merged.drop_duplicates(subset=[ENERGY_KEY_COL], keep="first")  

    # 6) Join by model to bring in E_relaxed
    merged = (
        merged
        .merge(energy_df[[ENERGY_KEY_COL, ENERGY_RELAXED_COL]],
               on=ENERGY_KEY_COL, how="left", validate="one_to_one")
    )

    # Drop known outliers by model, if any (explicit list)
    if EXCLUDE_MODELS:
        merged = merged[~merged[ENERGY_KEY_COL].isin(EXCLUDE_MODELS)]

    # Save pre-outlier table
    pre_csv = RESULTS_DIR / "merged_no_ATG.csv"
    merged.to_csv(pre_csv, index=False)
    print(f" merged_no_ATG.csv written ({len(merged):,} rows)")

    # IQR outlier detection on E_relaxed 
    work = merged.dropna(subset=[ENERGY_RELAXED_COL]).copy()
    # Apply optional simple filters BEFORE computing IQR (set to None to skip)
    if RELAXED_MIN_ABS is not None:
        work = work[work[ENERGY_RELAXED_COL].abs() > RELAXED_MIN_ABS]
    if RELAXED_MAX is not None:
        work = work[work[ENERGY_RELAXED_COL] < RELAXED_MAX]

    if work.empty:
        sys.exit(" No rows with numeric E_relaxed to compute IQR.")

    x = work[ENERGY_RELAXED_COL]
    mean_val = x.mean()
    q1, q3 = x.quantile(0.25), x.quantile(0.75)
    iqr = q3 - q1
    lower, upper = q1 - 1.5 * iqr, q3 + 1.5 * iqr

    print("\n Energy stats (E_relaxed):")
    print(f"mean = {mean_val:.6g}")
    print(f"Q1   = {q1:.6g}")
    print(f"Q3   = {q3:.6g}")
    print(f"IQR  = {iqr:.6g}")
    print(f"Lower fence = {lower:.6g}")
    print(f"Upper fence = {upper:.6g}")

    mask_inliers = (work[ENERGY_RELAXED_COL] >= lower) & (work[ENERGY_RELAXED_COL] <= upper)
    outliers_df = work.loc[~mask_inliers, [ENERGY_KEY_COL, ENERGY_RELAXED_COL]].copy()
    inliers_df  = work.loc[mask_inliers].copy()

    # Report & save excluded outliers
    excl_path = RESULTS_DIR / "excluded_outliers_by_IQR.csv"
    outliers_df.to_csv(excl_path, index=False)
    print(f"\n Excluded {len(outliers_df):,} outlier model(s) by IQR.")
    if not outliers_df.empty:
        print("Models excluded:")
        for m in outliers_df[ENERGY_KEY_COL].tolist():
            print("  -", m)
        print(f" Saved list to: {excl_path}")

    # Build final plotting table: keep inliers and preserve rl column
    filtered = merged.merge(
        inliers_df[[ENERGY_KEY_COL]],
        on=ENERGY_KEY_COL, how="inner"
    )

    post_csv = RESULTS_DIR / "merged_no_ATG_iqr_filtered.csv"
    filtered.to_csv(post_csv, index=False)
    print(f"\n merged_no_ATG_iqr_filtered.csv written ({len(filtered):,} rows)")

    # Plot rl vs E_relaxed (linear) using IQR-filtered data 
    out_png = RESULTS_DIR / f"{UTR_RL_COL}_vs_{ENERGY_RELAXED_COL}_vac_linear_ATG_outlier_filtered_12193.png"
    linear_scatter(filtered, x=ENERGY_RELAXED_COL, y=UTR_RL_COL, out_png=out_png)

if __name__ == "__main__":
    main()
