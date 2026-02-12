# !/usr/bin/env python
"""
run_predictions.py
------------------------
Batch-runs Boltz on 5'-UTR sequences read from a CSV.

Usage:
    python run_batch_predictions.py        # CPU only
    python run_batch_predictions.py --gpu  # try GPU, fall back to CPU if unavailable
"""

from __future__ import annotations
import argparse, os, shutil, subprocess, sys, time, yaml
import pandas as pd
from pathlib import Path
from datetime import datetime

# ---------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------
CSV_PATH      = Path("../Therascript_master/data/GSE114002/GSM3130435_egfp_unmod_1.csv")
YAML_ROOT     = Path("PATH_TO_YAML_INPUTS_FOLDER")
PRED_ROOT     = Path("PATH_TO_DESTINATION_OF_OUTPUTS")
BATCH_SIZE    = 100
SAMPLING_STEPS   = 50
RECYCLING_STEPS  = 1
OUTPUT_FORMAT    = "pdb"
MAX_SEQ_LEN      = 50           # trim UTRs to 50 nt



# ---------------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------------
def make_yaml(seq:str, key:str, out_dir:Path) -> Path:
    """Write a single-sequence YAML in <out_dir>/<key>.yaml ."""
    out_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = out_dir / f"{key}.yaml"
    content = {"sequences":[{"rna":{"id":"A","sequence":seq,"msa":"empty"}}]}
    with open(yaml_path, "w") as fh:
        yaml.dump(content, fh, default_flow_style=False)
    return yaml_path

def gpu_is_available() -> bool:
    """Return True iff PyTorch detects a CUDA device _and_ nvidia-smi is present"""
    try:
        import torch, subprocess, shutil
        has_torch_gpu = torch.cuda.is_available()
        has_driver    = shutil.which("nvidia-smi") is not None
        return has_torch_gpu and has_driver
    except Exception:
        return False

def run_boltz(yaml_dir:Path, out_dir:Path, use_gpu:bool) -> None:
    """Call boltz predict on a directory of YAML files."""
    accel = "gpu" if use_gpu else "cpu"
    cmd   = [
        "boltz", "predict", str(yaml_dir),
        "--out_dir", str(out_dir),
        "--accelerator", accel,
        "--sampling_steps", str(SAMPLING_STEPS),
        "--recycling_steps", str(RECYCLING_STEPS),
        "--output_format", OUTPUT_FORMAT,
        # "--override",                 # re-run even if some outputs exist
    ]
    if use_gpu:
        cmd += ["--devices", "1"]     # tell Lightning to use 1 GPU

    print(f"boltz predict  |  {yaml_dir.name}  →  {out_dir.name}  |  accel={accel}")
    try:
        subprocess.run(cmd, check=True)
        print("finished\n")
    except subprocess.CalledProcessError as e:
        print(f"boltz failed (return-code {e.returncode})")
        print("command:", " ".join(cmd))
        print("stderr output above.\n")

def handle_batch(batch_paths:list[Path], batch_no:int,
                 yaml_root:Path, pred_root:Path, use_gpu:bool)->None:
    # move listed YAMLs to Batch_<n>/ and run Boltz there 
    batch_dir = yaml_root / f"Batch_{batch_no:03d}"
    out_dir   = pred_root / f"Batch_{batch_no:03d}"
    batch_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    # move yaml files into the batch directory
    for p in batch_paths:
        p.rename(batch_dir / p.name)

    run_boltz(batch_dir, out_dir, use_gpu)

# ---------------------------------------------------------------------
# main loop
# ---------------------------------------------------------------------
def main(argv:list[str]|None=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gpu", action="store_true",
                        help="try to run on GPU; falls back to CPU if CUDA unavailable")
    args = parser.parse_args(argv)

    # data prep
    print("Loading CSV …")
    df = pd.read_csv(CSV_PATH).sort_values("total_reads", ascending=False).reset_index(drop=True)
    df_subset = df.loc[:19_999, "utr"]      # top 20,000 sequences
    print("Sequences loaded:", len(df_subset))

    # REPLACE T with U (only in sequence strings)
    df_subset = df_subset.str.replace('T', 'U', regex=False)
    print("Replaced all T’s with U’s in sequences.")

    # decide accelerator
    want_gpu = args.gpu
    have_gpu = gpu_is_available()
    use_gpu  = want_gpu and have_gpu
    if want_gpu and not have_gpu:
        print("--gpu requested but no usable CUDA device found, using CPU instead.\n")

    # timestamped root so repeated runs don’t overwrite
    run_stamp   = datetime.now().strftime("%Y%m%d-%H%M%S")
    yaml_root   = YAML_ROOT / f"run_{run_stamp}"
    pred_root   = PRED_ROOT / f"run_{run_stamp}"
    yaml_root.mkdir(parents=True, exist_ok=True)
    pred_root.mkdir(parents=True, exist_ok=True)

    start = time.time()
    batch_no, batch_yaml = 1, []

    for idx, seq in enumerate(df_subset, start=1):
        seq = seq[:MAX_SEQ_LEN]          # (optional) truncate
        key = f"seq_{idx:06d}"
        batch_yaml.append(make_yaml(seq, key, yaml_root))

        if len(batch_yaml) == BATCH_SIZE:
            handle_batch(batch_yaml, batch_no, yaml_root, pred_root, use_gpu)
            batch_no  += 1
            batch_yaml = []

    # leftover YAMLs
    if batch_yaml:
        handle_batch(batch_yaml, batch_no, yaml_root, pred_root, use_gpu)

    print(f"all batches finished in {time.time()-start:,.1f}s → {pred_root}")


# entry-point
if __name__ == "__main__":
    main()
