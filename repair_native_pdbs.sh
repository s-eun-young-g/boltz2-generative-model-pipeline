#!/usr/bin/env bash
set -euo pipefail

N_THREADS=8                            # <= your 8 vCPUs
PDB_DIR="$HOME/PATH_TO_NATIVE_PDBS"

export PYTHONUNBUFFERED=1

fix_one() {
    local pdb="$1"
    # NOTE: filename argument ($pdb) comes *before* the heredoc
    python3 - "$pdb" <<'PY'
import sys, pathlib, os               
from pdbfixer import PDBFixer
from openmm.app import PDBFile

pdb_path = pathlib.Path(sys.argv[1])
print(f"[PID {os.getpid()}] Fixing {pdb_path.name}", flush=True)

fixer = PDBFixer(filename=str(pdb_path))
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.0)

with open(pdb_path, "w") as fh:
    PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
PY
}

export -f fix_one

find "$PDB_DIR" -name '*.pdb' -print0 | \
  xargs -0 -n1 -P "$N_THREADS" bash -c 'fix_one "$0"'

echo “All PDBs repaired"
