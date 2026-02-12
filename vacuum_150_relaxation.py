# !/usr/bin/env python3
# ------------------------------------------------------------
# RNA relaxation pipeline (tLeap → OpenMM)
# SAFE version: skips completed, logs failures
# ------------------------------------------------------------
import subprocess
import tempfile
import warnings
from pathlib import Path
from time import time

import mdtraj as md
import pandas as pd
from openmm.app import *
from openmm import *
from openmm.unit import *

warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------
ROOT         = Path(__file__).resolve().parent
PDB_DIR      = ROOT / "PATH_TO_EDITED_PDBS"
RELAXED_DIR  = ROOT / "PATH_TO_DESTINATION_FOLDER"
RELAXED_DIR.mkdir(parents=True, exist_ok=True)
CSV_PATH     = RELAXED_DIR / "DESIRED_CSV_NAME"
FAILED_PATH  = RELAXED_DIR / "failed_pdbs.txt"


# ---------------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------------
def run_tleap(input_pdb: str, prmtop_path: str, inpcrd_path: str):
    """Build Amber topology/coordinates with tLeap."""
    tleap_input = f"""\
source leaprc.RNA.OL3
mol = loadpdb {input_pdb}
# Uncomment the next line if you want a periodic box
# setbox mol centers
saveamberparm mol {prmtop_path} {inpcrd_path}
quit
"""
    with tempfile.NamedTemporaryFile("w", suffix=".in", delete=False) as tf:
        tf.write(tleap_input)
        script_path = tf.name

    result = subprocess.run(
        ["tleap", "-f", script_path],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print("tLeap failed:\n", result.stdout, result.stderr)
        raise RuntimeError("tLeap execution failed.")
    else:
        print(f"tLeap succeeded for {input_pdb}")


def _method_to_string(method_const):
    """Return a readable name for an OpenMM non‑bonded method."""
    mapping = {
        NoCutoff:           "NoCutoff",
        CutoffNonPeriodic:  "CutoffNonPeriodic",
        CutoffPeriodic:     "CutoffPeriodic",
        Ewald:              "Ewald",
        PME:                "PME",
        LJPME:              "LJPME",
    }
    return mapping.get(method_const, str(method_const))


def minimize_with_openmm(prmtop_path: str, inpcrd_path: str,
                         output_pdb: str,
                         tolerance=10*kilojoules_per_mole/nanometer,
                         max_iters=150):
    """Minimise Amber system with OpenMM and write relaxed PDB."""
    prmtop = AmberPrmtopFile(prmtop_path)
    inpcrd = AmberInpcrdFile(inpcrd_path)

    # Choose method based on periodic box presence
    if prmtop.topology.getUnitCellDimensions() is None:
        nb_method = CutoffNonPeriodic
        nb_cutoff = 1.0 * nanometer
    else:
        nb_method = PME
        nb_cutoff = 1.0 * nanometer

    system = prmtop.createSystem(
        nonbondedMethod=nb_method,
        nonbondedCutoff=nb_cutoff,
        constraints=HBonds
    )

    integrator = LangevinIntegrator(300*kelvin, 1/picosecond,
                                    0.002*picoseconds)
    sim = Simulation(prmtop.topology, system, integrator)
    sim.context.setPositions(inpcrd.positions)

    print(f"  • minimising with {_method_to_string(nb_method)} …")
    sim.minimizeEnergy(tolerance, max_iters)

    with open(output_pdb, "w") as fh:
        PDBFile.writeFile(
            prmtop.topology,
            sim.context.getState(getPositions=True).getPositions(),
            fh
        )
    print(f"  • minimised PDB written to {output_pdb}")


def clashes_count(pdb_file, thresh=0.7):
    """Count close‑contact pairs (clashes) below `thresh` nm."""
    traj  = md.load(pdb_file)
    pairs = traj.topology.select_pairs("all", "all")
    dists = md.compute_distances(traj, pairs)
    return int((dists[0] < thresh).sum())


# ---------------------------------------------------------------------
# main loop
# ---------------------------------------------------------------------
# Load completed PDBs
done_pdbs = set(p.name for p in RELAXED_DIR.glob("*.pdb"))
print(f"Found {len(done_pdbs)} already relaxed PDBs in {RELAXED_DIR}")

stats = []

for pdb_file in sorted(PDB_DIR.glob("*.pdb")):
    if pdb_file.name in done_pdbs:
        print(f"Already relaxed: {pdb_file.name}, skipping.")
        continue

    print(f"\nProcessing {pdb_file.name}")
    prmtop_file = RELAXED_DIR / pdb_file.with_suffix(".prmtop").name
    inpcrd_file = RELAXED_DIR / pdb_file.with_suffix(".inpcrd").name
    relaxed_pdb = RELAXED_DIR / pdb_file.name

    try:
        # 1 – tLeap build
        run_tleap(str(pdb_file), str(prmtop_file), str(inpcrd_file))

        # 2 – OpenMM minimisation
        t0 = time()
        minimize_with_openmm(str(prmtop_file), str(inpcrd_file), str(relaxed_pdb))
        elapsed = time() - t0

        # 3 – clash stats
        before = clashes_count(str(pdb_file))
        after  = clashes_count(str(relaxed_pdb))

        result = dict(
            pdb_name=pdb_file.name,
            time_seconds=elapsed,
            clashes_before=before,
            clashes_after=after,
            clash_difference=before - after
        )
        stats.append(result)

        # Append to CSV immediately
        pd.DataFrame([result]).to_csv(CSV_PATH, mode='a',
                                      header=not CSV_PATH.exists(),
                                      index=False)

    except Exception as e:
        print(f"Error processing {pdb_file.name}: {e}")
        with open(FAILED_PATH, "a") as fail_log:
            fail_log.write(f"{pdb_file.name}\n")
        continue

# ---------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------
print(f"\n Results table saved to {CSV_PATH}")
if stats:
    print("Average minimisation time  : %.2f s" % pd.DataFrame(stats)["time_seconds"].mean())
    print("Average clash reduction    : %.2f"   % pd.DataFrame(stats)["clash_difference"].mean())
else:
    print(" No new PDBs processed in this run.")
