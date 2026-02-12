# !/usr/bin/env python3
"""
energy_compare.py  –  Compare single-point energies (OpenMM)
                       for native vs. relaxed PDB structures
                       that share *one and the same* AMBER prmtop.

If a native PDB is missing atoms the prmtop expects, coordinates are
patched in from the partner relaxed PDB; if they are absent there as
well (typical for stripped HO* hydrogens) a fallback coordinate is
synthesised ∼1 Å along +x from the bonded heavy atom.
"""

from pathlib import Path
from functools import lru_cache
import argparse, sys, time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from openmm import unit, LangevinIntegrator, Vec3   # Vec3 lives here
from openmm.app import (
    Simulation, AmberPrmtopFile, PDBFile,
    HBonds, PME, CutoffNonPeriodic,
)

# --------------------------------------------------------------------- #
# configuration                                                                
# --------------------------------------------------------------------- #

T           = 300*unit.kelvin
GAMMA       = 1/unit.picosecond
DT          = 0.002*unit.picoseconds
NB_CUTOFF   = 1.0*unit.nanometer

ROOT     = Path(__file__).resolve().parent
CSV_OUT  = ROOT/"DESIRED_CSV_NAME"
PARTIAL  = ROOT/"DESIRED_PARTIAL_CSV_NAME"
FIG_OUT  = ROOT/"DESIRED_HISTOGRAM_NAME"

START_FROM = 0            # lowest index you want to process

# ---------------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------------

def build_simulation(prmtop_path: Path) -> Simulation:
    """Return a bare Simulation (no positions yet)."""
    prmtop = AmberPrmtopFile(str(prmtop_path))
    nb_meth = PME if prmtop.topology.getUnitCellDimensions() else CutoffNonPeriodic
    system  = prmtop.createSystem(
        nonbondedMethod   = nb_meth,
        nonbondedCutoff   = NB_CUTOFF,
        constraints       = HBonds,
    )
    return Simulation(prmtop.topology, system,
                      LangevinIntegrator(T, GAMMA, DT))


@lru_cache(maxsize=None)
def get_sim(prmtop_path: Path) -> Simulation:
    """Cache-wrapper so every PRMTOP is parsed just once."""
    return build_simulation(prmtop_path)


def iter_atoms_with_keys(pdb: PDBFile):
    """Yield ((chain, resSeq, iCode, atomName), position) tuples."""
    for at, pos in zip(pdb.topology.atoms(), pdb.positions):
        key = (at.residue.chain.id, at.residue.id, at.residue.insertionCode, at.name)
        yield key, pos


def _vector_along_x(atom):
    """1 Å unit-vector pointing along +x  – used for dummy H coords."""
    return Vec3(0.1, 0.0, 0.0) * unit.nanometer   # 0.1 nm ≈ 1 Å


def patched_positions(topology, pdb_primary: Path, pdb_fallback: Path, verbose=False):
    """
    Return positions list matching *topology* atom order.
    • start from fallback PDB
    • overwrite coords that *are* present in primary
    • fill anything still missing with a reasonable guess
    """
    prim = PDBFile(str(pdb_primary))
    if len(prim.positions) == topology.getNumAtoms():
        return prim.positions            # nothing missing – fast-path

    fall = PDBFile(str(pdb_fallback))
    pos_dict = {k: p for k, p in iter_atoms_with_keys(fall)}

    # overwrite with whatever the primary *does* have
    for k, p in iter_atoms_with_keys(prim):
        pos_dict[k] = p

    filled, n_copied, n_synth = [], 0, 0
    primary_keys = {k for k, _ in iter_atoms_with_keys(prim)}

    for at in topology.atoms():
        k = (at.residue.chain.id, at.residue.id, at.residue.insertionCode, at.name)
        if k in pos_dict:
            filled.append(pos_dict[k])
            if k not in primary_keys:
                n_copied += 1          # came from fallback PDB
        else:
            # synthesize a plausible coordinate (only hydrogens should hit this)
            res_atoms = list(at.residue.atoms())
            ref       = res_atoms[0] if res_atoms else at
            ref_key   = (ref.residue.chain.id, ref.residue.id,
                         ref.residue.insertionCode, ref.name)
            ref_pos = pos_dict.get(ref_key, Vec3(0.0, 0.0, 0.0) * unit.nanometer)
            filled.append(ref_pos + _vector_along_x(at))
            n_synth += 1

    if verbose and (n_copied or n_synth):
        print(f"    patched {n_copied:3d} atoms (copied) + {n_synth:3d} (synthetic)")
    return filled


def single_point_energy(sim: Simulation, positions):
    sim.context.setPositions(positions)
    state = sim.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def pair_iterator(native_dir: Path, relaxed_dir: Path):
    """Yield (stem, native_pdb, relaxed_pdb, prmtop) tuples
       starting from model *START_FROM*."""
    native  = {p.stem.replace('_native',''):  p for p in native_dir .glob('*.pdb')}
    relaxed = {p.stem.replace('_relaxed',''): p for p in relaxed_dir.glob('*.pdb')}

    for stem in sorted(native.keys() & relaxed.keys()):
        # --- NEW FILTER ----------------------------------------------------
        # Extract the numeric part and skip until we reach START_FROM
        import re
        match = re.search(r'\d{6}', stem)
        idx = int(match.group()) if match else None
        if idx is not None and idx < START_FROM:
            continue
        # -------------------------------------------------------------------

        prmtop = relaxed[stem].with_suffix('.prmtop')
        if prmtop.exists():
            yield stem, native[stem], relaxed[stem], prmtop

# ---------------------------------------------------------------------
# main loop
# ---------------------------------------------------------------------

def main(native_dir: Path, relaxed_dir: Path, verbose=False):
    rows = []
    start_time = time.time()

    for stem, nat_pdb, rel_pdb, prmtop in pair_iterator(native_dir, relaxed_dir):
        if verbose:
            print(f"\n  {stem}")

        try:
            sim = get_sim(prmtop)
        except Exception as e:
            print(f"  {stem}: system build failed → {e}")
            continue

        try:
            pos_nat = patched_positions(sim.topology, nat_pdb, rel_pdb, verbose)
            pos_rel = PDBFile(str(rel_pdb)).positions
            e_nat   = single_point_energy(sim, pos_nat)
            e_rel   = single_point_energy(sim, pos_rel)
        except Exception as e:
            print(f"  {stem}: energy calc failed → {e}")
            continue

        rows.append(dict(model=stem, E_native=e_nat,
                         E_relaxed=e_rel, Delta_E=e_rel-e_nat))

        # Write partial results
        pd.DataFrame(rows).sort_values('Delta_E').to_csv(PARTIAL, index=False)

        elapsed = time.time() - start_time
        print(f"  {stem} complete – elapsed: {elapsed:.1f} sec")

        if verbose:
            print(f"    E_native  = {e_nat/1000:10,.1f} ×10³ kJ"
                  f" | E_relaxed = {e_rel/1000:10,.1f}"
                  f" | Δ = {(e_rel-e_nat)/1000:10,.1f}")

    if not rows:
        sys.exit("  No energies computed (all pairs failed).")

    df = pd.DataFrame(rows).sort_values('Delta_E')
    df.to_csv(CSV_OUT, index=False)
    print(f"\n  energies → {CSV_OUT}  ({len(df)} models)")
    print(f"  checkpoint → {PARTIAL}")

    # ------------- histograms -------------------------------------------
    plt.style.use('ggplot')
    fig, ax = plt.subplots(1, 2, figsize=(12,5))

    ax[0].hist(df['E_native']/1000,  bins=30, alpha=.6, label='native')
    ax[0].hist(df['E_relaxed']/1000, bins=30, alpha=.6, label='relaxed')
    ax[0].set_xlabel('Potential energy (10³ kJ mol⁻¹)')
    ax[0].set_title('Absolute energies'); ax[0].legend()

    ax[1].hist(df['Delta_E']/1000, bins=30, color='steelblue')
    ax[1].axvline(0, color='k', lw=1)
    ax[1].set_xlabel('ΔE (relaxed − native) / 10³ kJ mol⁻¹')
    ax[1].set_title('Energy change')

    fig.tight_layout(); fig.savefig(FIG_OUT, dpi=300)
    print(f"  histograms → {FIG_OUT}")

# --------------------------------------------------------------------- #
# CLI                                                                   
# --------------------------------------------------------------------- #

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="OpenMM energy comparison")
    ap.add_argument("--native_dir",  type=Path, required=True,
                    help="directory with *_native.pdb files")
    ap.add_argument("--relaxed_dir", type=Path, required=True,
                    help="directory with *_relaxed.pdb + *.prmtop files")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    if not args.native_dir.is_dir() or not args.relaxed_dir.is_dir():
        sys.exit("  supplied directories do not exist")

    main(args.native_dir, args.relaxed_dir, args.verbose)