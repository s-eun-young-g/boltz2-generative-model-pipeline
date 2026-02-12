#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob               # ignore empty globs

for pdb in ~/PATH_TO_NATIVE_PDBS/*.pdb \
           ~/PATH_TO_RELAXED_PDBS/*.pdb; do

    base=${pdb%.pdb}            
    echo "tleap on $(basename "$pdb")"

    tleap -s -f - <<EOF
source leaprc.protein.ff14SB
source leaprc.RNA.OL3 

addAtomTypes { { "H5T" "H" " " } { "H3T" "H" " " } }

mol = loadPDB "$pdb"
setBox mol vdw 10.0
saveAmberParm mol "${base}.prmtop" "${base}.rst7"
quit
EOF
done
echo “All prmtops/rst7 files done."
