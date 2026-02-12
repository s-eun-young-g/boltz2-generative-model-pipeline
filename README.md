# 1. Environment/dependencies installed

I established my environment and installed dependencies as follows:

make env
conda create -y -n boltz-env python=3.10

activate
conda activate boltz-env

core scientific stack
conda install -y -c conda-forge numpy pandas scipy matplotlib

structure I/O
conda install -y -c conda-forge biopython pdbfixer

MD engine
conda install -y -c conda-forge openmm

AMBER tools (tleap, parmchk2, etc.)
conda install -y -c conda-forge ambertools

See the full list of packages installed in my environment in packages.csv.

# 2. Dataset used

I used the egfp_unmod_1 dataset found here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114002

# 3. Running Boltz-2 predictions

I ran my predictions using the most up-to-date version of Boltz-2 as of June 30th, 2025. To run the predictions, I used the script uploaded as run_predictions.py. 

Some notes: this script finds and replaces all instances of the nucleotide T with the nucleotide U to avoid errors in Boltz’ processing of the sequences. Furthermore, a 50 nucleotide cap is included for both performance and runtime reasons, but more importantly because it is around the maximum sequence length for which Boltz remains reliable and accurate. I included this cap to avoid potential misfolding of any longer sequences.

# 4. Running energy minimizations

Before I was able to run minimizations, I needed to modify the pdb files output by Boltz first. I deleted the first 3 atoms of every PDB file, as tleap was unable to recognize the 5′-phosphate cap that remained on the first residues (as the force field used, OL3.RNA, only has templates for de-phosphorylated 5’ ends).

I then cleaned up the formatting of all the pdb files, as I encountered spacing issues.

I then used the script vacuum_150_relaxation.py to relax the edited pdbs in a vacuum.

# 5. Preparing files for energy calculations

Now, I had two sets of pdbs: the edited native pdbs and the relaxed versions of these pdbs. Before running energy calculations, I needed to prepare these files using two scripts. The first script, fix_native_pdbs.sh, repairs all of the edited native pdbs I removed the first 3 atoms of (repairing the minimized pdbs is not necessary, as it was done as part of the relaxation process). The second script, make_prmtops.sh, generate prmtops/rst7 files for all the repaired native and relaxed pdbs.

# 6. Running energy calculations

I used the script compare_energies.py to compute the energies of the native and relaxed pdbs, as well as their difference.

# 7. Plotting results

I then used rl_vs_e_plot.py to plot the ribosomal load of each sequence to its relaxed energy. This script filters out sequences containing AUG/ATG and filters out outliers according to the IQR method.

# 8. Results

You can find a full list of all ribosome load-energy pairs in rl_E_relaxed_pairs.csv and my final plot, rl_vs_E_relaxed_ATG_outlier_filtered.png, uploaded in this directory.
