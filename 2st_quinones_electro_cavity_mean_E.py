# -*- coding: utf-8 -*-
"""
Quinone redox: energies + error bars + 3D overlay plots
- Compares 5 quinone couples at B3LYP-D3BJ/def2-SVP + PCM(MeCN)
- Produces Excel/CSV results and PNG figures.
"""

import os
import psi4
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from openpyxl import Workbook

# ---------------------------
# Global settings
# ---------------------------
FARADAY = 96485.33212       # C/mol
METHOD = "B3LYP-D3BJ"
BASIS = "def2-SVP"

# PCM template with escaped braces for str.format()
PCM_TEMPLATE = """
Units = Angstrom
Medium {{
    SolverType = IEFPCM
    Solvent = Acetonitrile
}}
Cavity {{
    RadiiSet = UFF
    Type = GePol
    Area = {area}
    Mode = Implicit
}}
"""

# Psi4 runtime settings
psi4.set_memory('12 GB')      # keep headroom on an 8 GB laptop
psi4.set_num_threads(8)      # adjust to your CPU
psi4.core.set_output_file("psi4.out", False)

# ---------------------------
# RDKit helpers
# ---------------------------
def make_3d_mol(smiles: str) -> Chem.Mol:
    """Return an RDKit Mol with Hs, 3D coords (ETKDG) and UFF minimization."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit failed to parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) < 0:
        # Fall back to older embedding heuristics
        AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    AllChem.UFFOptimizeMolecule(mol)
    return mol

def rdkit_mol_to_psi4_geometry_from_mol(mol: Chem.Mol, charge: int = 0, multiplicity: int = 1) -> str:
    """Convert an RDKit Mol (with 3D conformer) to a Psi4 geometry string."""
    conf = mol.GetConformer()
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    geom_lines = [f"{s} {c.x:.6f} {c.y:.6f} {c.z:.6f}" for s, c in zip(symbols, coords)]
    return f"{charge} {multiplicity}\n" + "\n".join(geom_lines) + "\n"

def rdkit_mol_to_psi4_geometry(smiles: str, charge: int = 0, multiplicity: int = 1) -> str:
    """Convenience shim (builds Mol then converts to Psi4 geom)."""
    return rdkit_mol_to_psi4_geometry_from_mol(make_3d_mol(smiles), charge, multiplicity)

# ---------------------------
# Psi4 energy (PCM)
# ---------------------------
def psi4_energy_in_solvent(geom: str, area: float = 0.1) -> float:
    """
    Compute single-point energy with PCM solvent (Acetonitrile).
    Returns energy in J/mol.
    """
    # Create PCM block dynamically based on cavity area
    PCM_BLOCK = PCM_TEMPLATE.format(area=area)

    psi4.set_options({
        "scf_type": "df",
        "pcm": True,
        "pcm_scf_type": "total",
        "PCM__INPUT": PCM_BLOCK,  # Pass formatted block
    })

    mol = psi4.geometry(geom)
    e_h = psi4.energy(f"{METHOD}/{BASIS}", molecule=mol)  # Hartree
    HARTREE_TO_J_MOL = 4.3597447222071e-18 * 6.02214076e23
    return e_h * HARTREE_TO_J_MOL

# ---------------------------
# Redox computation removed ", 0.05"
# ---------------------------
def compute_redox(smiles_ox: str, smiles_red: str,
                  n_electrons: int = 2, repeats: int = 3, areas=(0.05, 0.1)) -> tuple[float, float]:
    """
    Repeat energies across 'repeats' and multiple cavity 'areas' to estimate variability.
    Returns mean(E_abs), std(E_abs).
    """
    geom_ox = rdkit_mol_to_psi4_geometry(smiles_ox, 0, 1)
    geom_red = rdkit_mol_to_psi4_geometry(smiles_red, 0, 1)
    E_abs_values = []
    for _ in range(repeats):
        for area in areas:
            G_ox = psi4_energy_in_solvent(geom_ox, area)
            G_red = psi4_energy_in_solvent(geom_red, area)
            deltaG = G_red - G_ox
            E_abs = -deltaG / (n_electrons * FARADAY)
            E_abs_values.append(E_abs)
    return float(np.mean(E_abs_values)), float(np.std(E_abs_values))

# ---------------------------
# Quinone sets
# ---------------------------
# ("Anthraquinone/Anthrahydroquinone",          "O=C1c2ccccc2C(=O)c3ccccc13", "Oc1ccc2c(c1O)ccc3ccccc23"),
#    ("2-Me-Benzoquinone/2-Me-Hydroquinone",       "Cc1ccc(=O)c(=O)c1",           "Cc1ccc(O)cc1O"),
#    ("Lawsone/Reduced Lawsone",                   "O=C1C(=O)c2ccccc2C1O",        "Oc1ccc2c(c1O)ccc(=O)c2O"),
quinone_pairs = [
    ("Benzoquinone/Hydroquinone",                 "O=C1C=CC(=O)C=C1",           "Oc1ccc(O)cc1"),
    ("Naphthoquinone/Dihydroxynaphthalene",       "O=C1C=CC(=O)c2ccccc12",       "Oc1ccc2ccccc2c1O"),
   
]

# ---------------------------
# Energy table + plots removec ", 0.05"
# ---------------------------
results = []
for name, ox, red in quinone_pairs:
    mean_E, std_E = compute_redox(ox, red, n_electrons=2, repeats=3, areas=(0.05, 0.1))
    results.append({
        "System": name,
        "SMILES_ox": ox,
        "SMILES_red": red,
        "E_abs_mean": mean_E,
        "E_abs_std": std_E,
        "E_vs_SHE_mean": mean_E - 4.44,     # simple reference shift; solvent/conditions dependent
        "E_vs_SHE_std": std_E               # same std as E_abs (constant shift)
    })

df = pd.DataFrame(results)
df.to_excel("quinone_redox_with_error.xlsx", index=False)
df.to_csv("quinone_redox_with_error.csv", index=False)
print("Saved: quinone_redox_with_error.xlsx / .csv")


# ---------------------------
# Save arrays for plotting
# ---------------------------
mean_E = df["E_abs_mean"].to_numpy()
std_E = df["E_abs_std"].to_numpy()
np.savez("cavity.npz", mean_E=mean_E, std_E=std_E)
print("Saved cavity.npz with mean_E and std_E for plotting.")
