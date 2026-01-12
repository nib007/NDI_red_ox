
"""
Quinone redox: read results and plot
- Reads cavity_2q.npz and quinone_redox_with_error_2q.xlsx
- Produces bar chart and 3D overlay plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from matplotlib.lines import Line2D

# ---------------------------
# Load data
# ---------------------------
excel_file = "quinone_redox_with_error_2q.xlsx"
npz_file = "cavity_2q.npz"

df = pd.read_excel(excel_file)
data = np.load(npz_file)
mean_E = data["mean_E"]
std_E = data["std_E"]

print("Loaded data:")
print(df)

# ---------------------------
# Bar chart with error bars
# ---------------------------
plt.figure(figsize=(8, 5), dpi=140)
x = np.arange(len(df))
w = 0.35
plt.bar(x - w/2, df["E_abs_mean"], yerr=df["E_abs_std"], capsize=5,
        color="#1f77b4", label="E_abs (V)")
plt.bar(x + w/2, df["E_vs_SHE_mean"], yerr=df["E_vs_SHE_std"], capsize=5,
        color="#ff7f0e", alpha=0.8, label="E_vs_SHE (V)")
plt.xticks(x, df["System"], rotation=20, ha="right")
plt.ylabel("Potential (V)")
plt.title("Quinone Redox Potentials (PCM, MeCN)")
plt.legend()
plt.tight_layout()
plt.savefig("quinone_redox_errorbars.png", dpi=300)
plt.show()
print("Saved: quinone_redox_errorbars.png")

# ---------------------------
# RDKit helpers for 3D overlay
# ---------------------------
def make_3d_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    return mol

def plot_mol3d(ax, mol, color="#1f77b4", alpha=0.9):
    conf = mol.GetConformer()
    pts = np.array([[conf.GetAtomPosition(i).x,
                     conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=30, color=color, alpha=alpha)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        xs, ys, zs = [pts[i, 0], pts[j, 0]], [pts[i, 1], pts[j, 1]], [pts[i, 2], pts[j, 2]]
        ax.plot(xs, ys, zs, color=color, linewidth=1.5, alpha=alpha)

def align_mol(mobile, ref):
    try:
        rdMolAlign.AlignMol(mobile, ref)
    except:
        pass

def mol_points(mol):
    c = mol.GetConformer()
    return np.array([[c.GetAtomPosition(i).x, c.GetAtomPosition(i).y, c.GetAtomPosition(i).z]
                     for i in range(mol.GetNumAtoms())])

def set_equal_aspect_3d(ax, pts):
    mins, maxs = pts.min(axis=0), pts.max(axis=0)
    centers = (maxs + mins) / 2
    ranges = (maxs - mins)
    max_range = ranges.max() * 0.6
    ax.set_xlim(centers[0] - max_range, centers[0] + max_range)
    ax.set_ylim(centers[1] - max_range, centers[1] + max_range)
    ax.set_zlim(centers[2] - max_range, centers[2] + max_range)

# ---------------------------
# Overlay plots for two quinones
# ---------------------------
quinone_pairs = [
    ("Benzoquinone/Hydroquinone", "O=C1C=CC(=O)C=C1", "Oc1ccc(O)cc1"),
    ("Naphthoquinone/Dihydroxynaphthalene", "O=C1C=CC(=O)c2ccccc12", "Oc1ccc2ccccc2c1O"),
]

for name, ox, red in quinone_pairs:
    mol_ox = make_3d_mol(ox)
    mol_red = make_3d_mol(red)
    align_mol(mol_red, mol_ox)

    fig = plt.figure(figsize=(7, 6), dpi=140)
    ax = fig.add_subplot(111, projection="3d")
    plot_mol3d(ax, mol_ox, color="#1f77b4", alpha=0.95)
    plot_mol3d(ax, mol_red, color="#ff7f0e", alpha=0.8)

    proxy = [Line2D([0],[0], color="#1f77b4", lw=3),
             Line2D([0],[0], color="#ff7f0e", lw=3)]
    ax.legend(proxy, ["Oxidized", "Reduced"], loc="upper left")

    both = np.vstack([mol_points(mol_ox), mol_points(mol_red)])
    set_equal_aspect_3d(ax, both)
    ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)"); ax.set_zlabel("Z (Å)")
    ax.set_title(f"Overlay: {name}")
    plt.tight_layout()
    outpng = f"overlay_{name.replace('/', '_').replace(' ', '_')}.png"
    plt.savefig(outpng, dpi=300)
    plt.show()
    print(f"Saved overlay: {outpng}")

print("Done.")
