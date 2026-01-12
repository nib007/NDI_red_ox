
# -*- coding: utf-8 -*-
"""
Quinone redox: energies (G in water) + pH transform + CV simulation
- Computes E° vs SHE at pH=5 for BQ/HQ and other pairs (B3LYP-D3BJ/def2-SVP + PCM(Water))
- Produces Excel/CSV results and PNG CV figures at selected scan rates.

Notes:
- For higher accuracy, consider a larger basis (def2-TZVP) and/or COSMO/SMD; calibrate vs a reference couple.
- The CV simulator uses semi-infinite diffusion with Butler–Volmer kinetics (quasi-reversible).
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
FARADAY = 96485.33212        # C/mol
R = 8.314462618              # J/mol/K
T = 298.15                   # K

METHOD = "B3LYP-D3BJ"
BASIS = "def2-SVP"           # consider def2-TZVP for production runs

# Absolute SHE potential (IUPAC recommendation ~4.44 V at 25 °C)
EABS_SHE = 4.44

# pH settings for PCET (default BQ/HQ: m=2 H+, n=2 e-)
PH_TARGET = 5.0
M_PROTONS = 2
N_ELECTRONS = 2

# PCM template with escaped braces for str.format(), water solvent
PCM_TEMPLATE = """
Units = Angstrom
Medium {{
    SolverType = IEFPCM
    Solvent = Water
}}
Cavity {{
    RadiiSet = UFF
    Type = GePol
    Area = {area}
    Mode = Implicit
}}
"""

# Psi4 runtime settings
psi4.set_memory('12 GB')
psi4.set_num_threads(8)
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
# Psi4: Gibbs free energy (PCM, water)
# ---------------------------
def psi4_free_energy_in_water(geom: str, area: float = 0.1,
                              method: str = METHOD, basis: str = BASIS) -> float:
    """
    Optimize + frequency in water PCM to get Gibbs free energy (J/mol).
    Falls back to electronic energy if G is not exposed (older Psi4 versions).
    """
    PCM_BLOCK = PCM_TEMPLATE.format(area=area)
    psi4.set_options({
        "scf_type": "df",
        "pcm": True,
        "pcm_scf_type": "total",
        "PCM__INPUT": PCM_BLOCK,
    })
    mol = psi4.geometry(geom)

    # Optimize geometry in solvent
    print("Starting geometry optimisation...")
    e_opt_h = psi4.optimize(f"{method}/{basis}", molecule=mol)
    print("Geometry optimisation complete.")

    # Frequency for thermal corrections (RRHO)
    try:
        psi4.frequency(f"{method}/{basis}", molecule=mol)
        # Psi4 thermochemistry stores variable; try to fetch Gibbs free energy
        G_h = psi4.variable("GIBBS FREE ENERGY")  # Hartree, if available
        HARTREE_TO_J_MOL = 4.3597447222071e-18 * 6.02214076e23
        return float(G_h) * HARTREE_TO_J_MOL
    except Exception:
        # Fallback: use optimized electronic energy (less accurate)
        HARTREE_TO_J_MOL = 4.3597447222071e-18 * 6.02214076e23
        return float(e_opt_h) * HARTREE_TO_J_MOL

# ---------------------------
# Thermodynamics: pH transform & referencing
# ---------------------------
def pH_transform_and_reference(E_abs: float,
                               m_protons: int = M_PROTONS,
                               n_electrons: int = N_ELECTRONS,
                               pH: float = PH_TARGET,
                               Eabs_SHE: float = EABS_SHE) -> float:
    """
    Convert absolute potential (vacuum) to vs SHE at target pH:
    E_vs_SHE(pH) = (E_abs - Eabs_SHE) - (m/n)*0.05916*pH   [V at 25 °C]
    """
    E_vs_SHE_pH0 = E_abs - Eabs_SHE
    E_vs_SHE_pH = E_vs_SHE_pH0 - (m_protons / n_electrons) * 0.05916 * pH
    return E_vs_SHE_pH

# ---------------------------
# Redox computation with error bars over cavity areas
# ---------------------------
def compute_redox_pH(smiles_ox: str, smiles_red: str,
                     n_electrons: int = N_ELECTRONS,
                     m_protons: int = M_PROTONS,
                     pH: float = PH_TARGET,
                     repeats: int = 2,
                     areas=(0.05, 0.10)) -> dict:
    """
    Compute mean/std for E_abs and E_vs_SHE at given pH by varying PCM cavity area.
    """
    geom_ox = rdkit_mol_to_psi4_geometry(smiles_ox, 0, 1)
    geom_red = rdkit_mol_to_psi4_geometry(smiles_red, 0, 1)
    E_abs_values = []
    E_vs_SHE_values = []
    for _ in range(repeats):
        for area in areas:
            G_ox = psi4_free_energy_in_water(geom_ox, area)
            G_red = psi4_free_energy_in_water(geom_red, area)
            deltaG = G_red - G_ox
            E_abs = -deltaG / (n_electrons * FARADAY)  # V (absolute scale)
            E_vs_SHE = pH_transform_and_reference(E_abs, m_protons, n_electrons, pH, EABS_SHE)
            E_abs_values.append(E_abs)
            E_vs_SHE_values.append(E_vs_SHE)
    return {
        "E_abs_mean": float(np.mean(E_abs_values)),
        "E_abs_std": float(np.std(E_abs_values)),
        "E_vs_SHE_mean": float(np.mean(E_vs_SHE_values)),
        "E_vs_SHE_std": float(np.std(E_vs_SHE_values)),
    }

# ---------------------------
# Quinone sets
# ---------------------------
quinone_pairs = [
    ("Benzoquinone/Hydroquinone",           "O=C1C=CC(=O)C=C1",      "Oc1ccc(O)cc1"),
    ("Naphthoquinone/Dihydroxynaphthalene", "O=C1C=CC(=O)c2ccccc12", "Oc1ccc2ccccc2c1O"),
]

# ---------------------------
# Build energy table
# ---------------------------
results = []
for name, ox, red in quinone_pairs:
    stats = compute_redox_pH(ox, red, n_electrons=N_ELECTRONS, m_protons=M_PROTONS,
                             pH=PH_TARGET, repeats=2, areas=(0.05, 0.10))
    results.append({
        "System": name,
        "SMILES_ox": ox,
        "SMILES_red": red,
        **stats
    })

df = pd.DataFrame(results)
df.to_excel("quinone_redox_pH5_water.xlsx", index=False)
df.to_csv("quinone_redox_pH5_water.csv", index=False)
print("Saved: quinone_redox_pH5_water.xlsx / .csv")

# ---------------------------
# CV simulation (semi-infinite diffusion + Butler–Volmer)
# ---------------------------
F = FARADAY

def triangular_wave(E_start, E_vertex, scan_rate, dt):
    # Forward
    t_f = abs((E_vertex - E_start) / scan_rate)
    n_f = int(np.ceil(t_f / dt))
    E_f = np.linspace(E_start, E_vertex, n_f, endpoint=False)
    # Reverse
    E_r = np.linspace(E_vertex, E_start, n_f, endpoint=False)
    E = np.concatenate([E_f, E_r])
    t = np.arange(len(E)) * dt
    return E, t

def simulate_cv(E0, n=2, alpha=0.5, k0=0.02,     # kinetics (cm/s)
                D_O=7e-6, D_R=7e-6,               # diffusion (cm^2/s)
                C_bulk=0.0007,                    # mol/cm^3 (0.7 M)
                A=0.07,                           # cm^2 (3 mm GC disk)
                scan_rate=0.1,                    # V/s
                E_start=-0.4, E_vertex=0.6,
                L=0.03, Nx=400, dt=1e-3):
    """
    Quasi-reversible CV for O + n e- <=> R at pH5, vs SHE.
    """
    x = np.linspace(0.0, L, Nx)
    dx = x[1] - x[0]

    CO = np.full(Nx, C_bulk)  # start all oxidized
    CR = np.zeros(Nx)

    lam_O = D_O * dt / (dx*dx)
    lam_R = D_R * dt / (dx*dx)

    def crank_nicolson_mats(lam, N):
        A = np.diag((1+2*lam)*np.ones(N)) + np.diag(-lam*np.ones(N-1),1) + np.diag(-lam*np.ones(N-1),-1)
        B = np.diag((1-2*lam)*np.ones(N)) + np.diag(lam*np.ones(N-1),1) + np.diag(lam*np.ones(N-1),-1)
        return A, B

    A_O, B_O = crank_nicolson_mats(lam_O, Nx)
    A_R, B_R = crank_nicolson_mats(lam_R, Nx)

    E, t = triangular_wave(E_start, E_vertex, scan_rate, dt)
    currents = []

    from numpy.linalg import solve

    for k in range(len(t)):
        Ek = E[k]
        eta = Ek - E0  # overpotential

        # Butler–Volmer current density (A/cm^2) at surface
        j = n*F*k0*(CO[0]*np.exp(-alpha*n*F*eta/(R*T)) - CR[0]*np.exp((1-alpha)*n*F*eta/(R*T)))

        # RHS with Neumann flux at x=0
        rhs_O = B_O.dot(CO)
        rhs_R = B_R.dot(CR)

        # Flux BC: -(D dC/dx)|0 = j/(nF)
        rhs_O[0] += 2*lam_O * dx * ( - j/(n*F*D_O) )
        rhs_R[0] += 2*lam_R * dx * ( + j/(n*F*D_R) )

        # Dirichlet at x=L
        rhs_O[-1] = C_bulk
        rhs_R[-1] = 0.0
        A_O[-1,:] = 0.0; A_O[-1,-1] = 1.0
        A_R[-1,:] = 0.0; A_R[-1,-1] = 1.0

        CO = solve(A_O, rhs_O)
        CR = solve(A_R, rhs_R)
        CO = np.clip(CO, 0.0, None)
        CR = np.clip(CR, 0.0, None)

        currents.append(j * A)

    return E, np.array(currents), t

# Choose BQ/HQ E° vs SHE at pH5 from our table
E0_bq_pH5 = df.loc[df["System"] == "Benzoquinone/Hydroquinone", "E_vs_SHE_mean"].values[0]

# Scan rates to plot (V/s): 10–500 mV/s
scan_rates = [0.01, 0.05, 0.10, 0.20, 0.50]
plt.figure(figsize=(7,5))
for v in scan_rates:
    E, i, t = simulate_cv(E0=E0_bq_pH5, n=2, alpha=0.5, k0=0.02,
                          D_O=7e-6, D_R=7e-6, C_bulk=0.0007,
                          A=0.07, scan_rate=v,
                          E_start=E0_bq_pH5-0.3, E_vertex=E0_bq_pH5+0.3,
                          L=0.03, Nx=400, dt=1e-3)
    plt.plot(E, i*1e3, lw=2, label=f"{int(v*1000)} mV/s")

plt.xlabel("Potential (V vs SHE)")
plt.ylabel("Current (mA)")
plt.title("Simulated CV of BQ/HQ at 5 in water (0.7 M)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("cv_bq_pH5_water.png", dpi=200)
print("Saved: cv_bq_pH5_water.png")
