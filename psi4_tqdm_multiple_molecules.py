
from tqdm import tqdm
import psi4

# Define molecules for different redox states
molecules = [mol1, mol2, mol3]  # e.g., quinone, semiquinone, hydroquinone

print("Starting calculations...")
for mol in tqdm(molecules, desc="Computing energies", unit="state"):
    psi4.optimize('B3LYP/6-31G(d)', molecule=mol)
    energy = psi4.energy('B3LYP/6-31G(d)', molecule=mol)
    print(f"Energy for {mol.name()} = {energy:.6f} Hartree")

# Define molecules from your quinone pairs
mol1_ = "Benzoquinone/Hydroquinone"
mol1_ox = "O=C1C=CC(=O)C=C1"
mol1_red = "Oc1ccc(O)cc1"

mol2 = "Naphthoquinone/Dihydroxynaphthalene"
mol2_ox = "O=C1C=CC(=O)c2ccccc12"
mol2_red = "Oc1ccc2ccccc2c1O"

mol3 = "Anthraquinone/Anthrahydroquinone"
mol3_ox = "O=C1c2ccccc2C(=O)c3ccccc13"
mol3_red = "Oc1ccc2c(c1O)ccc3ccccc23"

mol4 = "2-Me-Benzoquinone/2-Me-Hydroquinone"
mol4_ox = "Cc1ccc(=O)c(=O)c1"
mol4_red = "Cc1ccc(O)cc1O"

mol5 = "Lawsone/red Lawsone"
mol5_ox = "O=C1C(=O)c2ccccc2C1O"
mol5_red = "Oc1ccc2c(c1O)ccc(=O)c2O"

quinone_pairs = [
    ("Benzoquinone/Hydroquinone",           "O=C1C=CC(=O)C=C1",      "Oc1ccc(O)cc1"),
    ("Naphthoquinone/Dihydroxynaphthalene", "O=C1C=CC(=O)c2ccccc12", "Oc1ccc2ccccc2c1O"),
    ("Anthraquinone/Anthrahydroquinone",          "O=C1c2ccccc2C(=O)c3ccccc13", "Oc1ccc2c(c1O)ccc3ccccc23"),
    ("2-Me-Benzoquinone/2-Me-Hydroquinone",       "Cc1ccc(=O)c(=O)c1",           "Cc1ccc(O)cc1O"),
    ("Lawsone/Reduced Lawsone",                   "O=C1C(=O)c2ccccc2C1O",        "Oc1ccc2c(c1O)ccc(=O)c2O"),
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