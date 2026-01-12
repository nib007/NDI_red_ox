
# Install dependencies if needed:
# pip install rdkit deepchem mordred matplotlib pandas

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
import deepchem as dc
import matplotlib.pyplot as plt
import pandas as pd

# Quinone pairs
quinone_pairs = [
    ("Benzoquinone/Hydroquinone", "O=C1C=CC(=O)C=C1", "Oc1ccc(O)cc1"),
    ("Naphthoquinone/Dihydroxynaphthalene", "O=C1C=CC(=O)c2ccccc12", "Oc1ccc2ccccc2c1O"),
    ("Anthraquinone/Anthrahydroquinone", "O=C1c2ccccc2C(=O)c3ccccc13", "Oc1ccc2c(c1O)ccc3ccccc23"),
    ("2-Me-Benzoquinone/2-Me-Hydroquinone", "Cc1ccc(=O)c(=O)c1", "Cc1ccc(O)cc1O"),
    ("Lawsone/Reduced Lawsone", "O=C1C(=O)c2ccccc2C1O", "Oc1ccc2c(c1O)ccc(=O)c2O")
]

# DeepChem model (GraphConv for solubility)
model = dc.models.GraphConvModel(n_tasks=1, mode='regression')

# Featurizer
featurizer = dc.feat.ConvMolFeaturizer()

# pH range
pH_values = [3,4,5,6,7,8,9,10]

# Collect data
data_rows = []

for name, ox_smiles, red_smiles in quinone_pairs:
    ox_mol = Chem.MolFromSmiles(ox_smiles)
    red_mol = Chem.MolFromSmiles(red_smiles)

    # RDKit descriptors
    ox_logP = Crippen.MolLogP(ox_mol)
    red_logP = Crippen.MolLogP(red_mol)
    ox_MW = Descriptors.MolWt(ox_mol)
    red_MW = Descriptors.MolWt(red_mol)

    # DeepChem prediction
    ox_feat = featurizer.featurize([ox_mol])
    red_feat = featurizer.featurize([red_mol])
    ox_sol = model.predict_on_batch(ox_feat)[0][0]
    red_sol = model.predict_on_batch(red_feat)[0][0]

    # Simulate pH effect (placeholder heuristic)
    for pH in pH_values:
        factor = 1 + (7 - pH) * 0.1 if pH < 7 else 1 - (pH - 7) * 0.05
        ox_adj = ox_sol * factor
        red_adj = red_sol * factor
        data_rows.append([name, pH, ox_adj, red_adj, ox_logP, red_logP, ox_MW, red_MW])

# Save to CSV
df = pd.DataFrame(data_rows, columns=["Pair", "pH", "Oxidised_Sol", "Reduced_Sol", "Ox_LogP", "Red_LogP", "Ox_MW", "Red_MW"])
df.to_csv("quinone_sol.csv", index=False)
print("Data saved to quinone_sol.csv")

# Plot from saved file
df_plot = pd.read_csv("quinone_sol.csv")
plt.figure(figsize=(10,6))
for pair in df_plot["Pair"].unique():
    subset = df_plot[df_plot["Pair"] == pair]
    plt.plot(subset["pH"], subset["Oxidised_Sol"], label=f"{pair} (oxidised)")
    plt.plot(subset["pH"], subset["Reduced_Sol"], linestyle='--', label=f"{pair} (reduced)")

plt.xlabel("pH")
plt.ylabel("Predicted Solubility (logS)")
plt.title("Solubility vs pH for Quinone Pairs")
plt.legend()
plt.grid(True)
plt.savefig("quinone_sol_plot.png", dpi=300)
plt.show()
