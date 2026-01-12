
# Required libraries
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

# Placeholder for results
results = []

# pH range
pH_values = [3,4,5,6,7,8,9,10]

for name, ox_smiles, red_smiles in quinone_pairs:
    ox_mol = Chem.MolFromSmiles(ox_smiles)
    red_mol = Chem.MolFromSmiles(red_smiles)

    # RDKit descriptors
    ox_logP = Crippen.MolLogP(ox_mol)
    red_logP = Crippen.MolLogP(red_mol)
    ox_MW = Descriptors.MolWt(ox_mol)
    red_MW = Descriptors.MolWt(red_mol)

    # Predict solubility using DeepChem (simplified)
    # Normally you'd featurize and predict:
    featurizer = dc.feat.ConvMolFeaturizer()
    ox_feat = featurizer.featurize([ox_mol])
    red_feat = featurizer.featurize([red_mol])
    ox_sol = model.predict_on_batch(ox_feat)
    red_sol = model.predict_on_batch(red_feat)

    # Simulate pH effect (placeholder: hydroquinone more soluble at low pH)
    sol_by_pH = []
    for pH in pH_values:
        # Simple heuristic: adjust solubility by pH (for demo)
        factor = 1 + (7 - pH) * 0.1 if pH < 7 else 1 - (pH - 7) * 0.05
        sol_by_pH.append((ox_sol[0][0]*factor, red_sol[0][0]*factor))

    results.append((name, sol_by_pH))

# Plot
plt.figure(figsize=(10,6))
for name, sol_by_pH in results:
    ox_vals = [x[0] for x in sol_by_pH]
    red_vals = [x[1] for x in sol_by_pH]
    plt.plot(pH_values, ox_vals, label=f"{name} (oxidised)")
    plt.plot(pH_values, red_vals, linestyle='--', label=f"{name} (reduced)")

plt.xlabel("pH")
plt.ylabel("Predicted Solubility (logS)")
plt.title("Solubility vs pH for Quinone Pairs")
plt.legend()
plt.grid(True)
plt.show()
