
# Install dependencies if needed:
# pip install deepchem mordred rdkit scikit-learn matplotlib pandas

import deepchem as dc
# from deepchem.feat import MorganFingerprintGenerator
from mordred import Calculator, descriptors
from rdkit import Chem
from mordred import Calculator, descriptors
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

# ---------------------------
# Step 1: Load ESOL (Delaney) dataset using DeepChem
# ---------------------------

tasks, datasets, transformers = dc.molnet.load_delaney(featurizer='ECFP')
train_dataset, valid_dataset, test_dataset = datasets
print(f"Delaney (ESOL) dataset loaded: {len(train_dataset)} molecules")

# Convert ESOL dataset to Pandas for Mordred descriptors
esol_smiles = [x for x in train_dataset.ids]
esol_y = train_dataset.y.flatten()

# ---------------------------
# Step 2: Compute Mordred descriptors for ESOL
# ---------------------------
calc = Calculator(descriptors, ignore_3D=True)
esol_mols = [Chem.MolFromSmiles(s) for s in esol_smiles if Chem.MolFromSmiles(s) is not None]
print("Calculating Mordred descriptors for ESOL...")
esol_desc = [calc(mol) for mol in esol_mols]
esol_df = pd.DataFrame(esol_desc)
esol_df["solubility"] = esol_y[:len(esol_mols)]  # Align lengths

# Drop columns with NaN
esol_df = esol_df.dropna(axis=1)

# ---------------------------
# Step 3: Train QSAR model
# ---------------------------
X = esol_df.drop(columns=["solubility"])
y = esol_df["solubility"]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

model = RandomForestRegressor(n_estimators=200, random_state=42)
model.fit(X_train, y_train)
print(f"Model trained. R² on test set: {r2_score(y_test, model.predict(X_test)):.3f}")

# ---------------------------
# Step 4: Quinone pairs
# ---------------------------
quinone_pairs = [
    ("Benzoquinone/Hydroquinone", "O=C1C=CC(=O)C=C1", "Oc1ccc(O)cc1"),
    ("Naphthoquinone/Dihydroxynaphthalene", "O=C1C=CC(=O)c2ccccc12", "Oc1ccc2ccccc2c1O"),
    ("Anthraquinone/Anthrahydroquinone", "O=C1c2ccccc2C(=O)c3ccccc13", "Oc1ccc2c(c1O)ccc3ccccc23"),
    ("2-Me-Benzoquinone/2-Me-Hydroquinone", "Cc1ccc(=O)c(=O)c1", "Cc1ccc(O)cc1O"),
    ("Lawsone/Reduced Lawsone", "O=C1C(=O)c2ccccc2C1O", "Oc1ccc2c(c1O)ccc(=O)c2O")
]

pH_values = [3,4,5,6,7,8,9,10]
results = []

# ---------------------------
# Step 5: Predict solubility for quinones
# ---------------------------
print("Predicting solubility for quinone pairs...")
for name, ox_smiles, red_smiles in quinone_pairs:
    ox_mol = Chem.MolFromSmiles(ox_smiles)
    red_mol = Chem.MolFromSmiles(red_smiles)

    if ox_mol is None or red_mol is None:
        print(f"Invalid SMILES for {name}. Skipping...")
        continue

    ox_desc = pd.DataFrame([calc(ox_mol)])
    red_desc = pd.DataFrame([calc(red_mol)])

    # Align columns with training set
    ox_desc = ox_desc[X.columns].fillna(0)
    red_desc = red_desc[X.columns].fillna(0)

    ox_pred = model.predict(ox_desc)[0]
    red_pred = model.predict(red_desc)[0]

    # Apply simple pH adjustment heuristic
    for pH in pH_values:
        factor = 1 + (7 - pH) * 0.1 if pH < 7 else 1 - (pH - 7) * 0.05
        ox_adj = ox_pred * factor
        red_adj = red_pred * factor
        results.append([name, pH, ox_adj, red_adj])

# ---------------------------
# Step 6: Save results
# ---------------------------
df = pd.DataFrame(results, columns=["Pair", "pH", "Oxidised_Sol", "Reduced_Sol"])
df.to_csv("quinone_sol.csv", index=False)
print("Results saved to quinone_sol.csv")

# ---------------------------
# Step 7: Plot
# ---------------------------
plt.figure(figsize=(10,6))
for pair in df["Pair"].unique():
    subset = df[df["Pair"] == pair]
    plt.plot(subset["pH"], subset["Oxidised_Sol"], label=f"{pair} (oxidised)")
    plt.plot(subset["pH"], subset["Reduced_Sol"], linestyle='--', label=f"{pair} (reduced)")

plt.xlabel("pH")
plt.ylabel("Predicted Solubility (logS)")
plt.title("Solubility vs pH for Quinone Pairs (QSAR Model)")
plt.legend()
plt.grid(True)
plt.savefig("quinone_sol_plot.png", dpi=300)
plt.show()
