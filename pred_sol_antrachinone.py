
# Install dependencies if needed:
# pip install deepchem mordred rdkit scikit-learn matplotlib pandas

import deepchem as dc
from mordred import Calculator, descriptors
from rdkit import Chem
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

# ---------------------------
# Utility: Validate SMILES
# ---------------------------
def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES: {smiles}")
    return mol

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
esol_mols = [validate_smiles(s) for s in esol_smiles if validate_smiles(s) is not None]
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
# Feature Importance
# ---------------------------
importances = pd.Series(model.feature_importances_, index=X.columns)
top20 = importances.sort_values(ascending=False).head(20)

# Print to console
print("\nTop 20 Most Important Features:")
print(top20)

# Save to CSV for easy reading
top20.to_csv("top20_features.csv", header=["Importance"])
print("Top 20 features saved to top20_features.csv")

# ---------------------------
# Plot Top 20 Features
# ---------------------------
plt.figure(figsize=(10, 6))
top20.sort_values().plot(kind='barh', color='skyblue')
plt.title("Top 20 Most Important Features")
plt.xlabel("Importance")
plt.ylabel("Feature")
plt.tight_layout()
plt.savefig("top20_features.png", dpi=300)
plt.show()
print("Bar plot saved as top20_features.png")


# ---------------------------
# Step 4: Quinone pairs
# ---------------------------
quinone_pairs = [
    ("Benzoquinone/Hydroquinone", "O=C1C=CC(=O)C=C1", "Oc1ccc(O)cc1"),

    # Anthracene-based systems
    ("Anthracene/Anthcenediol", 
     "[N+](C)(CCCN3C(=O)C1=CC(=C4C2=C(C=C(C(=C12)C3=O)Br)C(=O)N(CCC[N+](C)(C)C)C4=O)Br)(C)C", 
     "[N+](C)(C)(CCC[N+]4=C(C1=C(Br)C=C3C2=C1C(=CC(=C2C(=[N+](C3=O)CC[N+](C)(C)C)O[H])Br)C4=O)O[H])C"),

    ("Anthracediene/Anthracedinediol", 
     "[N+](C)(C=CCN3C(=O)C1=CC(=C4C2=C(C=C(C(=C12)C3=O)Br)C(=O)N(C=CC[N+](C)(C)C)C4=O)Br)(C)C", 
     "[N+](C)(C)(C=C[N+]4=C(C1=C(Br)C=C3C2=C1C(=CC(=C2C(=[N+](C3=O)C=C[N+](C)(C)C)O[H])Br)C4=O)O[H])C"),

    ("AnthraceAmidiene/AnthraceAmidinediol", 
     "[N+](C)(C=NCN3C(=O)C1=CC(=C4C2=C(C=C(C(=C12)C3=O)Br)C(=O)N(C=NC[N+](C)(C)C)C4=O)Br)(C)C", 
     "[N+](C)(C)(C=NC[N+]4=C(C1=C(Br)C=C3C2=C1C(=CC(=C2C(=[N+](C3=O)CN=C[N+](C)(C)C)O[H])Br)C4=O)O[H])C"),

    ("AnthraceneI/AnthraceneIdiol", 
     "[N+](C)(CCCN3C(=O)C1=CC(=C4C2=C(C=C(C(=C12)C3=O)I)C(=O)N(CCC[N+](C)(C)C)C4=O)I)(C)C", 
     "[N+](C)(C)(CCC[N+]4=C(C1=C(I)C=C3C2=C1C(=CC(=C2C(=[N+](C3=O)CC[N+](C)(C)C)O[H])I)C4=O)O[H])C")
]

pH_values = [3,4,5,6,7,8,9,10]
results = []

# ---------------------------
# Step 5: Predict solubility for quinones
# ---------------------------
print("Predicting solubility for quinone pairs...")
for name, ox_smiles, red_smiles in quinone_pairs:
    ox_mol = validate_smiles(ox_smiles)
    red_mol = validate_smiles(red_smiles)

    if ox_mol is None or red_mol is None:
        print(f"Skipping {name} due to invalid SMILES.")
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
df.to_csv("Anthracene_5st_sol.csv", index=False)
print("Results saved to Anthracene_5st_sol.csv")

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
