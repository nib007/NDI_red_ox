
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
import multiprocessing

# ---------------------------
# Utility: Validate SMILES
# ---------------------------
def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol

if __name__ == "__main__":
    # ---------------------------
    # Step 1: Load ESOL (Delaney) dataset using DeepChem
    # ---------------------------
    tasks, datasets, transformers = dc.molnet.load_delaney(featurizer='ECFP')
    train_dataset, valid_dataset, test_dataset = datasets  # ✅ Unpack datasets
    print(f"Delaney (ESOL) dataset loaded: {len(train_dataset)} molecules")

    # Convert ESOL dataset to Pandas for Mordred descriptors
    esol_smiles = [x for x in train_dataset.ids]
    esol_y = train_dataset.y.flatten()

    # ---------------------------
    # Step 2: Compute Mordred descriptors (parallel for speed)
    # ---------------------------
    calc = Calculator(descriptors, ignore_3D=True)
    nproc = multiprocessing.cpu_count()  # Use all available cores
    print(f"Using {nproc} CPU cores for Mordred calculation...")

    # Validate SMILES and create RDKit molecules
    esol_mols = [validate_smiles(s) for s in esol_smiles if validate_smiles(s) is not None]
    print(f"Valid molecules: {len(esol_mols)}")

    # Compute descriptors in parallel
    print("Calculating Mordred descriptors for ESOL...")
    esol_df = calc.pandas(esol_mols, nproc=nproc)  # ✅ Parallel computation
    esol_df["solubility"] = esol_y[:len(esol_mols)]  # Align lengths

    # Drop columns with NaN
    esol_df = esol_df.dropna(axis=1)

    # ---------------------------
    # Step 3: Train QSAR model
    # ---------------------------
    X = esol_df.drop(columns=["solubility"])
    y = esol_df["solubility"]

    # Save descriptor names for later mapping
    pd.DataFrame(X.columns, columns=["Descriptor"]).to_csv("descriptor_names.csv", index=False)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    model = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
    model.fit(X_train, y_train)
    print(f"Model trained. R² on test set: {r2_score(y_test, model.predict(X_test)):.3f}")

    # ---------------------------
    # Step 4: Feature Importance
    # ---------------------------
    importances = pd.Series(model.feature_importances_, index=X.columns)
    top20 = importances.sort_values(ascending=False).head(20)

    # Print to console
    print("\nTop 20 Most Important Features wtih name:")
    print(top20)

    # Save raw top 20 to CSV
    top20.to_csv("top20_features_named.csv", header=["Importance"])
    print("Top 20 features saved to top20_features_named.csv")

    # Create DataFrame with index, name, importance
    top20_df = pd.DataFrame({
        "Index": [X.columns.get_loc(col) for col in top20.index],
        "Descriptor": top20.index,
        "Importance": top20.values
    })

    # ---------------------------
    # Step 5: Add Descriptor Category & Meaning
    # ---------------------------
    descriptor_info = {}
    descriptor_category = {}

    for d in descriptors.__all__:
        desc_class = getattr(descriptors, d)
        if hasattr(desc_class, "__doc__"):
            descriptor_info[d] = desc_class.__doc__
        descriptor_category[d] = d.split("_")[0]

    top20_df["Category"] = top20_df["Descriptor"].apply(lambda name: descriptor_category.get(name.split("_")[0], "Unknown"))
    top20_df["Meaning"] = top20_df["Descriptor"].apply(lambda name: descriptor_info.get(name.split("_")[0], "Description not found"))

    # Save to CSV
    top20_df.to_csv("top20_features_named.csv", index=False)
    print("✅ Top 20 features saved to top20_features_named.csv")

    # ---------------------------
    # Step 6: Plot Top 20 Features
    # ---------------------------
    plt.figure(figsize=(10, 6))
    top20.sort_values().plot(kind='barh', color='skyblue')
    plt.title("Top 20 Most Important Features")
    plt.xlabel("Importance")
    plt.ylabel("Feature")
    plt.tight_layout()
    plt.savefig("top20_features_named.png", dpi=300)
    plt.show()
    print("✅ Bar plot saved as top20_features_named.png")
