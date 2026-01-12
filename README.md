
# Solubility Analysis of Anthracene and Naphthalenediimide Derivatives  
[![License: MIT](https://img.shields
_Comprehensive Evaluation of Structural Modifications, Oxidation States, and Electrochemical Implications_

---

## 📌 Overview
This repository contains computational and analytical work on **redox-active aromatic molecules**, focusing on **anthracene derivatives** and **naphthalenediimide (NDI) derivatives**.  
The goal: **Understand how structural modifications and oxidation states influence aqueous solubility and electrode interactions.**

---

## 🔬 Systems Studied
- **Anthracene (Antracene)**  
  - Oxidized and reduced forms  
  - Variants: amidine, SO₂N groups, sulfonate (SO₃) adducts  

- **Naphthalenediimide (NDI)**  
  - Oxidized and reduced forms  
  - Variants: halogens (Br, I), diamines, quaternary ammonium (N⁺(CH₃)), sulfonate groups  

---

## ⚙️ Methodology
- **Computational modeling** using `fastsolv` (ver 250528)
- Predicted **logS vs temperature** (250–350 °C)
- Primary solvent: **water (H₂O)**; alternative solvent: **CC#N**
- Error estimates included (standard deviation)

---

## 📈 Key Findings
| Series       | Observation |
|--------------|-------------|
| **Anthracene** | Solubility remains low (mM range); slight improvement with temperature |
| **NDI**       | Baseline solubility ≈ 10 mM; modifications (Br, I, SO₃, N⁺) do not significantly increase solubility |
| **Ox vs Red** | Reduced forms slightly more soluble than oxidized |

---

## 💡 Interpretation
- Large aromatic cores limit aqueous solubility despite polar substituents.
- Increased polarity may improve **electrode surface affinity**, not bulk solubility.
- Future strategies:
  - Polymeric/dendritic architectures
  - Mixed solvents or ionic liquids
  - Surface functionalization for electrode compatibility

---

## ✅ Next Steps
- Experimental validation of predictions
- Study electrode adsorption and electron transfer kinetics
- Explore co-solvent systems or formulation strategies

---

## 📂 Repository Structure
