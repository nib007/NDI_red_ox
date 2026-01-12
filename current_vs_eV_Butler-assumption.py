"""Computed redox potentials (E_abs_mean and E_vs_SHE_mean) for quinones.
These represent the equilibrium potential for the 2H⁺ + 2e⁻ redox couple.


✅ What You Need for Current vs Potential
To plot current (I) vs potential (E), you need a model for electron transfer kinetics. Common approaches:


Butler–Volmer Equation (most realistic for electrochemistry):
I=I0[eαnF(E−Eeq)/RT−e−(1−α)nF(E−Eeq)/RT]I = I_0 \left[ e^{\alpha n F (E - E_{eq}) / RT} - e^{-(1-\alpha) n F (E - E_{eq}) / RT} \right]I=I0​[eαnF(E−Eeq​)/RT−e−(1−α)nF(E−Eeq​)/RT]
Where:

I0I_0I0​ = exchange current density
α\alphaα = transfer coefficient (≈ 0.5)
nnn = number of electrons (2)
FFF = Faraday constant
RRR = gas constant
TTT = temperature (298 K)
EeqE_{eq}Eeq​ = equilibrium potential (your computed value)



Simplified Linear Sweep:

Assume current increases linearly with overpotential:
I∝(E−Eeq)I \propto (E - E_{eq})I∝(E−Eeq​)"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load your Excel file
df = pd.read_excel("quinone_redox_with_error_2q.xlsx")

# Constants
F = 96485.33212  # C/mol
R = 8.314        # J/mol·K
T = 298          # K
alpha = 0.5
n = 2
I0 = 1e-6        # A (example exchange current)

# quinnones
quinone_pairs = [
    ("Benzoquinone/Hydroquinone", "O=C1C=CC(=O)C=C1", "Oc1ccc(O)cc1"),
    ("Naphthoquinone/Dihydroxynaphthalene", "O=C1C=CC(=O)c2ccccc12", "Oc1ccc2ccccc2c1O"),
]

# Potential range (V)
E_range = np.linspace(-0.5, 1.0, 200)  # adjust as needed

plt.figure(figsize=(8, 6), dpi=140)

for idx, row in df.iterrows():
    E_eq = row["E_vs_SHE_mean"]  # equilibrium potential vs SHE
    I = I0 * (np.exp(alpha * n * F * (E_range - E_eq) / (R * T)) -
              np.exp(-(1 - alpha) * n * F * (E_range - E_eq) / (R * T)))
    plt.plot(E_range, I * 1e6, label=row["System"])  # µA for readability

plt.axhline(0, color="black", linewidth=0.8)
plt.xlabel("Potential vs SHE (V)")
plt.ylabel("Current (µA)")
plt.title("Simulated Butler–Volmer Curves for Quinone Redox")
plt.legend()
plt.tight_layout()
plt.savefig("quinone_current_vs_potential.png", dpi=300)
plt.show()
