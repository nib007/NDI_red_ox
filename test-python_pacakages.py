
#python - <<'PY'
import psi4
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
print("OK:",
      "psi4", psi4.__version__,
      "| RDKit", Chem.rdBase.rdkitVersion,
      "| pandas", pd.__version__,
      "| matplotlib", plt.matplotlib.__version__,
      "| numpy", np.__version__)
# PY
