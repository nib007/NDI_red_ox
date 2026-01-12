Solubility Analysis of Anthracene and Naphthalenediimide Derivatives
Comprehensive Evaluation of Structural Modifications, Oxidation States, and Electrochemical Implications
Summary of Anthracene and Naphthalenediimide Solubility Studies
 Objective
The primary aim of this study is to design and evaluate redox-active aromatic molecules, specifically anthracene derivatives and naphthalenediimide (NDI) derivatives. The goal is to understand how structural modifications and variations in oxidation state influence their aqueous solubility and their potential interactions with electrode surfaces. 
Systems Studied
The investigation focused on two main classes of compounds:
1.	Anthracene (sometimes referred to as “Antracene”):
2.	Both oxidized and reduced forms were examined.
3.	Variants included those functionalized with amidine, SO₂N groups, and sulfonate (SO₃) adducts.
4.	Naphthalenediimide (NDI):
5.	Both oxidized and reduced forms were considered.
6.	Structural modifications included the addition of halogens (Br, I), diamines, quaternary ammonium (N⁺(CH₃)), and sulfonate groups.
Methodology
•	Computational modeling was employed using fastsolv python package (ver 250528) to predict solubility, expressed as logS, across a temperature range of 250–350 °C (Figure 2).
•	The primary solvent used in these simulations was water (H₂O), with some alternative solvent representations, such as CC#N, included for comparison.
•	The results incorporated error estimates, specifically standard deviation, to provide a measure of confidence in the predictions.
Key Findings
Anthracene Series
•	Oxidized vs Reduced: Both forms exhibited only modest improvements in solubility with increasing temperature, and solubility remained within the low millimolar range (Figure 3).
•	Functionalization: The introduction of amidine or SO₃ groups resulted in only slight predicted improvements in logS, which were insufficient to achieve high aqueous concentrations.
•	Trend: For all anthracene variants, solubility increased with temperature; however, the rate of increase was similar across the different modifications.
Naphthalenediimide Series
•	Baseline NDI: The maximum predicted solubility in water for unmodified NDI was approximately 10 mM.
•	Variants (Br, I, diamine, N⁺(CH₃), SO₃): These structural modifications did not significantly enhance solubility relative to the baseline NDI compound.
•	Oxidized vs Reduced: Reduced forms were observed to be slightly more soluble than their oxidized counterparts, though the differences were relatively modest.
Interpretation
•	Solubility Limitation: Despite the introduction of polar substituents, the large aromatic core of these molecules continues to constrain their aqueous solubility.
•	Electrode Interaction Potential: While increased polarity through groups such as sulfonate and quaternary ammonium does not substantially improve bulk solubility, it may enhance surface affinity and electrochemical interface behavior.
•	Design Implication: Future strategies for improving these systems may include:
•	Polymeric or dendritic architectures to better disperse the molecules in solution.
•	Mixed solvent systems or the use of ionic liquids to achieve higher solubility.
•	Surface functionalization approaches that prioritize electrode compatibility over bulk solubility.
Next Steps
•	Experimental validation of the computational predictions, including solubility measurements at both room temperature and elevated temperatures.
•	Further investigation into electrode adsorption behaviors and electron transfer kinetics for the studied derivatives.
•	Exploration of co-solvent systems or formulation strategies aimed at overcoming current solubility limitations. 
 
Figure 2: Solubility vs Temp 250-370Kelvin of variuos forms of NDI (tentatively referred to as antracen_Br, _I. SO2, diene).
Figure 3: SMILES Codes and structures:
CN(C)CCN1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCN(C)C)C3=O Naphthalenediimide_ox
CN(C)CC[n+]1c(O)c2ccc3c(O)[n+](CCN(C)C)c(O)c4ccc(c1O)c2c34 Naphthalenediimide_ox
Antracene :
[N+](C)(CCCN3C(=O)C1=CC(=C4C2=C(C=C(C(=C12)C3=O)Br)C(=O)N(CCC[N+](C)(C)C)C4=O)Br)(C)C Antracene_ox
[N+](C)(C)(CCC[N+]4=C(C1=C(Br)C=C3C2=C1C(=CC(=C2C(=[N+](C3=O)CC[N+](C)(C)C)O[H])Br)C4=O)O[H])C Antracene_red
 
Antracediene:
[N+](C)(C=CCN3C(=O)C1=CC(=C4C2=C(C=C(C(=C12)C3=O)Br)C(=O)N(C=CC[N+](C)(C)C)C4=O)Br)(C)C Antracediene_ox
[N+](C)(C)(C=C[N+]4=C(C1=C(Br)C=C3C2=C1C(=CC(=C2C(=[N+](C3=O)C=C[N+](C)(C)C)O[H])Br)C4=O)O[H])C Antracediene_red
 
AntraceAmidiene:
[N+](C)(C=NCN3C(=O)C1=CC(=C4C2=C(C=C(C(=C12)C3=O)Br)C(=O)N(C=NC[N+](C)(C)C)C4=O)Br)(C)C AntraceAmidiene_ox
[N+](C)(C)(C=NC[N+]4=C(C1=C(Br)C=C3C2=C1C(=CC(=C2C(=[N+](C3=O)CN=C[N+](C)(C)C)O[H])Br)C4=O)O[H])C AntraceAmidiene_red
 
AntraceneI:
[N+](C)(CCCN3C(=O)C1=CC(=C4C2=C(C=C(C(=C12)C3=O)I)C(=O)N(CCC[N+](C)(C)C)C4=O)I)(C)C AntraceneI_ox
[N+](C)(C)(CCC[N+]4=C(C1=C(I)C=C3C2=C1C(=CC(=C2C(=[N+](C3=O)CC[N+](C)(C)C)O[H])I)C4=O)O[H])C AntraceneI_red
 
New sulphate version:
[N+](C)(C=NCN3C(=O)C1=CC(=C4C2=C(C=C(C(=C12)C3=O)(S(=O)(=O)[O-]))C(=O)N(C=NC[N+](C)(C)C)C4=O)(S(=O)(=O)[O-]))(C)C AnthracenAmiddieneSO3_ox
[N+](C)(C)(C=NC[N+]4=C(C1=C(S(=O)(=O)[O-])C=C3C2=C1C(=CC(=C2C(=[N+](C3=O)CN=C[N+](C)(C)C)O[H])(S(=O)(=O)[O-]))C4=O)O[H])C AnthracenAmiddieneSO3_red
S(=O)(=O)N(CC)C
[N+](C)(C=NCN3C(=O)C1=CC(=C4C2=C(C=C(C(=C12)C3=O)( S(=O)(=O)N(CC)C))C(=O)N(C=NC[N+](C)(C)C)C4=O)( S(=O)(=O)N(CC)C))(C)C AntraceAmidieneSO2N_ox
[N+](C)(C)(C=NC[N+]4=C(C1=C(S(=O)(=O)N(CC)C)C=C3C2=C1C(=CC(=C2C(=[N+](C3=O)CN=C[N+](C)(C)C)O[H])( S(=O)(=O)N(CC)C))C4=O)O[H])C AntraceAmidieneSO2N_red

Daylight smiles:
O=S1(OC(C2C3C1=CC1=C4C=3C(=CC3S(=O)(=O)OC(=O)C(=[N+](C/N=C/[N+](C)(C)C)C1=O)C=34)C(=O)[N+]=2/C=N/C[N+](C)(C)C)=O)=O AntraceAmidieneSO2N_ox
O1S(=O)(=O)C2C=C3C(=[N+](/C=C/C[N+](C)(C)C)C4C(OS(=O)(=O)C5=C(C6C(=C3C=45)C=2C(=[N+](C/N=C/[N+](C)(C)C)C=6O[H])C1=O)O[H])=O)O[H]
AntraceAmidieneSO2N_red
AntraceAmidieneSO2N_ox   
 
AntraceAmidieneSO2N_red:
 
 
Anthracene_sulphate_ox:
O=c1c2c(O)c(Br)c4c(=O)n(CCS(=O)(=O)O)c(=O)c3c(O)c(Br)c(c(=O)n1CCCS(=O)(=O)O)c2c34
 

