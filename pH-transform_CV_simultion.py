
# -*- coding: utf-8 -*-
"""
Cyclic voltammetry simulator for O + n e- <=> R (planar electrode)
- 1D diffusion (Fick's 2nd law, semi-infinite)
- Butler–Volmer boundary at x=0
- Triangular potential sweep (CV)
References:
  - Randles–Ševčík relation for reversible ip vs v (sanity check).  (IUPAC Gold Book; Pine Research) 
  - Nicholson–Shain/Bard–Faulkner theory for quasi-reversible systems. (Didactic review) 
"""

import numpy as np
import matplotlib.pyplot as plt

F = 96485.33212  # C/mol
R = 8.314462618  # J/mol/K
T = 298.15       # K

def pourbaix_shift(E0_pH0, m=2, n=2, pH=14.0):
    """
    Nernst/Legendre pH transform for m H+ and n e- at 25 °C.
    E(pH) = E(pH=0) - (m/n) * 0.0591 * pH   [V]
    """
    return E0_pH0 - (m/n) * 0.0591 * pH

def triangular_wave(E_start, E_vertex, scan_rate, dt):
    """
    Generate triangular potential sweep in time: start -> vertex -> start (one cycle).
    Returns E(t) and t.
    """
    # Forward branch
    t_forward = abs((E_vertex - E_start) / scan_rate)
    n_forward = int(np.ceil(t_forward / dt))
    E_f = np.linspace(E_start, E_vertex, n_forward, endpoint=False)
    # Reverse branch
    t_reverse = t_forward
    n_reverse = n_forward
    E_r = np.linspace(E_vertex, E_start, n_reverse, endpoint=False)
    E = np.concatenate([E_f, E_r])
    t = np.arange(len(E)) * dt
    return E, t

def simulate_cv(E0, n=1, alpha=0.5, k0=0.01,  # kinetics
                D_O=1e-5, D_R=1e-5,           # cm^2/s (typical small organics ~1e-5)
                C_bulk=1e-3,                  # mol/cm^3 (1 mM)
                A=0.07,                       # cm^2 (e.g., 3 mm disk)
                scan_rate=0.1,                # V/s
                E_start=-0.2, E_vertex=0.6,  # V vs SHE
                L=0.03,                       # cm (domain; ~300 um is fine for CV timescales)
                Nx=400,                       # grid points
                dt=1e-3):                     # s time step
    """
    Quasi-reversible CV with Butler–Volmer boundary and 1D diffusion.
    Returns potential E, current i (A), and time t.
    """
    x = np.linspace(0.0, L, Nx)
    dx = x[1]-x[0]

    # Initial conditions: all oxidized in bulk (O = C_bulk, R = 0)
    CO = np.full(Nx, C_bulk)
    CR = np.zeros(Nx)

    # Precompute diffusion matrices (Crank–Nicolson)
    lam_O = D_O * dt / (dx*dx)
    lam_R = D_R * dt / (dx*dx)

    def crank_nicolson_step(C, lam):
        # Tridiagonal solver A*C^{k+1} = B*C^{k} + boundary terms later
        # Here we assemble matrices once per species using lam
        A = np.diag((1+2*lam)*np.ones(Nx)) + np.diag(-lam*np.ones(Nx-1),1) + np.diag(-lam*np.ones(Nx-1),-1)
        B = np.diag((1-2*lam)*np.ones(Nx)) + np.diag(lam*np.ones(Nx-1),1) + np.diag(lam*np.ones(Nx-1),-1)
        return A, B

    A_O, B_O = crank_nicolson_step(CO, lam_O)
    A_R, B_R = crank_nicolson_step(CR, lam_R)

    # Potential waveform
    E, t = triangular_wave(E_start, E_vertex, scan_rate, dt)

    currents = []

    # Time marching
    from numpy.linalg import solve

    for k in range(len(t)):
        Ek = E[k]
        eta = Ek - E0  # overpotential vs formal potential

        # Butler–Volmer current density using SURFACE concentrations (x=0)
        j = n*F*k0*(CO[0]*np.exp(-alpha*n*F*eta/(R*T)) - CR[0]*np.exp((1-alpha)*n*F*eta/(R*T)))  # A/cm^2

        # Flux boundary condition at x=0:  -D * dC/dx |0 = j / (nF)
        # Implement via ghost points or modify RHS (Neumann BC)
        # Here: finite difference with 2nd-order Neumann using ghost point -> C[-1] = C[1] to enforce dC/dx at boundary

        # Build RHS (B*C^k) and apply BCs
        rhs_O = B_O.dot(CO)
        rhs_R = B_R.dot(CR)

        # Neumann BC at x=0 from flux: (C1 - C0)/dx = - j/(nF*D)  -> modifies rhs
        rhs_O[0] += 2*lam_O * (CO[1] - CO[0])  # base CN
        rhs_R[0] += 2*lam_R * (CR[1] - CR[0])

        # Replace with flux term:
        rhs_O[0] += 2*lam_O * dx * ( - j/(n*F*D_O) )
        rhs_R[0] += 2*lam_R * dx * ( + j/(n*F*D_R) )

        # Dirichlet BC at x=L: CO(L)=C_bulk, CR(L)=0
        rhs_O[-1] = C_bulk
        rhs_R[-1] = 0.0
        A_O[-1,:] = 0.0; A_O[-1,-1] = 1.0
        A_R[-1,:] = 0.0; A_R[-1,-1] = 1.0

        # Solve for next concentrations
        CO = solve(A_O, rhs_O)
        CR = solve(A_R, rhs_R)
        # Prevent negatives from numerical noise
        CO = np.clip(CO, 0.0, None)
        CR = np.clip(CR, 0.0, None)

        # Current = j * A
        currents.append(j * A)

    currents = np.array(currents)  # A

    return E, currents, t

def demo_bq_hq_cv():
    # Example: use E°(pH=0)= +0.693 V vs SHE (literature for BQ/HQ), then shift to pH=14
    E0_pH0 = 0.693
    E0_pH14 = pourbaix_shift(E0_pH0, m=2, n=2, pH=14.0)  # ~ -0.13 V
    # Simulate a quasi-reversible 2e- process
    E, i, t = simulate_cv(E0=E0_pH14, n=2, alpha=0.5, k0=0.02,
                          D_O=8e-6, D_R=8e-6, C_bulk=2e-3,
                          A=0.07, scan_rate=0.1,
                          E_start=-0.4, E_vertex=0.4, L=0.03, Nx=400, dt=1e-3)
    plt.figure(figsize=(6,4))
    plt.plot(E, i*1e3, lw=2)  # mA
    plt.xlabel("Potential (V vs. SHE)")
    plt.ylabel("Current (mA)")
    plt.title("Simulated CV of BQ/HQ at pH 14 (quasi-reversible, n=2)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    demo_bq_hq_cv()
``
