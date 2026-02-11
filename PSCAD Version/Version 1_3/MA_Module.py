# -*- coding: utf-8 -*-
"""
MA_Module.py
-------------------------------------------------------------------------------
Modal Analysis (MA) Tool for dq0 MIMO Systems
Technical University of Catalonia (UPC)
ETSEIB – CITCEA-UPC
Developed by: Luis Angel Garcia Reyes, MSc

This module performs modal analysis of the combined dq0 admittance:
    Y_full(f) = Ysys1(f) + Ysys2(f)

It computes:
    - Eigenvalues of Y(f)
    - Modal impedances (1 / λ)
    - Participation factors (left * right eigenvectors)
    - Critical mode identification
    - Dominant subsystem per frequency
    - Eigenvalue plot in the complex plane
    - Heatmap of participation factors
    - Summary table
    - PDF export to outcomes/

GNU GPL v3.0 License
-------------------------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# =============================================================================
# PUBLIC ENTRY POINT
# =============================================================================

def Run_Modal_Analysis(fd0, Ysys1, Ysys2, subsys_labels, outcomes_dir):
    fd0 = np.asarray(fd0)
    N = len(fd0)

    # -------------------------------------------------------------------------
    # 1. Build combined admittance Y_full(f) = Y1(f) + Y2(f)
    # -------------------------------------------------------------------------
    Y_full = np.zeros((2, 2, N), dtype=complex)

    for k, f in enumerate(fd0):
        Y_full[:, :, k] = Ysys1[f] + Ysys2[f]

    # -------------------------------------------------------------------------
    # 2. Allocate modal quantities
    # -------------------------------------------------------------------------
    eigenvalues = np.zeros((2, N), dtype=complex)
    modal_impedance = np.zeros((2, N))
    participation_factors = np.zeros((2, 2, N))
    critical_mode_idx = np.zeros(N, dtype=int)
    dominant_subsystem = np.empty(N, dtype=object)
    Z_modes = np.zeros((2, N), dtype=complex)  # eigenvalores modales de impedancia

    # -------------------------------------------------------------------------
    # 3. Figura: Re(Z_modo) vs frecuencia (MA_Eigenvalues.pdf)
    # -------------------------------------------------------------------------
    Path(outcomes_dir).mkdir(exist_ok=True)
    fig1, ax1 = plt.subplots(figsize=(8, 6))

    ax1.grid(True, which="both", linestyle="--", linewidth=0.7)
    ax1.set_xscale("log")
    ax1.set_xlabel("f (Hz)")
    ax1.set_ylabel(r"Re($Z_{\mathrm{mode}}$)")
    ax1.set_title("Real part of modal impedances vs frequency")
    ax1.axhline(0, linestyle="--", color="r", linewidth=2, label="Re(Z) = 0 boundary")

    # -------------------------------------------------------------------------
    # 4. Modal analysis loop
    # -------------------------------------------------------------------------
    for i in range(N):
        Yf = Y_full[:, :, i]

        eigvals, V = np.linalg.eig(Yf)   # Right eigenvectors
        W = np.linalg.inv(V)             # Left eigenvectors

        eigenvalues[:, i] = eigvals

        # Modal impedance (compleja)
        Z_modes[:, i] = 1.0 / eigvals
        modal_impedance[:, i] = np.abs(Z_modes[:, i])

        # Participation factors = |W' .* V|
        participation_factors[:, :, i] = np.abs(W.T * V)

        # Critical mode = smallest |λ|
        idx = np.argmin(np.abs(eigvals))
        critical_mode_idx[i] = idx

        # Dominant subsystem
        pf = participation_factors[idx, :, i]
        if pf[0] > pf[1]:
            dominant_subsystem[i] = subsys_labels[0]
        elif pf[1] > pf[0]:
            dominant_subsystem[i] = subsys_labels[1]
        else:
            dominant_subsystem[i] = "Undetermined"

    # Plot de Re(Z_modo) vs frecuencia (dos modos)
    ax1.plot(fd0, Z_modes[0, :].real, "-", color="b", linewidth=2.0, label=r"Mode 1")
    ax1.plot(fd0, Z_modes[1, :].real, "-", color="g", linewidth=2.0, label=r"Mode 2")
    ax1.legend(loc="best", fontsize=9)
    fig1.tight_layout()

    pdf1 = Path(outcomes_dir) / "MA_Eigenvalues.pdf"
    fig1.savefig(pdf1, format="pdf", bbox_inches="tight")
    plt.close(fig1)

    # -------------------------------------------------------------------------
    # 5. Participation factors: curvas normalizadas vs frecuencia
    #     (MA_ParticipationFactors.pdf)
    # -------------------------------------------------------------------------
    pf_norm = np.zeros((N, 2))
    for i in range(N):
        pf = participation_factors[critical_mode_idx[i], :, i]
        max_pf = np.max(pf)
        pf_norm[i, :] = pf / max_pf if max_pf > 0 else 0.0

    fig2, ax2 = plt.subplots(figsize=(8, 6))
    ax2.grid(True, which="both", linestyle="--", linewidth=0.7)

    ax2.set_xscale("log")
    ax2.set_xlabel("f (Hz)")
    ax2.set_ylabel("Normalized participation factor")
    ax2.set_title("Participation factors of critical mode vs frequency")

    h_pf1, = ax2.plot(fd0, pf_norm[:, 0], "-", color="b", linewidth=2.0, label="System 1")
    h_pf2, = ax2.plot(fd0, pf_norm[:, 1], "-", color="m", linewidth=2.0, label="System 2")

    ax2.legend(handles=[h_pf1, h_pf2], loc="best", fontsize=9)
    fig2.tight_layout()

    pdf2 = Path(outcomes_dir) / "MA_ParticipationFactors.pdf"
    fig2.savefig(pdf2, format="pdf", bbox_inches="tight")
    plt.close(fig2)

    # -------------------------------------------------------------------------
    # 6. Console summary table
    # -------------------------------------------------------------------------
    print("---> Impedance Modal Analysis results:")
    print("---")
    print("\nModal Summary by Frequency:")
    print("-" * 98)
    print("|   f (Hz)  | Crit |      Eigenvalue       |  |Z_m|   |   PF1   |   PF2   | Dominant     |")
    print("-" * 98)

    for i in range(N):
        idx = critical_mode_idx[i]
        lam = eigenvalues[idx, i]
        Zm = modal_impedance[idx, i]
        pf = participation_factors[idx, :, i]

        print(f"| {fd0[i]:9.3f} |  {idx+1}   | "
              f"{lam.real:8.4f}{lam.imag:+8.4f}j | "
              f"{Zm:7.3f} | "
              f"{pf[0]:7.4f} | {pf[1]:7.4f} | "
              f"{dominant_subsystem[i]:12} |")

    print("-" * 98)
    print(f"---> Eigenvalue/impedance plot saved to: {pdf1}")
    print(f"---> Participation factors plot saved to: {pdf2}")
    print("---> Modal analysis completed.\n")

    return {
        "eigenvalues": eigenvalues,
        "modal_impedance": modal_impedance,
        "participation_factors": participation_factors,
        "critical_mode_idx": critical_mode_idx,
        "dominant_subsystem": dominant_subsystem,
    }
