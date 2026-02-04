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
    """
    Perform dq0 modal analysis for the combined admittance Y = Y1 + Y2.

    Parameters
    ----------
    fd0 : array-like
        Frequency vector (Hz).
    Ysys1 : dict
        Dictionary mapping freq -> 2×2 complex matrix.
    Ysys2 : dict
        Dictionary mapping freq -> 2×2 complex matrix.
    subsys_labels : tuple of str
        Names of the two subsystems (default: ("Sys1", "Sys2")).
    outcomes_dir : str
        Folder where results will be saved.

    Returns
    -------
    dict with:
        - eigenvalues
        - modal_impedance
        - participation_factors
        - critical_mode_idx
        - dominant_subsystem
    """

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

    # -------------------------------------------------------------------------
    # 3. Plot eigenvalues in complex plane
    # -------------------------------------------------------------------------
    Path(outcomes_dir).mkdir(exist_ok=True)
    fig1, ax1 = plt.subplots(figsize=(8, 6))

    ax1.grid(True)
    ax1.set_xlabel("Re(λ)")
    ax1.set_ylabel("Im(λ)")
    ax1.set_title("Eigenvalues of Y(f) across Frequencies")
    ax1.axvline(0, linestyle="--", color="r", linewidth=2, label="Re(λ)=0 boundary")

    # -------------------------------------------------------------------------
    # 4. Modal analysis loop
    # -------------------------------------------------------------------------
    for i in range(N):
        Yf = Y_full[:, :, i]

        # Right eigenvectors V and eigenvalues D
        eigvals, V = np.linalg.eig(Yf)
        W = np.linalg.inv(V)  # Left eigenvectors

        eigenvalues[:, i] = eigvals

        # Modal impedance = |1 / λ|
        modal_impedance[:, i] = np.abs(1.0 / eigvals)

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

        # Plot eigenvalues
        ax1.plot(eigvals.real, eigvals.imag, "o", markersize=7,
                 label=f"f = {fd0[i]:.2f} Hz" if i == 0 else None)

    ax1.legend(loc="best", fontsize=8)
    fig1.tight_layout()

    # Save eigenvalue plot
    pdf1 = Path(outcomes_dir) / "MA_Eigenvalues.pdf"
    fig1.savefig(pdf1, format="pdf", bbox_inches="tight")
    plt.close(fig1)

    # -------------------------------------------------------------------------
    # 5. Heatmap of participation factors for critical mode
    # -------------------------------------------------------------------------
    heatmap_data = np.zeros((N, 2))
    row_labels = [f"f = {f:.2f}" for f in fd0]

    for i in range(N):
        heatmap_data[i, :] = participation_factors[critical_mode_idx[i], :, i]

    fig2, ax2 = plt.subplots(figsize=(7, 8))
    im = ax2.imshow(heatmap_data, cmap="viridis", aspect="auto")

    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(["Input 1", "Input 2"])
    ax2.set_yticks(range(N))
    ax2.set_yticklabels(row_labels)

    ax2.set_title("Participation Factors of Critical Mode")
    fig2.colorbar(im, ax=ax2)

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
    print(f"---> Eigenvalue plot saved to: {pdf1}")
    print(f"---> Participation heatmap saved to: {pdf2}")
    print("---> Modal analysis completed.\n")

    return {
        "eigenvalues": eigenvalues,
        "modal_impedance": modal_impedance,
        "participation_factors": participation_factors,
        "critical_mode_idx": critical_mode_idx,
        "dominant_subsystem": dominant_subsystem,
    }

