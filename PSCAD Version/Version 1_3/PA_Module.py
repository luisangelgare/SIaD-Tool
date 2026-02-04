# -*- coding: utf-8 -*-
"""
PA_Module.py
-------------------------------------------------------------------------------
Passivity Assessment Tool for dq0 MIMO Systems
Technical University of Catalonia (UPC)
ETSEIB – CITCEA-UPC
Developed by: Luis Angel Garcia Reyes, MSc

This module evaluates passivity of two dq0 admittance systems by computing:
    - P(f) = Y(f) + Y(f)ᴴ
    - Minimum eigenvalue λ_min(f)
    - Passivity violation ranges
    - Summary table
    - High-quality PDF plot

GNU GPL v3.0 License
-------------------------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# =============================================================================
# PUBLIC ENTRY POINT
# =============================================================================

def Run_Passivity_Analysis(fd0, Ysys1, Ysys2, outcomes_dir):
    """
    Perform passivity assessment for dq0 admittance matrices.

    Parameters
    ----------
    fd0 : array-like
        Frequency vector (Hz).
    Ysys1 : dict
        Dictionary mapping freq -> 2×2 complex admittance matrix (System 1).
    Ysys2 : dict
        Dictionary mapping freq -> 2×2 complex admittance matrix (System 2).
    outcomes_dir : str
        Folder where results will be saved.

    Returns
    -------
    dict with:
        - lambda_min_sys1
        - lambda_min_sys2
        - non_passive_ranges_sys1
        - non_passive_ranges_sys2
    """

    fd0 = np.asarray(fd0)
    N = len(fd0)

    lambda_min_sys1 = np.zeros(N)
    lambda_min_sys2 = np.zeros(N)

    # -------------------------------------------------------------------------
    # 1. Compute minimum eigenvalues of P = Y + Yᴴ
    # -------------------------------------------------------------------------
    for k, f in enumerate(fd0):
        Y1 = Ysys1[f]
        Y2 = Ysys2[f]
    
        P1 = Y1 + Y1.conj().T
        P2 = Y2 + Y2.conj().T
    
        # NOTE: Passivity index definition adapted from Z-Tool (Cifuentes et al., CC BY-NC-ND 4.0)
        # PI(f) = 0.5 * min eigenvalue of (Y + Yᴴ)
        lambda_min_sys1[k] = 0.5 * np.min(np.linalg.eigvals(P1).real)
        lambda_min_sys2[k] = 0.5 * np.min(np.linalg.eigvals(P2).real)


    # -------------------------------------------------------------------------
    # 2. Detect non-passive ranges
    # -------------------------------------------------------------------------
    ranges1 = _detect_passivity_ranges(fd0, lambda_min_sys1, "System 1")
    ranges2 = _detect_passivity_ranges(fd0, lambda_min_sys2, "System 2")

    # -------------------------------------------------------------------------
    # 3. Plot and save PDF
    # -------------------------------------------------------------------------
    Path(outcomes_dir).mkdir(exist_ok=True)
    fig = _plot_passivity(fd0, lambda_min_sys1, lambda_min_sys2)

    pdf_path = Path(outcomes_dir) / "PA_Passivity.pdf"
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    plt.close(fig)

    # -------------------------------------------------------------------------
    # 4. Console summary
    # -------------------------------------------------------------------------
    print("---> Passivity analysis running...")
    print("\n---> Passivity Violation Summary:")
    print("---------------------------------------------")
    print("| System   | Frequency Range [Hz] | Status   |")
    print("---------------------------------------------")

    _print_ranges(ranges1, "System 1")
    _print_ranges(ranges2, "System 2")

    print("---------------------------------------------")
    print(f"---> Passivity plot saved to: {pdf_path}")

    return {
        "lambda_min_sys1": lambda_min_sys1,
        "lambda_min_sys2": lambda_min_sys2,
        "non_passive_ranges_sys1": ranges1,
        "non_passive_ranges_sys2": ranges2,
    }


# =============================================================================
# INTERNAL UTILITIES
# =============================================================================

def _detect_passivity_ranges(fd0, eigvals, system_name):
    """
    Identify frequency intervals where λ_min(f) < 0 (non-passive).
    """
    non_passive = eigvals < 0
    ranges = []
    start = None

    for i in range(len(fd0)):
        if non_passive[i] and start is None:
            start = fd0[i]
        elif not non_passive[i] and start is not None:
            ranges.append((start, fd0[i - 1]))
            start = None

    if start is not None:
        ranges.append((start, fd0[-1]))

    return ranges


def _print_ranges(ranges, system_name):
    """
    Print passivity violation ranges in table format.
    """
    if not ranges:
        print(f"| {system_name:<8} |     -     | Passive   |")
        return

    for r in ranges:
        print(f"| {system_name:<8} | {r[0]:5.1f} to {r[1]:5.1f} | Non-Passive |")


def _plot_passivity(fd0, lam1, lam2):
    """
    Create passivity plot showing λ_min(f) for both systems.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.grid(True)
    ax.set_xlabel("Frequency (Hz)", fontsize=12)
    ax.set_ylabel("Minimum Eigenvalue λ_min", fontsize=12)
    ax.set_title("Passivity Evaluation via Minimum Eigenvalues", fontsize=14)

    ax.plot(fd0, lam1, "-", color="b", linewidth=3, label="System 1")
    ax.plot(fd0, lam2, "-", color="g", linewidth=3, label="System 2")
    ax.plot(fd0, np.zeros_like(fd0), "--", color="r", linewidth=2, label="Passive Limit")

    ax.legend(fontsize=10)
    fig.tight_layout()

    return fig
