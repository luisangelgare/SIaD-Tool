# -*- coding: utf-8 -*-
"""
PM_Module.py
-------------------------------------------------------------------------------
Phase Margin (PM) computation for dq0 MIMO systems.
Technical University of Catalonia (UPC)
ETSEIB – CITCEA-UPC
Developed by: Luis Angel Garcia Reyes, MSc

This module computes:
    - Crossing frequency (|Yc/Yg| = 1)
    - Phase margin PM = 180° + angle(Yc/Yg)
for each element of the dq0 admittance matrix.

GNU GPL v3.0 License
-------------------------------------------------------------------------------
"""

import numpy as np
from pathlib import Path
import pandas as pd


# =============================================================================
# PUBLIC ENTRY POINT
# =============================================================================

def Run_Phase_Margins(fd0, Ysys1, Ysys2, scanner_selector):
    """
    Compute phase margins for dq0 admittance matrices.

    Parameters
    ----------
    fd0 : array-like
        Frequency vector (Hz).
    Ysys1 : dict
        Dictionary mapping freq -> 2×2 complex converter admittance.
    Ysys2 : dict
        Dictionary mapping freq -> 2×2 complex grid admittance.

    Returns
    -------
    pandas.DataFrame
        Table with crossing frequencies and phase margins.
    """

    fd0 = np.asarray(fd0)
    N = len(fd0)

    # Build 2×2×N arrays
    Yc = np.zeros((2, 2, N), dtype=complex)
    Yg = np.zeros((2, 2, N), dtype=complex)

    for k, f in enumerate(fd0):
        Yc[:, :, k] = Ysys1[f]
        Yg[:, :, k] = Ysys2[f]

    # Compute PM for each element
    cf_dd, pm_dd = _PMs(fd0, Yc[0, 0, :], Yg[0, 0, :])
    cf_dq, pm_dq = _PMs(fd0, Yc[0, 1, :], Yg[0, 1, :])
    cf_qd, pm_qd = _PMs(fd0, Yc[1, 0, :], Yg[1, 0, :])
    cf_qq, pm_qq = _PMs(fd0, Yc[1, 1, :], Yg[1, 1, :])

    if scanner_selector == 2:
        cases = ["qq", "qd", "dq", "dd"]
        crossing = [cf_dd, cf_dq, cf_qd, cf_qq]
        margins = [pm_dd, pm_dq, pm_qd, pm_qq]
    if scanner_selector == 3:
        cases = ["pp", "pn", "np", "nn"]
        crossing = [cf_dd, cf_dq, cf_qd, cf_qq]
        margins = [pm_dd, pm_dq, pm_qd, pm_qq]

    df = pd.DataFrame({
        "Case": cases,
        "Crossing_Frequency [Hz]": crossing,
        "PM [deg]": margins
    })

    print("\n---> Phase Margin Results (dq0)")
    print(df)
    print("----------------------------------\n")

    return df


# =============================================================================
# INTERNAL: PM COMPUTATION
# =============================================================================

def _PMs(fd0, Yc_vec, Yg_vec):
    """
    Compute crossing frequency and phase margin for a single dq element.

    Parameters
    ----------
    fd0 : array-like
        Frequency vector.
    Yc_vec : array-like
        Converter admittance element Yc(i,j,f).
    Yg_vec : array-like
        Grid admittance element Yg(i,j,f).

    Returns
    -------
    crossing_freq : float or NaN
    PM_deg : float or NaN
    """

    L = Yc_vec / Yg_vec  # Open-loop element

    mag = np.abs(L)
    phase = np.angle(L, deg=True)

    # Find |L| = 1 crossing
    idx = np.where(np.diff(np.sign(mag - 1)) != 0)[0]

    if len(idx) == 0:
        return np.nan, np.nan

    k = idx[0]

    # Linear interpolation of crossing frequency
    f1, f2 = fd0[k], fd0[k + 1]
    m1, m2 = mag[k], mag[k + 1]
    crossing_freq = f1 + (1 - m1) * (f2 - f1) / (m2 - m1)

    # Phase margin at crossing
    ph1, ph2 = phase[k], phase[k + 1]
    phase_cross = ph1 + (crossing_freq - f1) * (ph2 - ph1) / (f2 - f1)

    PM = 180 + phase_cross  # Correct sign convention

    return crossing_freq, PM
