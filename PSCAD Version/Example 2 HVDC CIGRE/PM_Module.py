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

# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def Run_Phase_Margins(fd0, Ysys1, Ysys2, scanner_selector, outcomes_dir,
                      pm_threshold_deg=30.0, gm_threshold_db=6.0):
    """
    Compute phase and gain margins for dq0 or pn admittance matrices.

    Includes:
        - PM (Phase Margin)
        - GM (Gain Margin)
        - Crossing frequencies
        - Alert column
        - PDF with graphical interpretation of PM & GM

    Parameters
    ----------
    fd0 : array-like
        Frequency vector (Hz).
    Ysys1 : dict
        Converter admittance matrices.
    Ysys2 : dict
        Grid admittance matrices.
    scanner_selector : int
        2 = dq0, 3 = pn.
    pm_threshold_deg : float
        Minimum acceptable PM.
    gm_threshold_db : float
        Minimum acceptable GM.
    outcomes_dir : str
        Folder to save PDF.

    Returns
    -------
    pandas.DataFrame
    """

    fd0 = np.asarray(fd0)
    N = len(fd0)

    # Build 2×2×N arrays
    Yc = np.zeros((2, 2, N), dtype=complex)
    Yg = np.zeros((2, 2, N), dtype=complex)

    for k, f in enumerate(fd0):
        Yc[:, :, k] = Ysys1[f]
        Yg[:, :, k] = Ysys2[f]

    # Open-loop elements L = Yc / Yg
    L_dd = Yc[0, 0, :] / Yg[0, 0, :]
    L_dq = Yc[0, 1, :] / Yg[0, 1, :]
    L_qd = Yc[1, 0, :] / Yg[1, 0, :]
    L_qq = Yc[1, 1, :] / Yg[1, 1, :]

    # Compute PM & GM for each element
    cf_dd, pm_dd, gm_dd = _margins_single_loop(fd0, L_dd)
    cf_dq, pm_dq, gm_dq = _margins_single_loop(fd0, L_dq)
    cf_qd, pm_qd, gm_qd = _margins_single_loop(fd0, L_qd)
    cf_qq, pm_qq, gm_qq = _margins_single_loop(fd0, L_qq)

    if scanner_selector == 2:
        cases = ["qq", "qd", "dq", "dd"]
    elif scanner_selector == 3:
        cases = ["pp", "pn", "np", "nn"]
    else:
        raise ValueError("scanner_selector must be 2 (dq0) or 3 (pn).")

    crossing = [cf_dd, cf_dq, cf_qd, cf_qq]
    pm_list  = [pm_dd, pm_dq, pm_qd, pm_qq]
    gm_list  = [gm_dd, gm_dq, gm_qd, gm_qq]

    # Alert column
    alerts = []
    for pm, gm in zip(pm_list, gm_list):
        if np.isnan(pm) and np.isnan(gm):
            alerts.append("No margin (no crossover)")
        elif np.isnan(pm):
            alerts.append("No PM (no |L|=1 crossing)")
        elif np.isnan(gm):
            alerts.append("No GM (no phase=-180° crossing)")
        elif (pm < pm_threshold_deg) or (gm < gm_threshold_db):
            alerts.append("WARNING: Margin below threshold")
        else:
            alerts.append("OK")

    df = pd.DataFrame({
        "Case": cases,
        "Crossing_Frequency_PM [Hz]": crossing,
        "PM [deg]": pm_list,
        "GM [dB]": gm_list,
        "Alert": alerts
    })

    # Replace NaN with readable text
    df["Crossing_Frequency_PM [Hz]"] = df["Crossing_Frequency_PM [Hz]"].apply(
        lambda x: "No gain crossover" if pd.isna(x) else round(x, 4)
    )
    df["PM [deg]"] = df["PM [deg]"].apply(
        lambda x: "Not defined (no |L|=1)" if pd.isna(x) else round(x, 2)
    )
    df["GM [dB]"] = df["GM [dB]"].apply(
        lambda x: "Not defined (no phase=-180°)" if pd.isna(x) else round(x, 2)
    )

    print("\n---> Phase & Gain Margin Results (dq0 / pn)")
    print(df)
    print("----------------------------------\n")

    # ---------------------------------------------------------------------
    # PDF plotting
    # ---------------------------------------------------------------------
    Path(outcomes_dir).mkdir(exist_ok=True)
    pdf_path = Path(outcomes_dir) / "PM_GM_Analysis.pdf"

    _plot_margins_bode(
        fd0,
        [L_dd, L_dq, L_qd, L_qq],
        cases,
        [cf_dd, cf_dq, cf_qd, cf_qq],
        [pm_dd, pm_dq, pm_qd, pm_qq],
        [gm_dd, gm_dq, gm_qd, gm_qq],
        pdf_path
    )

    print(f"---> Phase/Gain margin plot saved to: {pdf_path}")

    return df


# =============================================================================
# INTERNAL: PM & GM COMPUTATION
# =============================================================================

def _margins_single_loop(fd0, L_vec):
    """
    Compute PM and GM for a single open-loop element L(jw).
    """

    fd0 = np.asarray(fd0)
    L_vec = np.asarray(L_vec)

    mag = np.abs(L_vec)
    phase = np.angle(L_vec, deg=True)
    phase_unwrapped = np.unwrap(np.deg2rad(phase)) * 180.0 / np.pi

    # -------------------------
    # 1) Phase margin (gain crossover: |L| = 1)
    # -------------------------
    idx_mag = np.where(np.diff(np.sign(mag - 1.0)) != 0)[0]

    if len(idx_mag) == 0:
        crossing_freq_pm = np.nan
        PM_deg = np.nan
    else:
        k = idx_mag[0]
        f1, f2 = fd0[k], fd0[k + 1]
        m1, m2 = mag[k], mag[k + 1]

        crossing_freq_pm = f1 + (1.0 - m1) * (f2 - f1) / (m2 - m1)

        ph1, ph2 = phase_unwrapped[k], phase_unwrapped[k + 1]
        phase_cross = ph1 + (crossing_freq_pm - f1) * (ph2 - ph1) / (f2 - f1)

        PM_deg = 180.0 + phase_cross

    # -------------------------
    # 2) Gain margin (phase crossover: arg(L) = -180 deg)
    # -------------------------
    target_phase = -180.0
    phase_diff = phase_unwrapped - target_phase
    idx_phase = np.where(np.diff(np.sign(phase_diff)) != 0)[0]

    if len(idx_phase) == 0:
        GM_dB = np.nan
    else:
        k = idx_phase[0]
        f1, f2 = fd0[k], fd0[k + 1]
        ph1, ph2 = phase_unwrapped[k], phase_unwrapped[k + 1]

        crossing_freq_gm = f1 + (target_phase - ph1) * (f2 - f1) / (ph2 - ph1)

        m1, m2 = mag[k], mag[k + 1]
        mag_cross = m1 + (crossing_freq_gm - f1) * (m2 - m1) / (f2 - f1)

        GM_dB = -20.0 * np.log10(mag_cross)

    return crossing_freq_pm, PM_deg, GM_dB


# =============================================================================
# INTERNAL: PLOTTING
# =============================================================================

def _plot_margins_bode(fd0, L_list, case_labels,
                       cf_list, pm_list, gm_list, pdf_path):
    """
    Plot magnitude and phase with markers for PM and GM.
    Includes visual indicators when no crossing exists.
    """

    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(10, 9), sharex=True)

    colors = ["r", "g", "b", "k"]

    for L_vec, label, cf, pm, gm, color in zip(
            L_list, case_labels, cf_list, pm_list, gm_list, colors):

        mag = np.abs(L_vec)
        phase = np.angle(L_vec, deg=True)
        phase_unwrapped = np.unwrap(np.deg2rad(phase)) * 180.0 / np.pi

        # Magnitude
        ax_mag.semilogx(fd0, 20 * np.log10(mag), color=color,
                        linewidth=2.0, label=f"{label} |L|")

        # Phase
        ax_phase.semilogx(fd0, phase_unwrapped, color=color,
                          linewidth=2.0, label=f"{label} ∠L")

        # PM marker
        if not np.isnan(cf):
            mag_interp = np.interp(cf, fd0, mag)
            ax_mag.semilogx(cf, 20 * np.log10(mag_interp),
                            marker="o", color=color, markersize=7)
            ax_mag.text(cf, 20*np.log10(mag_interp),
                        f" PM={pm:.1f}°", color=color)

        else:
            ax_mag.text(fd0[len(fd0)//2], -10,
                        f"{label}: No PM (no |L|=1 crossing)",
                        color=color)

        # GM marker
        if not np.isnan(gm):
            # Find phase crossing again
            target_phase = -180.0
            phase_diff = phase_unwrapped - target_phase
            idx_phase = np.where(np.diff(np.sign(phase_diff)) != 0)[0]

            if len(idx_phase) > 0:
                k = idx_phase[0]
                f1, f2 = fd0[k], fd0[k + 1]
                ph1, ph2 = phase_unwrapped[k], phase_unwrapped[k + 1]
                crossing_freq_gm = f1 + (target_phase - ph1) * (f2 - f1) / (ph2 - ph1)

                m1, m2 = mag[k], mag[k + 1]
                mag_cross = m1 + (crossing_freq_gm - f1) * (m2 - m1) / (f2 - f1)

                ax_mag.semilogx(crossing_freq_gm, 20*np.log10(mag_cross),
                                marker="s", color=color, markersize=7)
                ax_mag.text(crossing_freq_gm, 20*np.log10(mag_cross),
                            f" GM={gm:.1f} dB", color=color)
        else:
            ax_mag.text(fd0[len(fd0)//2], -20,
                        f"{label}: No GM (no phase=-180° crossing)",
                        color=color)

    # Formatting
    ax_mag.grid(True, which="both", linestyle="--", alpha=0.5)
    ax_phase.grid(True, which="both", linestyle="--", alpha=0.5)

    ax_mag.set_ylabel("Magnitude [dB]")
    ax_phase.set_ylabel("Phase [deg]")
    ax_phase.set_xlabel("Frequency [Hz] (log scale)")

    ax_mag.legend()
    ax_phase.legend()

    fig.tight_layout()
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    plt.close(fig)

    return fig