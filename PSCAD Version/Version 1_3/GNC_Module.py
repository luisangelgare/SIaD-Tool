# -*- coding: utf-8 -*-
# GNC_Module.py
# =============================================================================
#  Generalized Nyquist Criterion (GNC) Module for dq0 MIMO Systems
#  ---------------------------------------------------------------------------
#  Computes open-loop eigenvalues of L(s) = Ysys1(s) * Zsys2(s),
#  orders eigenvalue trajectories, plots MIMO Nyquist in the complex plane,
#  checks encirclements of the critical point (-1, 0), and saves the figure.
#
#  Intended to be called from SIaD-Tool after Ysys1, Ysys2 are computed.
# =============================================================================

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def Run_GNC(fd0, Ysys1, Ysys2, outcomes_dir, scanner_selector, show_critical_point=True):
    """
    Run Generalized Nyquist Criterion (GNC) analysis for dq0 MIMO systems.

    Parameters
    ----------
    fd0 : array-like
        Frequency vector (Hz) used in the identification (length N).
    Ysys1 : dict
        Dictionary mapping frequency -> 2x2 complex admittance matrix of System 1.
        Example: Ysys1[f] = np.array([[Yqq1, Yqd1], [Ydq1, Ydd1]])
    Ysys2 : dict
        Dictionary mapping frequency -> 2x2 complex admittance matrix of System 2.
    outcomes_dir : str, optional
        Folder where the Nyquist PDF will be stored. Default is "outcomes".
    show_critical_point : bool, optional
        If True, the critical point (-1, 0) is shown in the plot.

    Returns
    -------
    stability_info : dict
        Dictionary with basic stability assessment:
        {
            "eigenvalues": E,              # 2 x N complex array
            "encirclements_lambda1": n1,   # winding number around -1 for eigenvalue 1
            "encirclements_lambda2": n2,   # winding number around -1 for eigenvalue 2
            "is_stable": bool              # True if no encirclements, False otherwise
        }
    """

    fd0 = np.asarray(fd0)
    N = len(fd0)

    # -------------------------------------------------------------------------
    # 1. Build Y and Z matrices and compute open-loop L = Ysys1 * inv(Ysys2)
    # -------------------------------------------------------------------------
    L = np.zeros((2, 2, N), dtype=complex)
    E = np.zeros((2, N), dtype=complex)  # eigenvalues

    for idx, f in enumerate(fd0):
        Y1 = Ysys1[f]  # 2x2 complex
        Y2 = Ysys2[f]  # 2x2 complex

        Z2 = np.linalg.inv(Y2)
        L[:, :, idx] = Y1 @ Z2
        E[:, idx] = np.linalg.eigvals(L[:, :, idx])

    # -------------------------------------------------------------------------
    # 2. Order eigenvalues to obtain smooth trajectories (exchange algorithm)
    # -------------------------------------------------------------------------
    E_ordered = _order_eigenvalues(E)
    enc = _count_encirclements(E_ordered, fd0)
    
    # -------------------------------------------------------------------------
    # Nyquist Stability Margin (NSM) according to Garcia Reyes (2025)
    # -------------------------------------------------------------------------
    nsm_value, nsm_mode, nsm_idx = compute_nsm(E_ordered)
    nsm_freq = fd0[nsm_idx]

    # -------------------------------------------------------------------------
    # 3. Compute encirclements of the critical point (-1, 0) using winding number
    # -------------------------------------------------------------------------
    n1 = _winding_number(E_ordered[0, :], critical_point=-1 + 0j)
    n2 = _winding_number(E_ordered[1, :], critical_point=-1 + 0j)

    # Simple GNC-based decision:
    # If there are no encirclements of (-1,0), we declare "stable" (assuming open-loop stable).
    is_stable = (n1 == 0) and (n2 == 0)

    # -------------------------------------------------------------------------
    # 4. Plot MIMO Nyquist and save to PDF
    # -------------------------------------------------------------------------
    Path(outcomes_dir).mkdir(exist_ok=True)
    fig = _plot_mimo_nyquist(fd0, E_ordered, show_critical_point=show_critical_point)
    
    if scanner_selector == 2:
        pdf_path = Path(outcomes_dir) / "GNC_dq0_Nyquist.pdf"
        fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
        plt.close(fig)
    if scanner_selector == 3:
        pdf_path = Path(outcomes_dir) / "GNC_0pn_Nyquist.pdf"
        fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
        plt.close(fig)
        
    # -------------------------------------------------------------------------
    # 5. Console report
    # -------------------------------------------------------------------------
    print("---> GNC analysis completed.")
    print(f"---> Nyquist plot saved to: {pdf_path}")
    print(f"---> NSM (Nyquist Stability Margin) = {nsm_value:.6f}")
    print(f"     Achieved by eigenvalue λ_{nsm_mode} at f = {nsm_freq:.3f} Hz")
    # print("---> Encirclements around (-1,0):")
    # print("CW counts per eigenvalue:", enc["cw_counts"])
    # print("CCW counts per eigenvalue:", enc["ccw_counts"])
    # print("Net encirclements:", enc["net_encirclements"])
    print(f"---> Encirclements of (-1,0): λ1 -> {n1}, λ2 -> {n2}")
    print("---> Critical frequencies:", enc["critical_frequencies"])
    if is_stable:
        print("---> According to the Generalized Nyquist Criterion: STABLE response.")
    else:
        print("---> According to the Generalized Nyquist Criterion: UNSTABLE response.")

    return {
        "eigenvalues": E_ordered,
        "encirclements_lambda1": n1,
        "encirclements_lambda2": n2,
        "is_stable": is_stable,
    }


# =============================================================================
#  INTERNAL UTILITIES
# =============================================================================

def compute_nsm(E_sorted):
    """
    Compute the Nyquist Stability Margin (NSM) from ordered eigenvalues.

    Definition (Garcia Reyes, 2025):
        NSM = min_{ω ∈ Ω} ( min_i |1 + λ_i(jω)| )

    where λ_i(jω) are the eigenvalues of the loop transfer matrix L(jω)
    over the frequency range Ω.

    Parameters
    ----------
    E_sorted : np.ndarray
        dim x N array of ordered eigenvalues (each row = mode, each column = frequency).

    Returns
    -------
    nsm_value : float
        Nyquist Stability Margin (minimum distance to the critical point -1+0j).
    nsm_mode : int
        Index of the mode (1-based) where NSM occurs.
    nsm_index : int
        Index of the frequency sample (0-based) where NSM occurs.
    """
    # Shift eigenvalues by +1 to measure distance to (-1, 0)
    shifted = E_sorted + 1.0  # complex array
    distances = np.abs(shifted)  # dim x N

    # Global minimum over modes and frequencies
    idx_flat = np.argmin(distances)
    dim, N = E_sorted.shape
    mode_idx, freq_idx = divmod(idx_flat, N)

    nsm_value = float(distances[mode_idx, freq_idx])
    return nsm_value, int(mode_idx + 1), int(freq_idx)


def _count_encirclements(E_sorted, fd0):
    """
    Count CW and CCW encirclements of the critical point (-1, 0)
    using geometric crossing detection.

    NOTE: The geometric crossing logic implemented here is adapted
    from Z-Tool (Cifuentes et al., CC BY-NC-ND 4.0), specifically
    the method used in the 'nyquist' function to detect crossings
    of the imaginary axis and determine CW/CCW orientation.

    Parameters
    ----------
    E_sorted : np.ndarray
        2 x N array of ordered eigenvalues.
    fd0 : array-like
        Frequency vector.

    Returns
    -------
    dict with:
        - cw_counts
        - ccw_counts
        - net_encirclements
        - critical_frequencies
    """

    dim, N = E_sorted.shape
    cw_counts = [0] * dim
    ccw_counts = [0] * dim
    critical_freqs = []

    # Loop over each eigenvalue trajectory
    for idx in range(dim):

        lam = E_sorted[idx, :]
        x = lam.real
        y = lam.imag

        for k in range(1, N):

            # Check if the imaginary part crosses zero
            crossed = (y[k-1] < 0 and y[k] > 0) or (y[k-1] > 0 and y[k] < 0)
            if not crossed:
                continue

            # Only consider crossings to the LEFT of -1
            if x[k] > -1 and x[k-1] > -1:
                continue

            # Compute cross product to determine CW or CCW
            # NOTE: This cross-product orientation test is adapted
            # from Z-Tool's nyquist() implementation.
            v1 = np.array([x[k] - x[k-1], y[k] - y[k-1]])
            v2 = np.array([-1 - x[k-1], 0 - y[k-1]])
            cross = v1[0]*v2[1] - v1[1]*v2[0]

            # Determine CW or CCW
            if cross < 0:
                cw_counts[idx] += 1
                critical_freqs.append((fd0[k], "CW", idx+1))
            else:
                ccw_counts[idx] += 1
                critical_freqs.append((fd0[k], "CCW", idx+1))

    net = sum(cw_counts) - sum(ccw_counts)

    return {
        "cw_counts": cw_counts,
        "ccw_counts": ccw_counts,
        "net_encirclements": net,
        "critical_frequencies": critical_freqs
    }


def _order_eigenvalues(E):
    """
    Order eigenvalues column-wise to obtain smooth trajectories.

    Parameters
    ----------
    E : np.ndarray
        2 x N complex array of eigenvalues (each column corresponds to a frequency).

    Returns
    -------
    E_ordered : np.ndarray
        2 x N complex array with ordered eigenvalues.
    """
    E_ordered = E.copy()
    rows, cols = E_ordered.shape

    for j in range(1, cols):
        # Distance between current λ1 and previous λ1
        dist1 = _distance(E_ordered[0, j], E_ordered[0, j - 1])
        # Distance between current λ1 and previous λ2
        dist2 = _distance(E_ordered[0, j], E_ordered[1, j - 1])

        # If λ1(j) is closer to λ2(j-1), swap λ1 and λ2 at column j
        if dist1 > dist2:
            temp = E_ordered[0, j]
            E_ordered[0, j] = E_ordered[1, j]
            E_ordered[1, j] = temp

    return E_ordered


def _distance(z1, z2):
    """Euclidean distance between two complex numbers."""
    return abs(z1 - z2)


def _winding_number(eig_traj, critical_point=-1 + 0j):
    """
    Compute the winding number of an eigenvalue trajectory around a critical point.

    Parameters
    ----------
    eig_traj : array-like
        1D complex array representing the eigenvalue trajectory λ(f).
    critical_point : complex, optional
        Point around which the winding number is computed (default: -1+0j).

    Returns
    -------
    n : int
        Winding number (encirclements count, positive for CCW, negative for CW).
    """
    eig_traj = np.asarray(eig_traj)
    z = eig_traj - critical_point  # shift so that critical point is at origin

    angles = np.angle(z)
    angles_unwrapped = np.unwrap(angles)

    total_angle_change = angles_unwrapped[-1] - angles_unwrapped[0]
    n = int(np.round(total_angle_change / (2 * np.pi)))

    return n


def _plot_mimo_nyquist(frequencies, eigenvalues, show_critical_point=True):
    ordered_eigs = eigenvalues
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)

    ax.grid(True, which="both")
    ax.set_xlabel(r"Re($\lambda$)", fontsize=14)
    ax.set_ylabel(r"Im($\lambda$)", fontsize=14)

    # ============================================================
    # AUTO‑AJUSTE DE ESCALAS SEGÚN LOS DATOS + PUNTO CRÍTICO
    # ============================================================
    re_vals = np.concatenate([
        ordered_eigs[0, :].real,
        ordered_eigs[0, :].real,
        ordered_eigs[1, :].real,
        ordered_eigs[1, :].real
    ])
    im_vals = np.concatenate([
        ordered_eigs[0, :].imag,
        -ordered_eigs[0, :].imag,
        ordered_eigs[1, :].imag,
        -ordered_eigs[1, :].imag
    ])

    # Incluir también el punto crítico (-1, 0) en el cálculo de límites
    re_min_data, re_max_data = re_vals.min(), re_vals.max()
    im_min_data, im_max_data = im_vals.min(), im_vals.max()

    re_min = min(re_min_data, -1.0)
    re_max = max(re_max_data, -1.0)
    im_min = min(im_min_data, 0.0)
    im_max = max(im_max_data, 0.0)

    margin = 0.05
    dx = (re_max - re_min) * margin if re_max > re_min else 1.0
    dy = (im_max - im_min) * margin if im_max > im_min else 1.0

    ax.set_xlim([re_min - dx, re_max + dx])
    ax.set_ylim([im_min - dy, im_max + dy])
    # ============================================================

    # Main trajectories
    h1, = ax.plot(
        ordered_eigs[0, :].real,
        ordered_eigs[0, :].imag,
        "-",
        color="r",
        linewidth=2.5,
        label=r"$\lambda_{1}\ [0,+\infty]$",
    )
    h2, = ax.plot(
        ordered_eigs[0, :].real,
        -ordered_eigs[0, :].imag,
        "-",
        color="m",
        linewidth=2.5,
        label=r"$\lambda_{1}\ [-\infty,0]$",
    )
    h3, = ax.plot(
        ordered_eigs[1, :].real,
        ordered_eigs[1, :].imag,
        "-",
        color="b",
        linewidth=2.5,
        label=r"$\lambda_{2}\ [0,+\infty]$",
    )
    h4, = ax.plot(
        ordered_eigs[1, :].real,
        -ordered_eigs[1, :].imag,
        "-",
        color="k",
        linewidth=2.5,
        label=r"$\lambda_{2}\ [-\infty,0]$",
    )

    if show_critical_point:
        ax.plot(-1, 0, marker="p", markersize=7, color="g", linewidth=2.5)
        ax.text(-1, 0, "(-1, 0)", color="g", fontsize=10)

    # ============================================================
    # ZOOM (mantiene escalas fijas)
    # ============================================================
    X1, X2 = -1, 1
    Y1, Y2 = -0.1, 0.1

    ax_inset = fig.add_axes([0.65, 0.65, 0.25, 0.25])
    ax_inset.grid(True, which="both")

    ax_inset.plot(ordered_eigs[0, :].real,  ordered_eigs[0, :].imag,  "-", color="r", linewidth=1.5)
    ax_inset.plot(ordered_eigs[0, :].real, -ordered_eigs[0, :].imag,  "-", color="m", linewidth=1.5)
    ax_inset.plot(ordered_eigs[1, :].real,  ordered_eigs[1, :].imag,  "-", color="b", linewidth=1.5)
    ax_inset.plot(ordered_eigs[1, :].real, -ordered_eigs[1, :].imag,  "-", color="k", linewidth=1.5)

    if show_critical_point:
        ax_inset.plot(-1, 0, marker="p", markersize=5, color="g", linewidth=1.5)

    ax_inset.set_xlim([X1, X2])
    ax_inset.set_ylim([Y1, Y2])

    # Leyenda ordenada: primero ramas negativas, luego positivas
    ax.legend(handles=[h2, h1, h4, h3], loc="upper left", fontsize=10)

    return fig



