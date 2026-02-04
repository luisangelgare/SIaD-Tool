# -*- coding: utf-8 -*-
# =============================================================================
#  SIaD-Tool: Stability & Interaction assessment in the frequency-Domain Tool
#  ---------------------------------------------------------------------------
#  Technical University of Catalonia (UPC)
#  ETSEIB – CITCEA-UPC
#  Developed by: Luis Angel Garcia Reyes, MSc
#
#  Description:
#      Comprehensive frequency-domain tool for small-signal stability
#      and interaction assessment in modern power systems.
#
#  License:
#      CC-BY-NC-ND 4.0
#      Attribution-NonCommercial-NoDerivatives 4.0 International 
#      Copyright (C) 2025 Luis Angel Garcia Reyes, UPC-MSCA-ADOreD
#      Email: luis.reyes@upc.edu
#      This work is licensed under the Creative Commons 
#      Attribution‑NonCommercial‑NoDerivatives 4.0 International License.
#      You are free to copy and redistribute this material in any medium or format,
#      provided that you give appropriate credit to the original author.
#
#      You may NOT use this material for commercial purposes.
#      You may NOT distribute modified versions of this material; only verbatim
#      copies are permitted under this license.
#
#      This work is provided “as is”, without any warranty of any kind, either 
#      expressed or implied, including but not limited to the warranties of 
#      merchantability, fitness for a particular purpose, or non‑infringement.
#
#      A copy of the full license is available at:
#      https://creativecommons.org/licenses/by-nc-nd/4.0/


#      This work has received funding from the ADOreD project
#      under the European Union’s Horizon Europe Research and 
#      Innovation Programme under the Marie Skłodowska-Curie 
#      Grant Agreement No. 101073554.
# =============================================================================

# FDScanning.py
#
# Frequency-domain scanning routine for SIaD-Tool.
# Runs steady-state initialization, snapshot, and frequency-response
# identification in ABC, dq0, or 0pn frames using PSCAD.
#
# Outputs:
#   - fd0   : frequency vector used for scanning
#   - Ysys1 : dictionary {f: Y-matrix} for system 1
#   - Ysys2 : dictionary {f: Y-matrix} for system 2

import os
from decimal import Decimal
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mhi.pscad

from Utilities import (
    project_file,
    disable_except,
    run_steady_state,
    run_perturbation,
    find_all_out_files,
    find_out_file,
    save_YsysABC_txt,
    save_Ysysqd_txt,
    save_Ysyspn_txt,
    ABC_plot,
    ABC_plot_response,
    dq0_plot_response,
    pn0_plot_response,
    dq0_plot,
    pn0_plot,
    cartz2pol,
    choose,
    parse_complex
)


def FDScanningTool(
    workspace_name,
    project_name,
    step_time,
    fs,
    f_max,
    f_points,
    Tinit,
    f_scale,
    ss_snap,
    Vbase,
    Sbase,
    f0,
    scanner_selector,
    scanner_type,
    signal_type,
    Vperturbation,
    fortran_ext,
    out_filename,
    results_folder,
):
    """
    FDScanning
    ----------
    Run frequency-domain scanning in PSCAD for SIaD-Tool.

    The function:
      1) Initializes PSCAD and time-domain parameters.
      2) Runs a steady-state simulation and stores a snapshot.
      3) Runs frequency-response identification in the selected frame:
         - scanner_selector = 1 -> ABC frame
         - scanner_selector = 2 -> dq0 frame
         - scanner_selector = 3 -> 0pn frame
      4) Builds Y-matrices for two systems (Ysys1, Ysys2).
      5) Saves results (txt + pdf) into `results_folder`.

    Parameters
    ----------
    workspace_name : str
        Name of the PSCAD workspace.
    project_name : str
        Name of the PSCAD project.
    step_time : float
        Simulation time step.
    fs : float
        Fundamental frequency (Hz) used as lower bound for fd0.
    f_max : float
        Maximum frequency (Hz) for the scan.
    f_points : int
        Number of frequency points in the scan.
    Tinit : float
        Initial time before observation window (s).
    Vbase : float
        Base voltage (line-to-line RMS).
    Sbase : float
        Base apparent power.
    f0 : float
        Base frequency (Hz).
    scanner_selector : int
        1 = ABC scanner, 2 = dq0 scanner, 3 = 0pn scanner.
    scanner_type : int
        Perturbation type selector (used only for printing).
    signal_type : int
        Signal type selector (used only for printing).
    Vperturbation : float
        Relative perturbation magnitude (per unit of steady-state voltage).
    fortran_ext : str
        Fortran extension used by PSCAD output files.
    out_filename : str
        Base name of PSCAD output files.
    results_folder : str
        Folder where results (txt + pdf) will be stored.

    Returns
    -------
    fd0 : np.ndarray
        Frequency vector used for scanning.
    Ysys1 : dict
        Dictionary mapping frequency -> Y-matrix (system 1).
    Ysys2 : dict
        Dictionary mapping frequency -> Y-matrix (system 2).
    """

    # =========================================================================
    #  INITIALIZATION & PSCAD SETUP
    # =========================================================================

    print("---")
    print("---> SIaD-Tool is running...")

    project_file_path = project_file(workspace_name)
    print("---> Loaded project:", project_file_path)

    # Time-domain parameters
    coef_decimal = Decimal("1E-6")
    delta_t = float(coef_decimal * Decimal(step_time))
    fsampling = np.round(1 / delta_t)

    Tobs0 = 1 / fs
    t_w1 = Tinit
    t_w2 = Tinit + Tobs0
    Tobs = t_w2
    samples = int(Tobs / delta_t)
    plot_step = step_time * 1

    # Frequency vector for scanning
    if f_scale == 1:
        # Logarithmic spacing
        raw = np.logspace(np.log10(fs), np.log10(f_max), f_points)
    
    elif f_scale == 2:
        # Linear spacing
        raw = np.linspace(fs, f_max, f_points)
    
    # Snap all values to nearest multiple of fs
    fd0 = np.round(raw / fs) * fs
    
    # Remove duplicates and enforce bounds
    fd0 = np.unique(fd0)
    fd0 = fd0[(fd0 >= fs) & (fd0 <= f_max)]


    # Steady-state parameters (initial guesses)
    dist_act = 1000 * 1
    Vdist = 0.0
    R_source = 1e-6 * 1
    Vpeak = (Vbase / np.sqrt(3)) * np.sqrt(2) * 1
    Ipeak = Sbase / Vpeak

    # Dictionary to store steady-state signals
    steady_state_signals = {}

    # Launch PSCAD and load project
    pscad = mhi.pscad.launch()
    pscad.load(project_file_path)
    project = pscad.project(project_name)

    # Set simulation parameters for steady-state run
    project.parameters(
        {
            "time_duration": Tobs,
            "time_step": step_time,
            "sample_step": plot_step,
        }
    )

    main = project.canvas("Main")
    
    # Create or retrieve Simulation Set for steady-state
    if "Steady_State" not in pscad.simulation_sets():
        simset = pscad.create_simulation_set("Steady_State")
        simset.add_tasks(project_name)
    else:
        simset = pscad.simulation_set("Steady_State")

    task = simset.tasks()[0]

    # Scanner timing for steady-state
    t_coupling = (Tobs + 10) * 1
    t_decoupling = (Tobs + 10) * 1
    t_ABC = t_decoupling * 1
    t_dq0 = t_decoupling * 1
    t_0pn = t_decoupling * 1

    # Search for the main Disturbance_Source scanner component
    base_def = f"{project_name}:Disturbance_Source"
    scanner_list = project.find_all(base_def)

    if not scanner_list:
        for i in range(1, 50):
            definition = f"{base_def}_{i}"
            scanner_list = project.find_all(definition)
            if scanner_list:
                break

    if not scanner_list:
        raise ValueError(
            f"No Disturbance_Source found for project '{project_name}'."
        )

    Scanner = scanner_list[0]

    # Configure scanner for steady-state run
    Scanner.parameters(
        Vbase=Vbase,
        Vq_ss=0,
        Vd_ss=0,
        fbase=f0,
        Vdist=Vdist,
        Rsource=R_source,
        dist_time=dist_act,
        theta_ss=0,
        t_coupling=t_coupling,
        t_decoupling=t_decoupling,
        t_abc=t_ABC,
        t_dq0=t_dq0,
        t_0pn=t_0pn,
    )

    # Search for ABC disturbance block
    base_def1 = f"{project_name}:ABC_Disturbance"
    scanner_list1 = project.find_all(base_def1)

    if not scanner_list1:
        for k in range(1, 50):
            definition1 = f"{base_def1}_{k}"
            scanner_list1 = project.find_all(definition1)
            if scanner_list1:
                break

    if not scanner_list1:
        raise ValueError(
            "No ABC_Disturbance block found in PSCAD project."
        )

    abc_block = scanner_list1[0]
    
    # Search for dq0 disturbance block
    base_def2 = f"{project_name}:dq0_Disturbance_1"
    scanner_list2 = project.find_all(base_def2)

    if not scanner_list2:
        for i in range(1, 50):
            definition2 = f"{base_def2}_{i}"
            scanner_list2 = project.find_all(definition2)
            if scanner_list2:
                break

    if not scanner_list2:
        raise ValueError(
            f"No dq0_Disturbance_1 found for project '{project_name}'."
        )

    dq0_block = scanner_list2[0]
    
    # Search for pn0 disturbance block
    base_def3 = f"{project_name}:pn0_Disturbance"
    scanner_list3 = project.find_all(base_def3)

    if not scanner_list3:
        for k in range(1, 50):
            definition3 = f"{base_def3}_{k}"
            scanner_list3 = project.find_all(definition3)
            if scanner_list3:
                break

    if not scanner_list3:
        raise ValueError(
            "No pn0_Disturbance block found in PSCAD project."
        )

    pn0_block = scanner_list3[0]

    # Configure layers and measurement blocks for steady-state
    if scanner_selector == 1:
        
        # Setting project components in Disable or Enable mode
        abc_block.enable(True)
        dq0_block.enable(False)
        pn0_block.enable(False)
        
        outputs_meas = {"Vqd2", "theta", "Vabc1", "Iabc1", "Vabc2", "Iabc2"}
        disable_except(project.find_all("master:pgb"), outputs_meas)
        for meter in project.find_all("master:multimeter"):
            meter.parameters(Dis=0)

    if scanner_selector == 2:
        
        # Setting project components in Disable or Enable mode
        abc_block.enable(False)
        dq0_block.enable(True)
        pn0_block.enable(False)
        
        outputs_meas = {"Vqd1", "Iqd1", "Vqd2", "Iqd2", "theta"}
        disable_except(project.find_all("master:pgb"), outputs_meas)
        for meter in project.find_all("master:multimeter"):
            meter.parameters(Dis=0)

    if scanner_selector == 3:

        # Setting project components in Disable or Enable mode
        abc_block.enable(False)
        dq0_block.enable(False)
        pn0_block.enable(True)
        
        outputs_meas = {"Vqd2", "theta", "V0pn1", "I0pn1", "V0pn2", "I0pn2"}
        disable_except(project.find_all("master:pgb"), outputs_meas)
        for meter in project.find_all("master:multimeter"):
            meter.parameters(Dis=0)
    
    
    # =========================================================================
    #  STAGE 1 — STEADY-STATE SIMULATION + SNAPSHOT
    # =========================================================================

    print("---> SIaD-Tool is obtaining the steady state...")

    snapshot_name = "SIaD_steady_snap"

    snapshot_file = run_steady_state(
        simset=simset,
        task=task,
        snapshot_name=snapshot_name,
        Tinit=Tinit,
        Tobs0=Tobs0,
        project_name=project_name,
        fortran_ext=fortran_ext,
        out_filename=out_filename,
        scanner_selector=scanner_selector,
        ss_snap=ss_snap,
        steady_state_store=steady_state_signals,
    )
    
    # =========================================================================
    #  STAGE 2 — FREQUENCY SCAN (PERTURBATION SIMULATIONS)
    # =========================================================================

    print("---> SIaD-Tool is running the system identification process...")
    print(
        f"---> Simulation set with "
        f"{choose(scanner_type, 'None', 'voltage', 'current')} perturbation, "
        f"{choose(signal_type, 'None', 'single-tone', 'PRBS', 'multi-tone')} signal, "
        f"{choose(scanner_selector, 'None', 'ABC', 'dq0', '0pn')} frame."
    )

    # Create or retrieve Simulation Set for perturbations
    if "Perturbation" not in pscad.simulation_sets():
        simset = pscad.create_simulation_set("Perturbation")
        simset.add_tasks(project_name)
    else:
        simset = pscad.simulation_set("Perturbation")

    task = simset.tasks()[0]

    # Scanner timing for perturbation runs
    t_coupling = 0.0 * 1
    t_decoupling = 0.0 * 1
    t_ABC = t_decoupling * 1
    t_dq0 = t_decoupling * 1
    t_0pn = t_decoupling * 1

    # Steady-state values from previous stage (rounded)
    Vq_ss = round(steady_state_signals["vq_ss"] * 1, 4)
    Vd_ss = round(steady_state_signals["vd_ss"] * 1, 4)
    theta_ss = round(steady_state_signals["theta_ss"] * 1, 4)
    Vmag_ss = np.sqrt(Vq_ss**2 + Vd_ss**2) * 1
    Vdist = round(Vperturbation * Vmag_ss * 1, 4)
    dist_act = 0.0 * 1
    Tdelay = 0.5 * 1

    # Update scanner parameters for perturbation runs
    Scanner.parameters(
        Vbase=Vbase,
        Vq_ss=Vq_ss,
        Vd_ss=Vd_ss,
        fbase=f0,
        Vdist=Vdist,
        Rsource=R_source,
        dist_time=dist_act,
        theta_ss=theta_ss,
        t_coupling=t_coupling,
        t_decoupling=t_decoupling,
        t_abc=t_ABC,
        t_dq0=t_dq0,
        t_0pn=t_0pn,
    )
    
    if ss_snap == 1:
        # Update project time for perturbation runs with snapshot
        project.parameters(
            {
                "time_duration": Tobs0 + Tdelay,
                "time_step": step_time,
                "sample_step": plot_step,
            }
        )

    # Dictionaries for Y-matrices and polar representation
    Ysys1 = {}
    Ysys2 = {}
    Y_mag1 = {}
    Y_rad1 = {}
    Y_mag2 = {}
    Y_rad2 = {}

    # Ensure results folder exists
    Path(results_folder).mkdir(exist_ok=True)

    # -------------------------------------------------------------------------
    #  ABC SCANNER — FREQUENCY RESPONSE EXTRACTION
    # -------------------------------------------------------------------------
    if scanner_selector == 1:
        # Configure measurement blocks for ABC identification
        outputs_meas = {"Vabc1", "Iabc1", "Vabc2", "Iabc2"}
        disable_except(project.find_all("master:pgb"), outputs_meas)
        for meter in project.find_all("master:multimeter"):
            meter.parameters(Dis=0)

        t_ABC = t_decoupling
        t_dq0 = 10000
        t_0pn = 10000

        # Output file addresses (two .out files)
        steady_out_file1 = find_all_out_files(
            project_name, fortran_ext, out_filename
        )

        for i in range(len(fd0)):
            # -------------------------------------------------------------
            # Injection in a-axis
            # -------------------------------------------------------------
            Scanner.parameters(
                Vdist=Vdist,
                fd=fd0[i],
                t_abc=t_ABC,
                t_dq0=t_dq0,
                t_0pn=t_0pn,
            )
            abc_block.parameters(
                dist_time=dist_act,
                fd=fd0[i],
                Vdist_a=Vdist,
                Vdist_b=0.0,
                Vdist_c=0.0,
            )

            run_perturbation(
                simset,
                task,
                ss_snap,
                snapshot_name,
                project_name,
                fortran_ext,
                out_filename,
            )

            # Load .out files
            df1 = pd.read_csv(
                steady_out_file1[0],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned1 = df1.dropna(axis=1, how="all")

            df2 = pd.read_csv(
                steady_out_file1[1],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned2 = df2.dropna(axis=1, how="all")

            # Assign column names
            df_cleaned1.columns = [
                "time",
                "v_a1",
                "v_b1",
                "v_c1",
                "i_a1",
                "i_b1",
                "i_c1",
                "v_a2",
                "v_b2",
                "v_c2",
                "i_a2",
            ]
            df_cleaned2.columns = [
                "time",
                "i_b2",
                "i_c2",
            ]

            # Set index and extract windows
            df_cleaned1 = df_cleaned1.set_index("time")
            df_cleaned2 = df_cleaned2.set_index("time")
            
            if ss_snap == 1:
                time_w1 = df_cleaned1.loc[Tdelay : Tobs0 + Tdelay]
                time_w2 = df_cleaned2.loc[Tdelay : Tobs0 + Tdelay]
                
            if ss_snap == 0:
                time_w1 = df_cleaned1.loc[Tinit:Tinit + Tobs]
                time_w2 = df_cleaned2.loc[Tinit:Tinit + Tobs]

            # Subtract steady-state values
            v_a1 = (
                time_w1["v_a1"].to_numpy()
                - steady_state_signals["v_a1_ss"]
            )
            v_b1 = (
                time_w1["v_b1"].to_numpy()
                - steady_state_signals["v_b1_ss"]
            )
            v_c1 = (
                time_w1["v_c1"].to_numpy()
                - steady_state_signals["v_c1_ss"]
            )
            i_a1 = (
                time_w1["i_a1"].to_numpy()
                - steady_state_signals["i_a1_ss"]
            )
            i_b1 = (
                time_w1["i_b1"].to_numpy()
                - steady_state_signals["i_b1_ss"]
            )
            i_c1 = (
                time_w1["i_c1"].to_numpy()
                - steady_state_signals["i_c1_ss"]
            )
            v_a2 = (
                time_w1["v_a2"].to_numpy()
                - steady_state_signals["v_a2_ss"]
            )
            v_b2 = (
                time_w1["v_b2"].to_numpy()
                - steady_state_signals["v_b2_ss"]
            )
            v_c2 = (
                time_w1["v_c2"].to_numpy()
                - steady_state_signals["v_c2_ss"]
            )
            i_a2 = (
                time_w1["i_a2"].to_numpy()
                - steady_state_signals["i_a2_ss"]
            )
            i_b2 = (
                time_w2["i_b2"].to_numpy()
                - steady_state_signals["i_b2_ss"]
            )
            i_c2 = (
                time_w2["i_c2"].to_numpy()
                - steady_state_signals["i_c2_ss"]
            )

            # FFT
            Va1 = np.fft.fft(v_a1) / len(v_a1)
            Vb1 = np.fft.fft(v_b1) / len(v_b1)
            Vc1 = np.fft.fft(v_c1) / len(v_c1)
            Ia1 = np.fft.fft(i_a1) / len(i_a1)
            Ib1 = np.fft.fft(i_b1) / len(i_b1)
            Ic1 = np.fft.fft(i_c1) / len(i_c1)

            Va2 = np.fft.fft(v_a2) / len(v_a2)
            Vb2 = np.fft.fft(v_b2) / len(v_b2)
            Vc2 = np.fft.fft(v_c2) / len(v_c2)
            Ia2 = np.fft.fft(i_a2) / len(i_a2)
            Ib2 = np.fft.fft(i_b2) / len(i_b2)
            Ic2 = np.fft.fft(i_c2) / len(i_c2)

            frequencies = np.fft.fftfreq(len(v_a1), d=Tobs0)
            positive_freqs = frequencies[frequencies >= 0]

            Va1 = Va1[: len(positive_freqs)]
            Vb1 = Vb1[: len(positive_freqs)]
            Vc1 = Vc1[: len(positive_freqs)]
            Ia1 = Ia1[: len(positive_freqs)]
            Ib1 = Ib1[: len(positive_freqs)]
            Ic1 = Ic1[: len(positive_freqs)]

            Va2 = Va2[: len(positive_freqs)]
            Vb2 = Vb2[: len(positive_freqs)]
            Vc2 = Vc2[: len(positive_freqs)]
            Ia2 = Ia2[: len(positive_freqs)]
            Ib2 = Ib2[: len(positive_freqs)]
            Ic2 = Ic2[: len(positive_freqs)]

            wd = round(fd0[i] / fs)

            Yaa1 = Ia1[wd] / Va1[wd]
            Yba1 = Ib1[wd] / Va1[wd]
            Yca1 = Ic1[wd] / Va1[wd]
            Yaa2 = Ia2[wd] / Va2[wd]
            Yba2 = Ib2[wd] / Va2[wd]
            Yca2 = Ic2[wd] / Va2[wd]

            # -------------------------------------------------------------
            # Injection in b-axis
            # -------------------------------------------------------------
            abc_block.parameters(
                Vdist_a=0.0, Vdist_b=Vdist, Vdist_c=0.0
            )

            run_perturbation(
                simset,
                task,
                ss_snap,
                snapshot_name,
                project_name,
                fortran_ext,
                out_filename,
            )

            df1 = pd.read_csv(
                steady_out_file1[0],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned1 = df1.dropna(axis=1, how="all")

            df2 = pd.read_csv(
                steady_out_file1[1],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned2 = df2.dropna(axis=1, how="all")

            df_cleaned1.columns = [
                "time",
                "v_a1",
                "v_b1",
                "v_c1",
                "i_a1",
                "i_b1",
                "i_c1",
                "v_a2",
                "v_b2",
                "v_c2",
                "i_a2",
            ]
            df_cleaned2.columns = [
                "time",
                "i_b2",
                "i_c2",
            ]

            df_cleaned1 = df_cleaned1.set_index("time")
            df_cleaned2 = df_cleaned2.set_index("time")

            if ss_snap == 1:
                time_w1 = df_cleaned1.loc[Tdelay : Tobs0 + Tdelay]
                time_w2 = df_cleaned2.loc[Tdelay : Tobs0 + Tdelay]
                
            if ss_snap == 0:
                time_w1 = df_cleaned1.loc[Tinit:Tinit + Tobs]
                time_w2 = df_cleaned2.loc[Tinit:Tinit + Tobs]

            v_a1 = (
                time_w1["v_a1"].to_numpy()
                - steady_state_signals["v_a1_ss"]
            )
            v_b1 = (
                time_w1["v_b1"].to_numpy()
                - steady_state_signals["v_b1_ss"]
            )
            v_c1 = (
                time_w1["v_c1"].to_numpy()
                - steady_state_signals["v_c1_ss"]
            )
            i_a1 = (
                time_w1["i_a1"].to_numpy()
                - steady_state_signals["i_a1_ss"]
            )
            i_b1 = (
                time_w1["i_b1"].to_numpy()
                - steady_state_signals["i_b1_ss"]
            )
            i_c1 = (
                time_w1["i_c1"].to_numpy()
                - steady_state_signals["i_c1_ss"]
            )
            v_a2 = (
                time_w1["v_a2"].to_numpy()
                - steady_state_signals["v_a2_ss"]
            )
            v_b2 = (
                time_w1["v_b2"].to_numpy()
                - steady_state_signals["v_b2_ss"]
            )
            v_c2 = (
                time_w1["v_c2"].to_numpy()
                - steady_state_signals["v_c2_ss"]
            )
            i_a2 = (
                time_w1["i_a2"].to_numpy()
                - steady_state_signals["i_a2_ss"]
            )
            i_b2 = (
                time_w2["i_b2"].to_numpy()
                - steady_state_signals["i_b2_ss"]
            )
            i_c2 = (
                time_w2["i_c2"].to_numpy()
                - steady_state_signals["i_c2_ss"]
            )

            Va1 = np.fft.fft(v_a1) / len(v_a1)
            Vb1 = np.fft.fft(v_b1) / len(v_b1)
            Vc1 = np.fft.fft(v_c1) / len(v_c1)
            Ia1 = np.fft.fft(i_a1) / len(i_a1)
            Ib1 = np.fft.fft(i_b1) / len(i_b1)
            Ic1 = np.fft.fft(i_c1) / len(i_c1)

            Va2 = np.fft.fft(v_a2) / len(v_a2)
            Vb2 = np.fft.fft(v_b2) / len(v_b2)
            Vc2 = np.fft.fft(v_c2) / len(v_c2)
            Ia2 = np.fft.fft(i_a2) / len(i_a2)
            Ib2 = np.fft.fft(i_b2) / len(i_b2)
            Ic2 = np.fft.fft(i_c2) / len(i_c2)

            Va1 = Va1[: len(positive_freqs)]
            Vb1 = Vb1[: len(positive_freqs)]
            Vc1 = Vc1[: len(positive_freqs)]
            Ia1 = Ia1[: len(positive_freqs)]
            Ib1 = Ib1[: len(positive_freqs)]
            Ic1 = Ic1[: len(positive_freqs)]

            Va2 = Va2[: len(positive_freqs)]
            Vb2 = Vb2[: len(positive_freqs)]
            Vc2 = Vc2[: len(positive_freqs)]
            Ia2 = Ia2[: len(positive_freqs)]
            Ib2 = Ib2[: len(positive_freqs)]
            Ic2 = Ic2[: len(positive_freqs)]

            Yab1 = Ia1[wd] / Vb1[wd]
            Ybb1 = Ib1[wd] / Vb1[wd]
            Ycb1 = Ic1[wd] / Vb1[wd]
            Yab2 = Ia2[wd] / Vb2[wd]
            Ybb2 = Ib2[wd] / Vb2[wd]
            Ycb2 = Ic2[wd] / Vb2[wd]

            # -------------------------------------------------------------
            # Injection in c-axis
            # -------------------------------------------------------------
            abc_block.parameters(
                Vdist_a=0.0, Vdist_b=0.0, Vdist_c=Vdist
            )

            run_perturbation(
                simset,
                task,
                ss_snap,
                snapshot_name,
                project_name,
                fortran_ext,
                out_filename,
            )

            df1 = pd.read_csv(
                steady_out_file1[0],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned1 = df1.dropna(axis=1, how="all")

            df2 = pd.read_csv(
                steady_out_file1[1],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned2 = df2.dropna(axis=1, how="all")

            df_cleaned1.columns = [
                "time",
                "v_a1",
                "v_b1",
                "v_c1",
                "i_a1",
                "i_b1",
                "i_c1",
                "v_a2",
                "v_b2",
                "v_c2",
                "i_a2",
            ]
            df_cleaned2.columns = [
                "time",
                "i_b2",
                "i_c2",
            ]

            df_cleaned1 = df_cleaned1.set_index("time")
            df_cleaned2 = df_cleaned2.set_index("time")

            if ss_snap == 1:
                time_w1 = df_cleaned1.loc[Tdelay : Tobs0 + Tdelay]
                time_w2 = df_cleaned2.loc[Tdelay : Tobs0 + Tdelay]
                
            if ss_snap == 0:
                time_w1 = df_cleaned1.loc[Tinit:Tinit + Tobs]
                time_w2 = df_cleaned2.loc[Tinit:Tinit + Tobs]

            v_a1 = (
                time_w1["v_a1"].to_numpy()
                - steady_state_signals["v_a1_ss"]
            )
            v_b1 = (
                time_w1["v_b1"].to_numpy()
                - steady_state_signals["v_b1_ss"]
            )
            v_c1 = (
                time_w1["v_c1"].to_numpy()
                - steady_state_signals["v_c1_ss"]
            )
            i_a1 = (
                time_w1["i_a1"].to_numpy()
                - steady_state_signals["i_a1_ss"]
            )
            i_b1 = (
                time_w1["i_b1"].to_numpy()
                - steady_state_signals["i_b1_ss"]
            )
            i_c1 = (
                time_w1["i_c1"].to_numpy()
                - steady_state_signals["i_c1_ss"]
            )
            v_a2 = (
                time_w1["v_a2"].to_numpy()
                - steady_state_signals["v_a2_ss"]
            )
            v_b2 = (
                time_w1["v_b2"].to_numpy()
                - steady_state_signals["v_b2_ss"]
            )
            v_c2 = (
                time_w1["v_c2"].to_numpy()
                - steady_state_signals["v_c2_ss"]
            )
            i_a2 = (
                time_w1["i_a2"].to_numpy()
                - steady_state_signals["i_a2_ss"]
            )
            i_b2 = (
                time_w2["i_b2"].to_numpy()
                - steady_state_signals["i_b2_ss"]
            )
            i_c2 = (
                time_w2["i_c2"].to_numpy()
                - steady_state_signals["i_c2_ss"]
            )

            Va1 = np.fft.fft(v_a1) / len(v_a1)
            Vb1 = np.fft.fft(v_b1) / len(v_b1)
            Vc1 = np.fft.fft(v_c1) / len(v_c1)
            Ia1 = np.fft.fft(i_a1) / len(i_a1)
            Ib1 = np.fft.fft(i_b1) / len(i_b1)
            Ic1 = np.fft.fft(i_c1) / len(i_c1)

            Va2 = np.fft.fft(v_a2) / len(v_a2)
            Vb2 = np.fft.fft(v_b2) / len(v_b2)
            Vc2 = np.fft.fft(v_c2) / len(v_c2)
            Ia2 = np.fft.fft(i_a2) / len(i_a2)
            Ib2 = np.fft.fft(i_b2) / len(i_b2)
            Ic2 = np.fft.fft(i_c2) / len(i_c2)

            Va1 = Va1[: len(positive_freqs)]
            Vb1 = Vb1[: len(positive_freqs)]
            Vc1 = Vc1[: len(positive_freqs)]
            Ia1 = Ia1[: len(positive_freqs)]
            Ib1 = Ib1[: len(positive_freqs)]
            Ic1 = Ic1[: len(positive_freqs)]

            Va2 = Va2[: len(positive_freqs)]
            Vb2 = Vb2[: len(positive_freqs)]
            Vc2 = Vc2[: len(positive_freqs)]
            Ia2 = Ia2[: len(positive_freqs)]
            Ib2 = Ib2[: len(positive_freqs)]
            Ic2 = Ic2[: len(positive_freqs)]

            Yac1 = Ia1[wd] / Vc1[wd]
            Ybc1 = Ib1[wd] / Vc1[wd]
            Ycc1 = Ic1[wd] / Vc1[wd]
            Yac2 = Ia2[wd] / Vc2[wd]
            Ybc2 = Ib2[wd] / Vc2[wd]
            Ycc2 = Ic2[wd] / Vc2[wd]

            # Assemble 3x3 Y-matrices for both systems
            Ysys10 = np.array(
                [
                    [Yaa1, Yab1, Yac1],
                    [Yba1, Ybb1, Ybc1],
                    [Yca1, Ycb1, Ycc1],
                ]
            )
            Ysys20 = np.array(
                [
                    [Yaa2, Yab2, Yac2],
                    [Yba2, Ybb2, Ybc2],
                    [Yca2, Ycb2, Ycc2],
                ]
            )

            Ysys1[fd0[i]] = Ysys10
            Ysys2[fd0[i]] = Ysys20

            Y_mag1[fd0[i]], Y_rad1[fd0[i]] = cartz2pol(Ysys10)
            Y_mag2[fd0[i]], Y_rad2[fd0[i]] = cartz2pol(Ysys20)

        # End of ABC scanning
        pscad.quit()
        
        # Save results in results_folder using existing utilities
        Path(results_folder).mkdir(exist_ok=True)
        
        # Save results in results_folder using existing utilities
        save_YsysABC_txt("Ysys1", Ysys1, fd0, results_folder)
        save_YsysABC_txt("Ysys2", Ysys2, fd0, results_folder)

        print("---> SIaD-Tool has finished.")
        print(
            "---> Results stored in the files: Ysys1.txt and Ysys2.txt."
        )

        fig1 = ABC_plot(fd0, Y_mag1, Y_rad1)
        fig1.savefig(f"{results_folder}/Ysys1.pdf", format="pdf", bbox_inches="tight")
        plt.close(fig1)

        fig2 = ABC_plot(fd0, Y_mag2, Y_rad2)
        fig2.savefig(f"{results_folder}/Ysys2.pdf", format="pdf", bbox_inches="tight")
        plt.close(fig2)

        return fd0, Ysys1, Ysys2

    # -------------------------------------------------------------------------
    #  dq0 SCANNER — FREQUENCY RESPONSE EXTRACTION
    # -------------------------------------------------------------------------
    if scanner_selector == 2:
        # Configure measurement blocks for dq0 identification
        outputs_meas = {"Vqd1", "Iqd1", "Vqd2", "Iqd2"}
        disable_except(project.find_all("master:pgb"), outputs_meas)
        for meter in project.find_all("master:multimeter"):
            meter.parameters(Dis=0)

        t_ABC = 10000 * 1
        t_dq0 = t_decoupling * 1
        t_0pn = 10000 * 1

        # dq0_block = project.component(953818110)

        # Single .out file for dq0
        steady_out_file1 = find_out_file(
            project_name, fortran_ext, out_filename
        )

        for i in range(len(fd0)):
            # -------------------------------------------------------------
            # Injection in q-axis
            # -------------------------------------------------------------
            Scanner.parameters(
                Vdist=Vdist,
                fd=fd0[i],
                t_abc=t_ABC,
                t_dq0=t_dq0,
                t_0pn=t_0pn,
            )
            dq0_block.parameters(Vdist_q=Vdist, Vdist_d=0)

            run_perturbation(
                simset,
                task,
                ss_snap,
                snapshot_name,
                project_name,
                fortran_ext,
                out_filename,
            )
            
            df = pd.read_csv(
                steady_out_file1,
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned = df.dropna(axis=1, how="all")

            df_cleaned.columns = [
                "time",
                "i_q1",
                "i_d1",
                "v_q1",
                "v_d1",
                "i_q2",
                "i_d2",
                "v_q2",
                "v_d2",
            ]
            df_cleaned = df_cleaned.set_index("time")

            if ss_snap == 1: 
                time_w = df_cleaned.loc[Tdelay : Tobs0 + Tdelay]
            
            if ss_snap == 0: 
                time_w = df_cleaned.loc[Tinit:Tinit + Tobs]

            v_q1 = (
                time_w["v_q1"].to_numpy()
                - steady_state_signals["v_q1_ss"]
            )
            v_d1 = (
                time_w["v_d1"].to_numpy()
                - steady_state_signals["v_d1_ss"]
            )
            i_q1 = (
                time_w["i_q1"].to_numpy()
                - steady_state_signals["i_q1_ss"]
            )
            i_d1 = (
                time_w["i_d1"].to_numpy()
                - steady_state_signals["i_d1_ss"]
            )
            v_q2 = (
                time_w["v_q2"].to_numpy()
                - steady_state_signals["v_q2_ss"]
            )
            v_d2 = (
                time_w["v_d2"].to_numpy()
                - steady_state_signals["v_d2_ss"]
            )
            i_q2 = (
                time_w["i_q2"].to_numpy()
                - steady_state_signals["i_q2_ss"]
            )
            i_d2 = (
                time_w["i_d2"].to_numpy()
                - steady_state_signals["i_d2_ss"]
            )

            Vq1 = np.fft.fft(v_q1) / len(v_q1)
            Vd1 = np.fft.fft(v_d1) / len(v_d1)
            Iq1 = np.fft.fft(i_q1) / len(i_q1)
            Id1 = np.fft.fft(i_d1) / len(i_d1)
            Vq2 = np.fft.fft(v_q2) / len(v_q2)
            Vd2 = np.fft.fft(v_d2) / len(v_d2)
            Iq2 = np.fft.fft(i_q2) / len(i_q2)
            Id2 = np.fft.fft(i_d2) / len(i_d2)

            frequencies = np.fft.fftfreq(len(v_q1), d=Tobs0)
            positive_freqs = frequencies[frequencies >= 0]

            Vq1 = Vq1[: len(positive_freqs)]
            Vd1 = Vd1[: len(positive_freqs)]
            Iq1 = Iq1[: len(positive_freqs)]
            Id1 = Id1[: len(positive_freqs)]
            Vq2 = Vq2[: len(positive_freqs)]
            Vd2 = Vd2[: len(positive_freqs)]
            Iq2 = Iq2[: len(positive_freqs)]
            Id2 = Id2[: len(positive_freqs)]

            wd = round(fd0[i] / fs)

            Yqq1 = Iq1[wd] / Vq1[wd]
            Ydq1 = Id1[wd] / Vq1[wd]
            Yqq2 = Iq2[wd] / Vq2[wd]
            Ydq2 = Id2[wd] / Vq2[wd]

            # -------------------------------------------------------------
            # Injection in d-axis
            # -------------------------------------------------------------
            dq0_block.parameters(Vdist_q=0, Vdist_d=Vdist)

            run_perturbation(
                simset,
                task,
                ss_snap,
                snapshot_name,
                project_name,
                fortran_ext,
                out_filename,
            )

            df = pd.read_csv(
                steady_out_file1,
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned = df.dropna(axis=1, how="all")

            df_cleaned.columns = [
                "time",
                "i_q1",
                "i_d1",
                "v_q1",
                "v_d1",
                "i_q2",
                "i_d2",
                "v_q2",
                "v_d2",
            ]
            df_cleaned = df_cleaned.set_index("time")
            
            if ss_snap == 1: 
                time_w = df_cleaned.loc[Tdelay : Tobs0 + Tdelay]
            
            if ss_snap == 0: 
                time_w = df_cleaned.loc[Tinit:Tinit + Tobs]

            v_q1 = (
                time_w["v_q1"].to_numpy()
                - steady_state_signals["v_q1_ss"]
            )
            v_d1 = (
                time_w["v_d1"].to_numpy()
                - steady_state_signals["v_d1_ss"]
            )
            i_q1 = (
                time_w["i_q1"].to_numpy()
                - steady_state_signals["i_q1_ss"]
            )
            i_d1 = (
                time_w["i_d1"].to_numpy()
                - steady_state_signals["i_d1_ss"]
            )
            v_q2 = (
                time_w["v_q2"].to_numpy()
                - steady_state_signals["v_q2_ss"]
            )
            v_d2 = (
                time_w["v_d2"].to_numpy()
                - steady_state_signals["v_d2_ss"]
            )
            i_q2 = (
                time_w["i_q2"].to_numpy()
                - steady_state_signals["i_q2_ss"]
            )
            i_d2 = (
                time_w["i_d2"].to_numpy()
                - steady_state_signals["i_d2_ss"]
            )

            Vq1 = np.fft.fft(v_q1) / len(v_q1)
            Vd1 = np.fft.fft(v_d1) / len(v_d1)
            Iq1 = np.fft.fft(i_q1) / len(i_q1)
            Id1 = np.fft.fft(i_d1) / len(i_d1)
            Vq2 = np.fft.fft(v_q2) / len(v_q2)
            Vd2 = np.fft.fft(v_d2) / len(v_d2)
            Iq2 = np.fft.fft(i_q2) / len(i_q2)
            Id2 = np.fft.fft(i_d2) / len(i_d2)

            Vq1 = Vq1[: len(positive_freqs)]
            Vd1 = Vd1[: len(positive_freqs)]
            Iq1 = Iq1[: len(positive_freqs)]
            Id1 = Id1[: len(positive_freqs)]
            Vq2 = Vq2[: len(positive_freqs)]
            Vd2 = Vd2[: len(positive_freqs)]
            Iq2 = Iq2[: len(positive_freqs)]
            Id2 = Id2[: len(positive_freqs)]

            Yqd1 = Iq1[wd] / Vd1[wd]
            Ydd1 = Id1[wd] / Vd1[wd]
            Yqd2 = Iq2[wd] / Vd2[wd]
            Ydd2 = Id2[wd] / Vd2[wd]

            # Assemble 2x2 dq-admittance matrices
            Ysys10 = np.array([[Yqq1, Yqd1], [Ydq1, Ydd1]])
            Ysys20 = np.array([[Yqq2, Yqd2], [Ydq2, Ydd2]])

            Ysys1[fd0[i]] = Ysys10
            Ysys2[fd0[i]] = Ysys20

            Y_mag1[fd0[i]], Y_rad1[fd0[i]] = cartz2pol(Ysys10)
            Y_mag2[fd0[i]], Y_rad2[fd0[i]] = cartz2pol(Ysys20)

        # End of dq0 scanning
        pscad.quit()

        # Save results in results_folder using existing utilities
        Path(results_folder).mkdir(exist_ok=True)

        save_Ysysqd_txt("Ysys1", Ysys1, fd0, results_folder)
        save_Ysysqd_txt("Ysys2", Ysys2, fd0, results_folder)

        print("---> SIaD-Tool has finished.")
        print(
            "---> Results stored in the files: Ysys1.txt and Ysys2.txt."
        )

        fig1 = dq0_plot(fd0, Y_mag1, Y_rad1)
        fig1.savefig(f"{results_folder}/Ysys1.pdf", format="pdf", bbox_inches="tight")
        plt.close(fig1)

        fig2 = dq0_plot(fd0, Y_mag2, Y_rad2)
        fig2.savefig(f"{results_folder}/Ysys2.pdf", format="pdf", bbox_inches="tight")
        plt.close(fig2)

        return fd0, Ysys1, Ysys2

    # -------------------------------------------------------------------------
    #  0pn SCANNER — FREQUENCY RESPONSE EXTRACTION
    # -------------------------------------------------------------------------
    if scanner_selector == 3:
        
        # Configure measurement blocks for 0pn identification
        outputs_meas = {"V0pn1", "I0pn1", "V0pn2", "I0pn2"}
        disable_except(project.find_all("master:pgb"), outputs_meas)
        for meter in project.find_all("master:multimeter"):
            meter.parameters(Dis=0)

        t_ABC = 10000
        t_dq0 = 10000
        t_0pn = t_decoupling


        # Output file addresses (three .out files)
        steady_out_file1 = find_all_out_files(
            project_name, fortran_ext, out_filename
        )

        for i in range(len(fd0)):
            # -------------------------------------------------------------
            # Injection in p-axis
            # -------------------------------------------------------------
            Scanner.parameters(
                Vdist=Vdist,
                fd=fd0[i],
                t_abc=t_ABC,
                t_dq0=t_dq0,
                t_0pn=t_0pn,
            )
            pn0_block.parameters(
                Vdist_p=Vdist, Vdist_n=0.0, Vdist_0=0.0
            )

            run_perturbation(
                simset,
                task,
                ss_snap,
                snapshot_name,
                project_name,
                fortran_ext,
                out_filename,
            )

            # Load .out files
            df1 = pd.read_csv(
                steady_out_file1[0],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned1 = df1.dropna(axis=1, how="all")

            df2 = pd.read_csv(
                steady_out_file1[1],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned2 = df2.dropna(axis=1, how="all")

            df3 = pd.read_csv(
                steady_out_file1[2],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned3 = df3.dropna(axis=1, how="all")

            # Assign column names
            df_cleaned1.columns = [
                "time",
                "R_v_01",
                "R_v_p1",
                "R_v_n1",
                "I_v_01",
                "I_v_p1",
                "I_v_n1",
                "R_i_01",
                "R_i_p1",
                "R_i_n1",
                "I_i_01",
            ]

            df_cleaned2.columns = [
                "time",
                "I_i_p1",
                "I_i_n1",
                "R_v_02",
                "R_v_p2",
                "R_v_n2",
                "I_v_02",
                "I_v_p2",
                "I_v_n2",
                "R_i_02",
                "R_i_p2",
            ]

            df_cleaned3.columns = [
                "time",
                "R_i_n2",
                "I_i_02",
                "I_i_p2",
                "I_i_n2",
            ]

            # Set index
            df_cleaned1 = df_cleaned1.set_index("time")
            df_cleaned2 = df_cleaned2.set_index("time")
            df_cleaned3 = df_cleaned3.set_index("time")

            # System 1: complex signals
            df_cleaned1["i_01"] = (
                df_cleaned1["R_i_01"] + 1j * df_cleaned1["I_i_01"]
            )
            df_cleaned1["i_n1"] = (
                df_cleaned1["R_i_n1"] + 1j * df_cleaned2["I_i_n1"]
            )
            df_cleaned1["i_p1"] = (
                df_cleaned1["R_i_p1"] + 1j * df_cleaned2["I_i_p1"]
            )

            df_cleaned1["v_01"] = (
                df_cleaned1["R_v_01"] + 1j * df_cleaned1["I_v_01"]
            )
            df_cleaned1["v_n1"] = (
                df_cleaned1["R_v_n1"] + 1j * df_cleaned1["I_v_n1"]
            )
            df_cleaned1["v_p1"] = (
                df_cleaned1["R_v_p1"] + 1j * df_cleaned1["I_v_p1"]
            )

            # System 2: complex signals
            df_cleaned2["i_02"] = (
                df_cleaned2["R_i_02"] + 1j * df_cleaned3["I_i_02"]
            )
            df_cleaned2["i_n2"] = (
                df_cleaned3["R_i_n2"] + 1j * df_cleaned3["I_i_n2"]
            )
            df_cleaned2["i_p2"] = (
                df_cleaned2["R_i_p2"] + 1j * df_cleaned3["I_i_p2"]
            )

            df_cleaned3["v_02"] = (
                df_cleaned2["R_v_02"] + 1j * df_cleaned2["I_v_02"]
            )
            df_cleaned3["v_n2"] = (
                df_cleaned2["R_v_n2"] + 1j * df_cleaned2["I_v_n2"]
            )
            df_cleaned3["v_p2"] = (
                df_cleaned2["R_v_p2"] + 1j * df_cleaned2["I_v_p2"]
            )

            if ss_snap == 1:
                # Extract windows
                time_w1 = df_cleaned1.loc[Tdelay : Tobs0 + Tdelay]
                time_w2 = df_cleaned2.loc[Tdelay : Tobs0 + Tdelay]
                time_w3 = df_cleaned3.loc[Tdelay : Tobs0 + Tdelay]
            
            if ss_snap == 0:
                # Extract windows
                time_w1 = df_cleaned1.loc[Tinit:Tinit + Tobs]
                time_w2 = df_cleaned2.loc[Tinit:Tinit + Tobs]
                time_w3 = df_cleaned3.loc[Tinit:Tinit + Tobs]

            # Subtract steady-state
            v_p1 = (
                time_w1["v_p1"].to_numpy()
                - steady_state_signals["v_p1_ss"]
            )
            v_n1 = (
                time_w1["v_n1"].to_numpy()
                - steady_state_signals["v_n1_ss"]
            )
            v_01 = (
                time_w1["v_01"].to_numpy()
                - steady_state_signals["v_01_ss"]
            )
            i_p1 = (
                time_w1["i_p1"].to_numpy()
                - steady_state_signals["i_p1_ss"]
            )
            i_n1 = (
                time_w1["i_n1"].to_numpy()
                - steady_state_signals["i_n1_ss"]
            )
            i_01 = (
                time_w1["i_01"].to_numpy()
                - steady_state_signals["i_01_ss"]
            )

            v_p2 = (
                time_w3["v_p2"].to_numpy()
                - steady_state_signals["v_p2_ss"]
            )
            v_n2 = (
                time_w3["v_n2"].to_numpy()
                - steady_state_signals["v_n2_ss"]
            )
            v_02 = (
                time_w3["v_02"].to_numpy()
                - steady_state_signals["v_02_ss"]
            )
            i_p2 = (
                time_w2["i_p2"].to_numpy()
                - steady_state_signals["i_p2_ss"]
            )
            i_n2 = (
                time_w2["i_n2"].to_numpy()
                - steady_state_signals["i_n2_ss"]
            )
            i_02 = (
                time_w2["i_02"].to_numpy()
                - steady_state_signals["i_02_ss"]
            )

            # FFT
            Vp1 = np.fft.fft(v_p1) / len(v_p1)
            Vn1 = np.fft.fft(v_n1) / len(v_n1)
            Ip1 = np.fft.fft(i_p1) / len(i_p1)
            In1 = np.fft.fft(i_n1) / len(i_n1)

            Vp2 = np.fft.fft(v_p2) / len(v_p2)
            Vn2 = np.fft.fft(v_n2) / len(v_n2)
            Ip2 = np.fft.fft(i_p2) / len(i_p2)
            In2 = np.fft.fft(i_n2) / len(i_n2)

            frequencies = np.fft.fftfreq(len(v_p1), d=Tobs0)
            positive_freqs = frequencies[frequencies >= 0]

            Vp1 = Vp1[: len(positive_freqs)]
            Vn1 = Vn1[: len(positive_freqs)]
            Ip1 = Ip1[: len(positive_freqs)]
            In1 = In1[: len(positive_freqs)]

            Vp2 = Vp2[: len(positive_freqs)]
            Vn2 = Vn2[: len(positive_freqs)]
            Ip2 = Ip2[: len(positive_freqs)]
            In2 = In2[: len(positive_freqs)]

            wd = round(fd0[i] / fs)

            Ypp1 = Ip1[wd] / Vp1[wd]
            Ynp1 = In1[wd] / Vp1[wd]

            Ypp2 = Ip2[wd] / Vp2[wd]
            Ynp2 = In2[wd] / Vp2[wd]

            # -------------------------------------------------------------
            # Injection in n-axis
            # -------------------------------------------------------------
            pn0_block.parameters(
                Vdist_p=0.0, Vdist_n=Vdist, Vdist_0=0.0
            )

            run_perturbation(
                simset,
                task,
                ss_snap,
                snapshot_name,
                project_name,
                fortran_ext,
                out_filename,
            )

            # Load .out files again
            df1 = pd.read_csv(
                steady_out_file1[0],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned1 = df1.dropna(axis=1, how="all")

            df2 = pd.read_csv(
                steady_out_file1[1],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned2 = df2.dropna(axis=1, how="all")

            df3 = pd.read_csv(
                steady_out_file1[2],
                delimiter=r"\s+",
                header=None,
                engine="python",
            )
            df_cleaned3 = df3.dropna(axis=1, how="all")

            # Assign column names
            df_cleaned1.columns = [
                "time",
                "R_v_01",
                "R_v_p1",
                "R_v_n1",
                "I_v_01",
                "I_v_p1",
                "I_v_n1",
                "R_i_01",
                "R_i_p1",
                "R_i_n1",
                "I_i_01",
            ]

            df_cleaned2.columns = [
                "time",
                "I_i_p1",
                "I_i_n1",
                "R_v_02",
                "R_v_p2",
                "R_v_n2",
                "I_v_02",
                "I_v_p2",
                "I_v_n2",
                "R_i_02",
                "R_i_p2",
            ]

            df_cleaned3.columns = [
                "time",
                "R_i_n2",
                "I_i_02",
                "I_i_p2",
                "I_i_n2",
            ]

            # Set index
            df_cleaned1 = df_cleaned1.set_index("time")
            df_cleaned2 = df_cleaned2.set_index("time")
            df_cleaned3 = df_cleaned3.set_index("time")

            # System 1: complex signals
            df_cleaned1["i_01"] = (
                df_cleaned1["R_i_01"] + 1j * df_cleaned1["I_i_01"]
            )
            df_cleaned1["i_n1"] = (
                df_cleaned1["R_i_n1"] + 1j * df_cleaned2["I_i_n1"]
            )
            df_cleaned1["i_p1"] = (
                df_cleaned1["R_i_p1"] + 1j * df_cleaned2["I_i_p1"]
            )

            df_cleaned1["v_01"] = (
                df_cleaned1["R_v_01"] + 1j * df_cleaned1["I_v_01"]
            )
            df_cleaned1["v_n1"] = (
                df_cleaned1["R_v_n1"] + 1j * df_cleaned1["I_v_n1"]
            )
            df_cleaned1["v_p1"] = (
                df_cleaned1["R_v_p1"] + 1j * df_cleaned1["I_v_p1"]
            )

            # System 2: complex signals
            df_cleaned2["i_02"] = (
                df_cleaned2["R_i_02"] + 1j * df_cleaned3["I_i_02"]
            )
            df_cleaned2["i_n2"] = (
                df_cleaned3["R_i_n2"] + 1j * df_cleaned3["I_i_n2"]
            )
            df_cleaned2["i_p2"] = (
                df_cleaned2["R_i_p2"] + 1j * df_cleaned3["I_i_p2"]
            )

            df_cleaned3["v_02"] = (
                df_cleaned2["R_v_02"] + 1j * df_cleaned2["I_v_02"]
            )
            df_cleaned3["v_n2"] = (
                df_cleaned2["R_v_n2"] + 1j * df_cleaned2["I_v_n2"]
            )
            df_cleaned3["v_p2"] = (
                df_cleaned2["R_v_p2"] + 1j * df_cleaned2["I_v_p2"]
            )

            if ss_snap == 1:
                # Extract windows
                time_w1 = df_cleaned1.loc[Tdelay : Tobs0 + Tdelay]
                time_w2 = df_cleaned2.loc[Tdelay : Tobs0 + Tdelay]
                time_w3 = df_cleaned3.loc[Tdelay : Tobs0 + Tdelay]
            
            if ss_snap == 0:
                # Extract windows
                time_w1 = df_cleaned1.loc[Tinit:Tinit + Tobs]
                time_w2 = df_cleaned2.loc[Tinit:Tinit + Tobs]
                time_w3 = df_cleaned3.loc[Tinit:Tinit + Tobs]

            # Subtract steady-state
            v_p1 = (
                time_w1["v_p1"].to_numpy()
                - steady_state_signals["v_p1_ss"]
            )
            v_n1 = (
                time_w1["v_n1"].to_numpy()
                - steady_state_signals["v_n1_ss"]
            )
            v_01 = (
                time_w1["v_01"].to_numpy()
                - steady_state_signals["v_01_ss"]
            )
            i_p1 = (
                time_w1["i_p1"].to_numpy()
                - steady_state_signals["i_p1_ss"]
            )
            i_n1 = (
                time_w1["i_n1"].to_numpy()
                - steady_state_signals["i_n1_ss"]
            )
            i_01 = (
                time_w1["i_01"].to_numpy()
                - steady_state_signals["i_01_ss"]
            )

            v_p2 = (
                time_w3["v_p2"].to_numpy()
                - steady_state_signals["v_p2_ss"]
            )
            v_n2 = (
                time_w3["v_n2"].to_numpy()
                - steady_state_signals["v_n2_ss"]
            )
            v_02 = (
                time_w3["v_02"].to_numpy()
                - steady_state_signals["v_02_ss"]
            )
            i_p2 = (
                time_w2["i_p2"].to_numpy()
                - steady_state_signals["i_p2_ss"]
            )
            i_n2 = (
                time_w2["i_n2"].to_numpy()
                - steady_state_signals["i_n2_ss"]
            )
            i_02 = (
                time_w2["i_02"].to_numpy()
                - steady_state_signals["i_02_ss"]
            )

            Vp1 = np.fft.fft(v_p1) / len(v_p1)
            Vn1 = np.fft.fft(v_n1) / len(v_n1)
            Ip1 = np.fft.fft(i_p1) / len(i_p1)
            In1 = np.fft.fft(i_n1) / len(i_n1)

            Vp2 = np.fft.fft(v_p2) / len(v_p2)
            Vn2 = np.fft.fft(v_n2) / len(v_n2)
            Ip2 = np.fft.fft(i_p2) / len(i_p2)
            In2 = np.fft.fft(i_n2) / len(i_n2)

            Vp1 = Vp1[: len(positive_freqs)]
            Vn1 = Vn1[: len(positive_freqs)]
            Ip1 = Ip1[: len(positive_freqs)]
            In1 = In1[: len(positive_freqs)]

            Vp2 = Vp2[: len(positive_freqs)]
            Vn2 = Vn2[: len(positive_freqs)]
            Ip2 = Ip2[: len(positive_freqs)]
            In2 = In2[: len(positive_freqs)]

            Ypn1 = Ip1[wd] / Vn1[wd]
            Ynn1 = In1[wd] / Vn1[wd]

            Ypn2 = Ip2[wd] / Vn2[wd]
            Ynn2 = In2[wd] / Vn2[wd]

            # Assemble 2x2 pn-admittance matrices
            Ysys10 = np.array([[Ypp1, Ypn1], [Ynp1, Ynn1]])
            Ysys20 = np.array([[Ypp2, Ypn2], [Ynp2, Ynn2]])

            Ysys1[fd0[i]] = Ysys10
            Ysys2[fd0[i]] = Ysys20

            Y_mag1[fd0[i]], Y_rad1[fd0[i]] = cartz2pol(Ysys10)
            Y_mag2[fd0[i]], Y_rad2[fd0[i]] = cartz2pol(Ysys20)

        # End of 0pn scanning
        pscad.quit()

        # Save results in results_folder using existing utilities
        Path(results_folder).mkdir(exist_ok=True)

        save_Ysyspn_txt("Ysys1", Ysys1, fd0, results_folder)
        save_Ysyspn_txt("Ysys2", Ysys2, fd0, results_folder)

        print("---> SIaD-Tool has finished.")
        print(
            "---> Results stored in the files: Ysys1.txt and Ysys2.txt."
        )

        fig1 = pn0_plot(fd0, Y_mag1, Y_rad1)
        fig1.savefig(f"{results_folder}/Ysys1.pdf", format="pdf", bbox_inches="tight")
        plt.close(fig1)

        fig2 = pn0_plot(fd0, Y_mag2, Y_rad2)
        fig2.savefig(f"{results_folder}/Ysys2.pdf", format="pdf", bbox_inches="tight")
        plt.close(fig2)

        return fd0, Ysys1, Ysys2

    # If no valid scanner_selector is provided
    raise ValueError("Invalid scanner_selector. Must be 1 (ABC), 2 (dq0) or 3 (0pn).")

