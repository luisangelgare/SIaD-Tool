
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


# =============================================================================
#  IMPORTS
# =============================================================================

import time
from FDScanning import FDScanningTool
from GNC_Module import Run_GNC
from MA_Module import Run_Modal_Analysis
from PM_Module import Run_Phase_Margins
from PA_Module import Run_Passivity_Analysis


# =============================================================================
#  USER-DEFINED PARAMETERS
# =============================================================================

workspace_name = 'workspace'                  # PSCAD workspace name
project_name = "DCS1_Released_20140630"       # PSCAD project name
fortran_ext = ".gf46"           # Fortran or intel version: ".gf46", ".gf81", etc.
out_filename = "meas_dq"        # Proposed name of the file .out (must be assigned here)
results_folder = "Results_dq"   # Proposed name of the results folder

Tinit = 4                     # Initialization time (s)
fs = 1                         # Resolution frequency in Hz
f_max = 500                    # Upper limit of the frequency band
f_points = 180                 # Number of frequency points
f_scale = 1                    # Frequency vector scale type: 1 for log and 2 for linear
step_time = 50                 # Micro-seconds fixed step time for PSCAD 

ss_snap = 1                    # Take a snapshot in 1st simulation: 1 yes and 0 not
Sbase = 100                    # Base power (MVA) 
Vbase = 230E3                  # Base voltage (V)
f0 = 50                        # Base frequency (Hz)
Vperturbation = 0.03           # Percentage of the nominal voltage for the perturbation

# ---> Available scanning strategies:
# 1 -> Voltage perturbation
# 2 -> Current perturbation (NOT AVAILABLE YET)
scanner_type = 1
# 1 -> Single-tone perturbation
# 2 -> PRBS perturbation (NOT AVAILABLE YET)
# 3 -> Multi-tone perturbation (NOT AVAILABLE YET)
signal_type = 1
# 1 -> ABC scan
# 2 -> dq0 scan
# 3 -> 0pn scan
scanner_selector = 2


# =============================================================================
#  SIaD-Tool Scanning framework
# =============================================================================

start_time = time.time()

[fd0, Ysys1, Ysys2] = FDScanningTool(
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
                        results_folder
                    )

# =============================================================================
#  SIaD-Tool REPORT
# =============================================================================


gnc_result = Run_GNC(fd0, Ysys1, Ysys2, 
                     outcomes_dir=results_folder, 
                     scanner_selector=scanner_selector)

ma_results = Run_Modal_Analysis(fd0,
                                Ysys1,
                                Ysys2,
                                subsys_labels=("Sys1", "Sys2"),
                                outcomes_dir=results_folder)

PM_table = Run_Phase_Margins(fd0, Ysys1, Ysys2, 
                             scanner_selector,
                             outcomes_dir=results_folder)

PA_results = Run_Passivity_Analysis(fd0,
                                    Ysys1,
                                    Ysys2,
                                    outcomes_dir=results_folder)


end_time = time.time()
execution_time = end_time - start_time
print("---")
print(f"---> Total execution time of the SIaD-Tool: {execution_time} seconds")

