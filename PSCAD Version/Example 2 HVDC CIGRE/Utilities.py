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
#  IMPORTS
# =============================================================================

import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
#  AUXILIARY FUNCTIONS
# =============================================================================

def find_out_file(project_name, fortran_ext, out_filename):
    folder = Path(f"{project_name}{fortran_ext}")
    matches = [f for f in folder.iterdir()
               if f.name.startswith(out_filename + "_") and f.name.endswith(".out")]
    if not matches:
        raise FileNotFoundError(f"No .out file found for prefix {out_filename}_ in {folder}")
    matches.sort()
    return matches[-1]  # archivo más reciente

def find_all_out_files(project_name, fortran_ext, out_filename):
    folder = Path(f"{project_name}{fortran_ext}")

    # Buscar archivos con patrón: out_filename_XX.out
    matches = []
    for f in folder.iterdir():
        name = f.name
        if name.startswith(out_filename + "_") and name.endswith(".out"):
            # Extraer el número final XX
            try:
                suffix = name.split("_")[-1].replace(".out", "")
                num = int(suffix)
                matches.append((num, f))
            except ValueError:
                continue  # ignorar archivos que no sigan el patrón

    if not matches:
        raise FileNotFoundError(
            f"No .out files found for prefix {out_filename}_ in {folder}"
        )

    # Ordenar por número (1,2,3,...)
    matches.sort(key=lambda x: x[0])

    # Devolver solo las rutas, en orden
    return [f for (_, f) in matches]


def run_steady_state(simset,
                    task,
                    snapshot_name,
                    Tinit,
                    Tobs0,
                    project_name,
                    fortran_ext,
                    out_filename,
                    scanner_selector,
                    ss_snap,
                    steady_state_store):
    """
    -------------------------------------------------------------------------
    Stage 1: Steady-State Simulation + Snapshot Capture
    -------------------------------------------------------------------------
    Executes a single PSCAD simulation using the configuration already
    defined in the project. A timed snapshot is taken at t = Tinit.
    The .out file is then processed to extract steady-state waveforms.
    -------------------------------------------------------------------------
    """
    if ss_snap == 1:
        
        # Configure snapshot parameters (minimal overrides)
        task.overrides(
            start_method=0,                   # Standard startup
            timed_snapshots=1,                # Enable single snapshot
            snapshot_file=snapshot_name + ".snp",
            snap_time=Tinit,                   # Snapshot time [s]
            save_channels_file=out_filename + ".out",
            save_channels=1
        )

    if ss_snap == 0:
        
        # Configure parameters (without taking a snapshot)
        task.overrides(
            start_method=0,                   # Standard startup
            timed_snapshots=0,                # Enable single snapshot
            save_channels_file=out_filename + ".out",
            save_channels=1
        )

    # Run simulation
    simset.run()
    
    # abc scanner: assign column names and extract window
    if scanner_selector == 1:
        out_file_path1 = find_all_out_files(project_name,
                                            fortran_ext,
                                            out_filename)
    
        # Load .out files
        df1 = pd.read_csv(out_file_path1[0], delimiter=r'\s+', header=None, engine='python')
        df_cleaned1 = df1.dropna(axis=1, how='all')
    
        df2 = pd.read_csv(out_file_path1[1], delimiter=r'\s+', header=None, engine='python')
        df_cleaned2 = df2.dropna(axis=1, how='all')
    
        # Assign column names
        df_cleaned1.columns = [
            'time', 'theta',
            'v_a1','v_b1','v_c1',       
            'i_a1','i_b1','i_c1',
            'v_q2','v_d2',
            'v_a2'
        ]
    
        df_cleaned2.columns = [
            'time', 
            'v_b2','v_c2',       
            'i_a2','i_b2','i_c2'
        ]
    
        # Set index and align by time
        df_cleaned1 = df_cleaned1.set_index('time')
        df_cleaned2 = df_cleaned2.set_index('time')
        
        # Extract windows
        time_w1 = df_cleaned1.loc[Tinit:Tinit + Tobs0]
        time_w2 = df_cleaned2.loc[Tinit:Tinit + Tobs0]
        # Steady-state exact sample (theta only exists in df_cleaned1)
        time_ss = df_cleaned1.loc[Tinit]
        
        # Store steady-state signals
        steady_state_store['i_a1_ss'] = time_w1['i_a1'].to_numpy()
        steady_state_store['i_b1_ss'] = time_w1['i_b1'].to_numpy()
        steady_state_store['i_c1_ss'] = time_w1['i_c1'].to_numpy()
    
        steady_state_store['v_a1_ss'] = time_w1['v_a1'].to_numpy()
        steady_state_store['v_b1_ss'] = time_w1['v_b1'].to_numpy()
        steady_state_store['v_c1_ss'] = time_w1['v_c1'].to_numpy()
    
        steady_state_store['i_a2_ss'] = time_w2['i_a2'].to_numpy()
        steady_state_store['i_b2_ss'] = time_w2['i_b2'].to_numpy()
        steady_state_store['i_c2_ss'] = time_w2['i_c2'].to_numpy()
    
        steady_state_store['v_a2_ss'] = time_w1['v_a2'].to_numpy()
        steady_state_store['v_b2_ss'] = time_w2['v_b2'].to_numpy()
        steady_state_store['v_c2_ss'] = time_w2['v_c2'].to_numpy()
    
        steady_state_store['vq_ss'] = float(time_ss['v_q2'])
        steady_state_store['vd_ss'] = float(time_ss['v_d2'])
        steady_state_store['theta_ss'] = float(time_ss['theta'])

    # dq0 scanner: assign column names and extract window
    if scanner_selector == 2:
        out_file_path1 = find_out_file(project_name,
                                       fortran_ext,
                                       out_filename)

        # Load .out file and extract steady-state window
        df = pd.read_csv(out_file_path1, delimiter=r'\s+', header=None, engine='python')
        df_cleaned = df.dropna(axis=1, how='all')
        df_cleaned.columns = [
            'time', 'theta', 'i_q1', 'i_d1', 'v_q1', 'v_d1',
            'i_q2', 'i_d2', 'v_q2', 'v_d2'
        ]
        df_cleaned = df_cleaned.set_index('time')
        # df_cleaned2 = df_cleaned2.set_index('time')
        time_w = df_cleaned.loc[Tinit:Tinit + Tobs0]
        time_ss = df_cleaned.loc[Tinit]

        # Store steady-state signals
        steady_state_store['v_q1_ss'] = time_w['v_q1'].to_numpy()
        steady_state_store['v_d1_ss'] = time_w['v_d1'].to_numpy()
        steady_state_store['i_q1_ss'] = time_w['i_q1'].to_numpy()
        steady_state_store['i_d1_ss'] = time_w['i_d1'].to_numpy()
        steady_state_store['v_q2_ss'] = time_w['v_q2'].to_numpy()
        steady_state_store['v_d2_ss'] = time_w['v_d2'].to_numpy()
        steady_state_store['i_q2_ss'] = time_w['i_q2'].to_numpy()
        steady_state_store['i_d2_ss'] = time_w['i_d2'].to_numpy()
        steady_state_store['vq_ss'] = float(time_ss['v_q2'])
        steady_state_store['vd_ss'] = float(time_ss['v_d2'])
        steady_state_store['iq_ss'] = float(time_ss['i_q2'])
        steady_state_store['id_ss'] = float(time_ss['i_d2'])
        steady_state_store['theta_ss'] = float(time_ss['theta'])
        

    # 0pn scanner: assign column names and extract window
    if scanner_selector == 3:
        out_file_path1 = find_all_out_files(project_name,
                                            fortran_ext,
                                            out_filename)
    
        # Load .out files
        df1 = pd.read_csv(out_file_path1[0], delimiter=r'\s+', header=None, engine='python')
        df_cleaned1 = df1.dropna(axis=1, how='all')
    
        df2 = pd.read_csv(out_file_path1[1], delimiter=r'\s+', header=None, engine='python')
        df_cleaned2 = df2.dropna(axis=1, how='all')
    
        df3 = pd.read_csv(out_file_path1[2], delimiter=r'\s+', header=None, engine='python')
        df_cleaned3 = df3.dropna(axis=1, how='all')
    
        # Assign column names
        df_cleaned1.columns = [
            'time', 'theta',
            'R_v_01','R_v_p1','R_v_n1',
            'I_v_01','I_v_p1','I_v_n1',       
            'R_i_01','R_i_p1','R_i_n1'
        ]
    
        df_cleaned2.columns = [
            'time', 
            'I_i_01','I_i_p1','I_i_n1',
            'v_q2','v_d2',
            'R_v_02','R_v_p2','R_v_n2',
            'I_v_02','I_v_p2' 
        ]
    
        df_cleaned3.columns = [
            'time',
            'I_v_n2',
            'R_i_02','R_i_p2','R_i_n2',
            'I_i_02','I_i_p2','I_i_n2' 
        ]
    
        # Set index and align by time
        df_cleaned1 = df_cleaned1.set_index('time')
        df_cleaned2 = df_cleaned2.set_index('time')
        df_cleaned3 = df_cleaned3.set_index('time')
   
        # --- System 1: complex signals ---
        df_cleaned1['i_01'] = df_cleaned1['R_i_01'] + 1j * df_cleaned2['I_i_01']
        df_cleaned1['i_n1'] = df_cleaned1['R_i_n1'] + 1j * df_cleaned2['I_i_n1']
        df_cleaned1['i_p1'] = df_cleaned1['R_i_p1'] + 1j * df_cleaned2['I_i_p1']
    
        df_cleaned1['v_01'] = df_cleaned1['R_v_01'] + 1j * df_cleaned1['I_v_01']
        df_cleaned1['v_n1'] = df_cleaned1['R_v_n1'] + 1j * df_cleaned1['I_v_n1']
        df_cleaned1['v_p1'] = df_cleaned1['R_v_p1'] + 1j * df_cleaned1['I_v_p1']
    
        # --- System 2: complex signals ---
        df_cleaned2['i_02'] = df_cleaned3['R_i_02'] + 1j * df_cleaned3['I_i_02']
        df_cleaned2['i_n2'] = df_cleaned3['R_i_n2'] + 1j * df_cleaned3['I_i_n2']
        df_cleaned2['i_p2'] = df_cleaned3['R_i_p2'] + 1j * df_cleaned3['I_i_p2']
    
        df_cleaned3['v_02'] = df_cleaned2['R_v_02'] + 1j * df_cleaned2['I_v_02']
        df_cleaned3['v_n2'] = df_cleaned2['R_v_n2'] + 1j * df_cleaned3['I_v_n2']
        df_cleaned3['v_p2'] = df_cleaned2['R_v_p2'] + 1j * df_cleaned2['I_v_p2']
    
        # Extract windows
        time_w1 = df_cleaned1.loc[Tinit:Tinit + Tobs0]
        time_w2 = df_cleaned2.loc[Tinit:Tinit + Tobs0]
        time_w3 = df_cleaned3.loc[Tinit:Tinit + Tobs0]
    
        # Steady-state exact sample (theta only exists in df_cleaned1)
        time_ss = df_cleaned1.loc[Tinit]
    
        # Store steady-state signals
        steady_state_store['i_01_ss'] = time_w1['i_01'].to_numpy()
        steady_state_store['i_n1_ss'] = time_w1['i_n1'].to_numpy()
        steady_state_store['i_p1_ss'] = time_w1['i_p1'].to_numpy()
    
        steady_state_store['v_01_ss'] = time_w1['v_01'].to_numpy()
        steady_state_store['v_n1_ss'] = time_w1['v_n1'].to_numpy()
        steady_state_store['v_p1_ss'] = time_w1['v_p1'].to_numpy()
    
        steady_state_store['i_02_ss'] = time_w2['i_02'].to_numpy()
        steady_state_store['i_n2_ss'] = time_w2['i_n2'].to_numpy()
        steady_state_store['i_p2_ss'] = time_w2['i_p2'].to_numpy()
    
        steady_state_store['v_02_ss'] = time_w3['v_02'].to_numpy()
        steady_state_store['v_n2_ss'] = time_w3['v_n2'].to_numpy()
        steady_state_store['v_p2_ss'] = time_w3['v_p2'].to_numpy()
    
        steady_state_store['vq_ss'] = float(df_cleaned2.loc[Tinit]['v_q2'].real)
        steady_state_store['vd_ss'] = float(df_cleaned2.loc[Tinit]['v_d2'].real)
        steady_state_store['theta_ss'] = float(time_ss['theta'].real)


    return snapshot_name + ".snp"



def run_perturbation(simset,
                      task,
                      ss_snap,
                      snapshot_name,
                      project_name,
                      fortran_ext,
                      out_filename):
    """
    -------------------------------------------------------------------------
    Stage 2: Simulation from Snapshot (Perturbation Runs)
    -------------------------------------------------------------------------
    Restarts PSCAD simulation from the previously generated snapshot.
    Only minimal overrides are applied to ensure PSCAD stability.
    -------------------------------------------------------------------------
    """
    if ss_snap == 1:
        task.overrides(
            start_method=1,                       # Restart from snapshot
            timed_snapshots=0,                    # Disable timed snapshots
            startup_inputfile=snapshot_name + ".snp",
            save_channels_file=out_filename + ".out",
            save_channels=1
        )
        
    if ss_snap == 0:
        task.overrides(
        start_method=0,                 # No snapshot, start from t=0
        timed_snapshots=0,              # No timed snapshots
        save_channels_file=out_filename + ".out",
        save_channels=1
        )

    simset.run()
    
    out_file_path = find_out_file(project_name,
                                   fortran_ext,
                                   out_filename)
    
    return out_file_path



def project_file(project_name: str):
    """
    Utility function to locate the PSCAD workspace or project file.
    Searches for .pswx first, then any matching extension.
    """

    try:
        base_dir = Path(__file__).resolve().parent
    except NameError:
        base_dir = Path.cwd()

    preferred = list(base_dir.glob(f"{project_name}.pswx"))
    if preferred:
        return preferred[0]

    candidates = list(base_dir.glob(f"{project_name}.*"))
    if not candidates:
        raise FileNotFoundError(f"File '{project_name}' does not exist in {base_dir}")

    return candidates[0]


def save_Ysysqd_txt(filename, Ysys, freqs, results_folder):
    with open(f"{results_folder}/{filename}.txt", "w") as f:
        f.write("f_Hz   Yqq   Yqd   Ydq   Ydd\n")
        for freq in freqs:
            mat = Ysys[freq]
            f.write(f"{freq}   {mat[0,0]}   {mat[0,1]}   {mat[1,0]}   {mat[1,1]}\n")
            

def save_Ysyspn_txt(filename, Ysys, freqs, results_folder):
    with open(f"{results_folder}/{filename}.txt", "w") as f:
        f.write("f_Hz   Ypp   Ypn   Ynp   Ynn\n")
        for freq in freqs:
            mat = Ysys[freq]
            f.write(f"{freq}   {mat[0,0]}   {mat[0,1]}   {mat[1,0]}   {mat[1,1]}\n")
            

def save_YsysABC_txt(filename, Ysys, freqs, results_folder):
    with open(f"{results_folder}/{filename}.txt", "w") as f:
        f.write("f_Hz   Yaa   Yab   Yac   Yba   Ybb   Ybc  Yca   Ycb   Ycc\n")
        for freq in freqs:
            mat = Ysys[freq]
            f.write(f"{freq}   {mat[0,0]}   {mat[0,1]}   {mat[0,2]}   {mat[1,0]}   {mat[1,1]}   {mat[1,2]}   {mat[2,0]}   {mat[2,1]}   {mat[2,2]}\n")
            


# =============================================================================
#  MATHEMATICAL UTILITIES
# =============================================================================

def cartz2pol(z):
    """Convert complex number to magnitude and angle."""
    return np.abs(z), np.angle(z)


def parse_complex(complex_str):
    """Convert PSCAD complex string to Python complex."""
    return complex(complex_str.replace('i', 'j'))


def choose(index, *options):
    """Utility for formatted printing."""
    return options[index]


def disable_except(blocks, keep_names):
    """Disable PSCAD blocks except those listed in keep_names."""
    for blk in blocks:
        if blk.parameters()['Name'] not in keep_names:
            blk.disable()
            
            
            
# =============================================================================
#  GRAPH UTILITIES
# =============================================================================
       

def ABC_plot_response(fsampling, Fs_mag, Fs_ang, fd0, Y_mag, Y_rad):
    """
    Función para graficar la magnitud en dB y fase en grados vs frecuencia, 
    junto con un escaneo de frecuencias para el sistema ABC.

    Parámetros:
    - fsampling: Vector de frecuencias de muestreo (Hz)
    - Fs_mag: Magnitud de la función de transferencia
    - Fs_ang: Ángulo (fase) de la función de transferencia en radianes
    - fd0: Frecuencia límite inferior y superior para graficar
    - Y_mag: Diccionario con las magnitudes Y_ABC de las matrices
    - Y_rad: Diccionario con las fases Y_ABC en radianes de las matrices
    """
    
    # Convertir magnitud a decibelios y ángulo a grados
    Fs_mag_dB = 20 * np.log10(Fs_mag)
    Fs_ang_deg = (180 / np.pi) * Fs_ang
    
    # Definir límites del eje X (frecuencia)
    x_low_axis = fd0[0]
    x_up_axis = fd0[-1]
    
    # Ajustar los límites del eje Y para magnitud y fase
    mag_min = np.min(Fs_mag_dB[1:])
    mag_max = np.max(Fs_mag_dB)
    phase_min = np.min(Fs_ang_deg)
    phase_max = np.max(Fs_ang_deg)
    
    # Configuración de estilos globales para las gráficas (incluyendo tamaño de textos)
    plt.rcParams.update({
        'axes.titlesize': 20,       # Tamaño de los títulos de los ejes
        'axes.labelsize': 20,       # Tamaño de las etiquetas de los ejes (X, Y)
        'xtick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje X
        'ytick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje Y
        'legend.fontsize': 20,      # Tamaño del texto de las leyendas
        'figure.titlesize': 20,     # Tamaño del título de la figura (si lo hubiera)
        'lines.linewidth': 2        # Ancho de las líneas de las curvas
    })
    
    # Crear subplots
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    # Top row - Magnitude in dB vs log(fsampling)
    axs[0, 0].semilogx(fsampling, Fs_mag_dB, label='Theoretical response', color='k')
    axs[0, 0].set_title(r'$Y_{AA}(s)$')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].set_xlim([x_low_axis, x_up_axis])
    axs[0, 0].set_ylim([mag_min - 3, mag_max + 3])
    axs[0, 0].grid(True, which='major', axis='both')
    axs[0, 0].grid(True, which='minor', axis='x')

    axs[0, 1].semilogx(fsampling, Fs_mag_dB, label='Theoretical response', color='k')
    axs[0, 1].set_title(r'$Y_{BB}(s)$')
    axs[0, 1].set_xlim([x_low_axis, x_up_axis])
    axs[0, 1].set_ylim([mag_min - 3, mag_max + 3])
    axs[0, 1].grid(True, which='major', axis='both')
    axs[0, 1].grid(True, which='minor', axis='x')

    axs[0, 2].semilogx(fsampling, Fs_mag_dB, label='Theoretical response', color='k')
    axs[0, 2].set_title(r'$Y_{CC}(s)$')
    axs[0, 2].set_xlim([x_low_axis, x_up_axis])
    axs[0, 2].set_ylim([mag_min - 3, mag_max + 3])
    axs[0, 2].grid(True, which='major', axis='both')
    axs[0, 2].grid(True, which='minor', axis='x')
    
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[1, 0].semilogx(fsampling, Fs_ang_deg, label='Phase (degrees)', color='k')
    axs[1, 0].set_xlabel('Frequency (Hz)')
    axs[1, 0].set_ylabel('Phase (Degrees)')
    axs[1, 0].set_xlim([x_low_axis, x_up_axis])
    axs[1, 0].set_ylim([phase_min - 3, phase_max + 3])
    axs[1, 0].grid(True, which='major', axis='both')
    axs[1, 0].grid(True, which='minor', axis='x')

    axs[1, 1].semilogx(fsampling, Fs_ang_deg, label='Phase (degrees)', color='k')
    axs[1, 1].set_xlabel('Frequency (Hz)')
    axs[1, 1].set_xlim([x_low_axis, x_up_axis])
    axs[1, 1].set_ylim([phase_min - 3, phase_max + 3])
    axs[1, 1].grid(True, which='major', axis='both')
    axs[1, 1].grid(True, which='minor', axis='x')

    axs[1, 2].semilogx(fsampling, Fs_ang_deg, label='Phase (degrees)', color='k')
    axs[1, 2].set_xlabel('Frequency (Hz)')
    axs[1, 2].set_xlim([x_low_axis, x_up_axis])
    axs[1, 2].set_ylim([phase_min - 3, phase_max + 3])
    axs[1, 2].grid(True, which='major', axis='both')
    axs[1, 2].grid(True, which='minor', axis='x')

    # Agregar puntos rojos para "ABC frequency scan"
    frequencies = list(Y_mag.keys())
    for freq in frequencies:
        Y_ABC_mag = Y_mag[freq]
        Y_ABC_angle = Y_rad[freq]

        axs[0, 0].semilogx(freq, 20 * np.log10(Y_ABC_mag[0, 0]), 'ro', label='ABC frequency scan')
        axs[0, 1].semilogx(freq, 20 * np.log10(Y_ABC_mag[1, 1]), 'ro', label='ABC frequency scan')
        axs[0, 2].semilogx(freq, 20 * np.log10(Y_ABC_mag[2, 2]), 'ro', label=f'{freq} Hz')

        axs[1, 0].semilogx(freq, (180 / np.pi) * Y_ABC_angle[0, 0], 'ro', label='ABC frequency scan')
        axs[1, 1].semilogx(freq, (180 / np.pi) * Y_ABC_angle[1, 1], 'ro', label='ABC frequency scan')
        axs[1, 2].semilogx(freq, (180 / np.pi) * Y_ABC_angle[2, 2], 'ro', label='ABC frequency scan')

    axs[0, 2].legend(['Theoretical response', 'ABC frequency scan'], loc='lower left')

    # Ajuste final
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()



def ABC_plot(fd0, Y_mag, Y_rad):

    x_low_axis = fd0[0]
    x_up_axis = fd0[-1]

    plt.rcParams.update({
        'axes.titlesize': 20,
        'axes.labelsize': 20,
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        'legend.fontsize': 20,
        'figure.titlesize': 20,
        'lines.linewidth': 2
    })

    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    # === ORDENAR FRECUENCIAS ===
    freqs = sorted(Y_mag.keys())

    # === CONSTRUIR VECTORES ===
    YAA_mag = [20*np.log10(Y_mag[f][0,0]) for f in freqs]
    YBB_mag = [20*np.log10(Y_mag[f][1,1]) for f in freqs]
    YCC_mag = [20*np.log10(Y_mag[f][2,2]) for f in freqs]

    YAA_ang = [(180/np.pi)*Y_rad[f][0,0] for f in freqs]
    YBB_ang = [(180/np.pi)*Y_rad[f][1,1] for f in freqs]
    YCC_ang = [(180/np.pi)*Y_rad[f][2,2] for f in freqs]

    # === MAGNITUD ===
    axs[0,0].semilogx(freqs, YAA_mag, 'r-')
    axs[0,0].set_title(r'$Y_{AA}(s)$')
    axs[0,0].set_ylabel('Magnitude (dB)')
    axs[0,0].set_xlim([x_low_axis, x_up_axis])
    axs[0,0].grid(True, which='both')

    axs[0,1].semilogx(freqs, YBB_mag, 'r-')
    axs[0,1].set_title(r'$Y_{BB}(s)$')
    axs[0,1].set_xlim([x_low_axis, x_up_axis])
    axs[0,1].grid(True, which='both')

    axs[0,2].semilogx(freqs, YCC_mag, 'r-')
    axs[0,2].set_title(r'$Y_{CC}(s)$')
    axs[0,2].set_xlim([x_low_axis, x_up_axis])
    axs[0,2].grid(True, which='both')

    # === FASE ===
    axs[1,0].semilogx(freqs, YAA_ang, 'r-')
    axs[1,0].set_xlabel('Frequency (Hz)')
    axs[1,0].set_ylabel('Phase (Degrees)')
    axs[1,0].set_xlim([x_low_axis, x_up_axis])
    axs[1,0].grid(True, which='both')

    axs[1,1].semilogx(freqs, YBB_ang, 'r-')
    axs[1,1].set_xlabel('Frequency (Hz)')
    axs[1,1].set_xlim([x_low_axis, x_up_axis])
    axs[1,1].grid(True, which='both')

    axs[1,2].semilogx(freqs, YCC_ang, 'r-')
    axs[1,2].set_xlabel('Frequency (Hz)')
    axs[1,2].set_xlim([x_low_axis, x_up_axis])
    axs[1,2].grid(True, which='both')

    axs[0,2].legend(['ABC frequency scan'], loc='lower left')

    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()

    return fig


def dq0_plot_response(fsampling, Ym_qq, Ya_qq, Ym_qd, Ya_qd, Ym_dq, Ya_dq, Ym_dd, Ya_dd, fd0, Y_mag, Y_rad):
    # Convertir magnitudes a dB
    Ym_qq = 20 * np.log10(Ym_qq)
    Ym_qd = 20 * np.log10(Ym_qd)
    Ym_dq = 20 * np.log10(Ym_dq)
    Ym_dd = 20 * np.log10(Ym_dd)
    # Definir límites del eje X (frecuencia)
    x_low_axis = fd0[0]
    x_up_axis = fd0[-1]
    # Ajustar los límites del eje Y para magnitud y fase
    mag_min_qq = np.nanmin(Ym_qq)
    mag_max_qq = np.nanmax(Ym_qq)
    phase_min_qq = np.nanmin(Ya_qq)
    phase_max_qq = np.nanmax(Ya_qq)
    
    mag_min_dq = np.nanmin(Ym_dq)
    mag_max_dq = np.nanmax(Ym_dq)
    phase_min_dq = np.nanmin(Ya_dq)
    phase_max_dq = np.nanmax(Ya_dq)
    
    mag_min_qd = np.nanmin(Ym_qd)
    mag_max_qd = np.nanmax(Ym_qd)
    phase_min_qd = np.nanmin(Ya_qd)
    phase_max_qd = np.nanmax(Ya_qd)
    
    mag_min_dd = np.nanmin(Ym_dd)
    mag_max_dd = np.nanmax(Ym_dd)
    phase_min_dd = np.nanmin(Ya_dd)
    phase_max_dd = np.nanmax(Ya_dd)
    
    # Configuración de estilos globales para las gráficas (incluyendo tamaño de textos)
    plt.rcParams.update({
        'axes.titlesize': 20,       # Tamaño de los títulos de los ejes
        'axes.labelsize': 20,       # Tamaño de las etiquetas de los ejes (X, Y)
        'xtick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje X
        'ytick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje Y
        'legend.fontsize': 20,      # Tamaño del texto de las leyendas
        'figure.titlesize': 20,     # Tamaño del título de la figura (si lo hubiera)
        'lines.linewidth': 2        # Ancho de las líneas de las curvas
    })
    
    # Crear subplots
    fig, axs = plt.subplots(4, 2, figsize=(15, 10))
    # Top row - Magnitude in dB vs log(fsampling)
    axs[0, 0].semilogx(fsampling, Ym_qq, label='Linear state space', color='k')
    axs[0, 0].set_title(r'$Y_{qq}(s)$')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].set_xlim([x_low_axis, x_up_axis])
    axs[0, 0].set_ylim([mag_min_qq - 3, mag_max_qq + 3])
    axs[0, 0].grid(True, which='major', axis='both')
    axs[0, 0].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[1, 0].semilogx(fsampling, Ya_qq, label='Linear state space', color='k')
    axs[1, 0].set_ylabel('Phase (Degrees)')
    axs[1, 0].set_xlim([x_low_axis, x_up_axis])
    axs[1, 0].set_ylim([phase_min_qq - 3, phase_max_qq + 3])
    axs[1, 0].grid(True, which='major', axis='both')
    axs[1, 0].grid(True, which='minor', axis='x')
    
    axs[0, 1].semilogx(fsampling, Ym_qd, label='Linear state space', color='k')
    axs[0, 1].set_title(r'$Y_{qd}(s)$')
    axs[0, 1].set_xlim([x_low_axis, x_up_axis])
    axs[0, 1].set_ylim([mag_min_qd - 3, mag_max_qd + 3])
    axs[0, 1].grid(True, which='major', axis='both')
    axs[0, 1].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[1, 1].semilogx(fsampling, Ya_qd, label='Linear state space', color='k')
    axs[1, 1].set_xlim([x_low_axis, x_up_axis])
    axs[1, 1].set_ylim([phase_min_qd - 3, phase_max_qd + 3])
    axs[1, 1].grid(True, which='major', axis='both')
    axs[1, 1].grid(True, which='minor', axis='x')
    
    axs[2, 0].semilogx(fsampling, Ym_dq, label='Linear state space', color='k')
    axs[2, 0].set_title(r'$Y_{dq}(s)$')
    axs[2, 0].set_ylabel('Magnitude (dB)')
    axs[2, 0].set_xlim([x_low_axis, x_up_axis])
    axs[2, 0].set_ylim([mag_min_dq - 3, mag_max_dq + 3])
    axs[2, 0].grid(True, which='major', axis='both')
    axs[2, 0].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[3, 0].semilogx(fsampling, Ya_dq, label='Linear state space', color='k')
    axs[3, 0].set_ylabel('Phase (Degrees)')
    axs[3, 0].set_xlabel('Frequency (Hz)')
    axs[3, 0].set_xlim([x_low_axis, x_up_axis])
    axs[3, 0].set_ylim([phase_min_dq - 3, phase_max_dq + 3])
    axs[3, 0].grid(True, which='major', axis='both')
    axs[3, 0].grid(True, which='minor', axis='x')
    
    axs[2, 1].semilogx(fsampling, Ym_dd, label='Linear state space', color='k')
    axs[2, 1].set_title(r'$Y_{qd}(s)$')
    axs[2, 1].set_xlim([x_low_axis, x_up_axis])
    axs[2, 1].set_ylim([mag_min_dd - 3, mag_max_dd + 3])
    axs[2, 1].grid(True, which='major', axis='both')
    axs[2, 1].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[3, 1].semilogx(fsampling, Ya_dd, label='Linear state space', color='k')
    axs[3, 1].set_xlabel('Frequency (Hz)')
    axs[3, 1].set_xlim([x_low_axis, x_up_axis])
    axs[3, 1].set_ylim([phase_min_dd - 3, phase_max_dd + 3])
    axs[3, 1].grid(True, which='major', axis='both')
    axs[3, 1].grid(True, which='minor', axis='x')

    # Agregar puntos rojos para "dq0 frequency scan"
    frequencies = list(Y_mag.keys())
    for freq in frequencies:
        Y_qd0_mag = Y_mag[freq]
        Y_qd0_angle = Y_rad[freq]

        axs[0, 0].semilogx(freq, 20 * np.log10(Y_qd0_mag[0, 0]), 'ro', label='dq0 frequency scan')
        axs[0, 1].semilogx(freq, 20 * np.log10(Y_qd0_mag[0, 1]), 'ro', label='dq0 frequency scan')
        axs[2, 0].semilogx(freq, 20 * np.log10(Y_qd0_mag[1, 0]), 'ro', label='dq0 frequency scan')
        axs[2, 1].semilogx(freq, 20 * np.log10(Y_qd0_mag[1, 1]), 'ro', label='dq0 frequency scan')

        axs[1, 0].semilogx(freq, (180 / np.pi) * Y_qd0_angle[0, 0], 'ro', label='dq0 frequency scan')
        axs[1, 1].semilogx(freq, (180 / np.pi) * Y_qd0_angle[0, 1], 'ro', label='dq0 frequency scan')
        axs[3, 0].semilogx(freq, (180 / np.pi) * Y_qd0_angle[1, 0], 'ro', label='dq0 frequency scan')
        axs[3, 1].semilogx(freq, (180 / np.pi) * Y_qd0_angle[1, 1], 'ro', label='dq0 frequency scan')

    axs[0, 1].legend(['Linear state space', 'dq0 frequency scan'], loc='lower left')

    # Ajuste final
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()
    
def dq0_plot(fd0, Y_mag, Y_rad):
    
    # Definir límites del eje X (frecuencia)
    x_low_axis = fd0[0]
    x_up_axis = fd0[-1]

    # Configuración de estilos globales para las gráficas (incluyendo tamaño de textos)
    plt.rcParams.update({
        'axes.titlesize': 20,       # Tamaño de los títulos de los ejes
        'axes.labelsize': 20,       # Tamaño de las etiquetas de los ejes (X, Y)
        'xtick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje X
        'ytick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje Y
        'legend.fontsize': 20,      # Tamaño del texto de las leyendas
        'figure.titlesize': 20,     # Tamaño del título de la figura (si lo hubiera)
        'lines.linewidth': 2        # Ancho de las líneas de las curvas
    })
    
    # Crear subplots
    fig, axs = plt.subplots(4, 2, figsize=(15, 10))
    
    # Agregar puntos rojos para "dq0 frequency scan"
    frequencies = list(Y_mag.keys())

    # Create empty lists to store data for plotting outside the loop
    Yqq_mag, Yqq_phase, Yqd_mag, Yqd_phase = [], [], [], []
    Ydq_mag, Ydq_phase, Ydd_mag, Ydd_phase = [], [], [], []
    
    for freq in frequencies:
        Y_qd0_mag = Y_mag[freq]
        Y_qd0_angle = Y_rad[freq]

        # Collect the data
        Yqq_mag.append(20 * np.log10(Y_qd0_mag[0, 0]))
        Yqq_phase.append((180 / np.pi) * Y_qd0_angle[0, 0])
        Yqd_mag.append(20 * np.log10(Y_qd0_mag[0, 1]))
        Yqd_phase.append((180 / np.pi) * Y_qd0_angle[0, 1])
        Ydq_mag.append(20 * np.log10(Y_qd0_mag[1, 0]))
        Ydq_phase.append((180 / np.pi) * Y_qd0_angle[1, 0])
        Ydd_mag.append(20 * np.log10(Y_qd0_mag[1, 1]))
        Ydd_phase.append((180 / np.pi) * Y_qd0_angle[1, 1])

    # Now plot the full curves
    axs[0, 0].semilogx(frequencies, Yqq_mag, 'r-', label='dq0 frequency scan')
    axs[0, 0].set_title(r'$Y_{qq}(s)$')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].set_xlim([x_low_axis, x_up_axis])
    axs[0, 0].grid(True, which='major', axis='both')
    axs[0, 0].grid(True, which='minor', axis='x')

    axs[1, 0].semilogx(frequencies, Yqq_phase, 'r-', label='dq0 frequency scan')
    axs[1, 0].set_ylabel('Phase (Degrees)')
    axs[1, 0].set_xlim([x_low_axis, x_up_axis])
    axs[1, 0].grid(True, which='major', axis='both')
    axs[1, 0].grid(True, which='minor', axis='x')

    axs[0, 1].semilogx(frequencies, Yqd_mag, 'r-', label='dq0 frequency scan')
    axs[0, 1].set_title(r'$Y_{qd}(s)$')
    axs[0, 1].set_xlim([x_low_axis, x_up_axis])
    axs[0, 1].grid(True, which='major', axis='both')
    axs[0, 1].grid(True, which='minor', axis='x')

    axs[1, 1].semilogx(frequencies, Yqd_phase, 'r-', label='dq0 frequency scan')
    axs[1, 1].set_xlim([x_low_axis, x_up_axis])
    axs[1, 1].grid(True, which='major', axis='both')
    axs[1, 1].grid(True, which='minor', axis='x')

    axs[2, 0].semilogx(frequencies, Ydq_mag, 'r-', label='dq0 frequency scan')
    axs[2, 0].set_title(r'$Y_{dq}(s)$')
    axs[2, 0].set_ylabel('Magnitude (dB)')
    axs[2, 0].set_xlim([x_low_axis, x_up_axis])
    axs[2, 0].grid(True, which='major', axis='both')
    axs[2, 0].grid(True, which='minor', axis='x')

    axs[3, 0].semilogx(frequencies, Ydq_phase, 'r-', label='dq0 frequency scan')
    axs[3, 0].set_ylabel('Phase (Degrees)')
    axs[3, 0].set_xlabel('Frequency (Hz)')
    axs[3, 0].set_xlim([x_low_axis, x_up_axis])
    axs[3, 0].grid(True, which='major', axis='both')
    axs[3, 0].grid(True, which='minor', axis='x')

    axs[2, 1].semilogx(frequencies, Ydd_mag, 'r-', label='dq0 frequency scan')
    axs[2, 1].set_title(r'$Y_{dd}(s)$')
    axs[2, 1].set_xlim([x_low_axis, x_up_axis])
    axs[2, 1].grid(True, which='major', axis='both')
    axs[2, 1].grid(True, which='minor', axis='x')

    axs[3, 1].semilogx(frequencies, Ydd_phase, 'r-', label='dq0 frequency scan')
    axs[3, 1].set_xlabel('Frequency (Hz)')
    axs[3, 1].set_xlim([x_low_axis, x_up_axis])
    axs[3, 1].grid(True, which='major', axis='both')
    axs[3, 1].grid(True, which='minor', axis='x')

    axs[0, 1].legend(['dq0 frequency scan'], loc='lower left')

    # Ajuste final
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()
    
    return fig




def pn0_plot_response(fsampling, Ym_pp, Ya_pp, Ym_pn, Ya_pn, Ym_np, Ya_np, Ym_nn, Ya_nn, fd0, Y_mag, Y_rad):
    # Convertir magnitudes a dB
    Ym_pp = 20 * np.log10(Ym_pp)
    Ym_pn = 20 * np.log10(Ym_pn)
    Ym_np = 20 * np.log10(Ym_np)
    Ym_nn = 20 * np.log10(Ym_nn)
    # Definir límites del eje X (frecuencia)
    x_low_axis = fd0[0]
    x_up_axis = fd0[-1]
    # Ajustar los límites del eje Y para magnitud y fase
    mag_min_pp = np.nanmin(Ym_pp)
    mag_max_pp = np.nanmax(Ym_pp)
    phase_min_pp = np.nanmin(Ya_pp)
    phase_max_pp = np.nanmax(Ya_pp)
    
    mag_min_np = np.nanmin(Ym_np)
    mag_max_np = np.nanmax(Ym_np)
    phase_min_np = np.nanmin(Ya_np)
    phase_max_np = np.nanmax(Ya_np)
    
    mag_min_pn = np.nanmin(Ym_pn)
    mag_max_pn = np.nanmax(Ym_pn)
    phase_min_pn = np.nanmin(Ya_pn)
    phase_max_pn = np.nanmax(Ya_pn)
    
    mag_min_nn = np.nanmin(Ym_nn)
    mag_max_nn = np.nanmax(Ym_nn)
    phase_min_nn = np.nanmin(Ya_nn)
    phase_max_nn = np.nanmax(Ya_nn)
    
    # Configuración de estilos globales para las gráficas (incluyendo tamaño de textos)
    plt.rcParams.update({
        'axes.titlesize': 20,       # Tamaño de los títulos de los ejes
        'axes.labelsize': 20,       # Tamaño de las etiquetas de los ejes (X, Y)
        'xtick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje X
        'ytick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje Y
        'legend.fontsize': 20,      # Tamaño del texto de las leyendas
        'figure.titlesize': 20,     # Tamaño del título de la figura (si lo hubiera)
        'lines.linewidth': 2        # Ancho de las líneas de las curvas
    })
    
    # Crear subplots
    fig, axs = plt.subplots(4, 2, figsize=(15, 10))
    # Top row - Magnitude in dB vs log(fsampling)
    axs[0, 0].semilogx(fsampling, Ym_pp, label='Theoretical response', color='k')
    axs[0, 0].set_title(r'$Y_{pp}(s)$')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].set_xlim([x_low_axis, x_up_axis])
    axs[0, 0].set_ylim([mag_min_pp - 3, mag_max_pp + 3])
    axs[0, 0].grid(True, which='major', axis='both')
    axs[0, 0].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[1, 0].semilogx(fsampling, Ya_pp, label='Theoretical response', color='k')
    axs[1, 0].set_ylabel('Phase (Degrees)')
    axs[1, 0].set_xlim([x_low_axis, x_up_axis])
    axs[1, 0].set_ylim([phase_min_pp - 3, phase_max_pp + 3])
    axs[1, 0].grid(True, which='major', axis='both')
    axs[1, 0].grid(True, which='minor', axis='x')
    
    axs[0, 1].semilogx(fsampling, Ym_pn, label='Theoretical response', color='k')
    axs[0, 1].set_title(r'$Y_{pn}(s)$')
    axs[0, 1].set_xlim([x_low_axis, x_up_axis])
    axs[0, 1].set_ylim([mag_min_pn - 3, mag_max_pn + 3])
    axs[0, 1].grid(True, which='major', axis='both')
    axs[0, 1].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[1, 1].semilogx(fsampling, Ya_pn, label='Theoretical response', color='k')
    axs[1, 1].set_xlim([x_low_axis, x_up_axis])
    axs[1, 1].set_ylim([phase_min_pn - 3, phase_max_pn + 3])
    axs[1, 1].grid(True, which='major', axis='both')
    axs[1, 1].grid(True, which='minor', axis='x')
    
    axs[2, 0].semilogx(fsampling, Ym_np, label='Theoretical response', color='k')
    axs[2, 0].set_title(r'$Y_{np}(s)$')
    axs[2, 0].set_ylabel('Magnitude (dB)')
    axs[2, 0].set_xlim([x_low_axis, x_up_axis])
    axs[2, 0].set_ylim([mag_min_np - 3, mag_max_np + 3])
    axs[2, 0].grid(True, which='major', axis='both')
    axs[2, 0].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[3, 0].semilogx(fsampling, Ya_np, label='Theoretical response', color='k')
    axs[3, 0].set_ylabel('Phase (Degrees)')
    axs[3, 0].set_xlabel('Frequency (Hz)')
    axs[3, 0].set_xlim([x_low_axis, x_up_axis])
    axs[3, 0].set_ylim([phase_min_np - 3, phase_max_np + 3])
    axs[3, 0].grid(True, which='major', axis='both')
    axs[3, 0].grid(True, which='minor', axis='x')
    
    axs[2, 1].semilogx(fsampling, Ym_nn, label='Theoretical response', color='k')
    axs[2, 1].set_title(r'$Y_{nn}(s)$')
    axs[2, 1].set_xlim([x_low_axis, x_up_axis])
    axs[2, 1].set_ylim([mag_min_nn - 3, mag_max_nn + 3])
    axs[2, 1].grid(True, which='major', axis='both')
    axs[2, 1].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[3, 1].semilogx(fsampling, Ya_nn, label='Theoretical response', color='k')
    axs[3, 1].set_xlabel('Frequency (Hz)')
    axs[3, 1].set_xlim([x_low_axis, x_up_axis])
    axs[3, 1].set_ylim([phase_min_nn - 3, phase_max_nn + 3])
    axs[3, 1].grid(True, which='major', axis='both')
    axs[3, 1].grid(True, which='minor', axis='x')

    # Agregar puntos rojos para "pn frequency scan"
    frequencies = list(Y_mag.keys())
    for freq in frequencies:
        Y_pn0_mag = Y_mag[freq]
        Y_pn0_angle = Y_rad[freq]

        axs[0, 0].semilogx(freq, 20 * np.log10(Y_pn0_mag[0, 0]), 'ro', label='pn frequency scan')
        axs[0, 1].semilogx(freq, 20 * np.log10(Y_pn0_mag[0, 1]), 'ro', label='pn frequency scan')
        axs[2, 0].semilogx(freq, 20 * np.log10(Y_pn0_mag[1, 0]), 'ro', label='pn frequency scan')
        axs[2, 1].semilogx(freq, 20 * np.log10(Y_pn0_mag[1, 1]), 'ro', label='pn frequency scan')

        axs[1, 0].semilogx(freq, (180 / np.pi) * Y_pn0_angle[0, 0], 'ro', label='pn frequency scan')
        axs[1, 1].semilogx(freq, (180 / np.pi) * Y_pn0_angle[0, 1], 'ro', label='pn frequency scan')
        axs[3, 0].semilogx(freq, (180 / np.pi) * Y_pn0_angle[1, 0], 'ro', label='pn frequency scan')
        axs[3, 1].semilogx(freq, (180 / np.pi) * Y_pn0_angle[1, 1], 'ro', label='pn frequency scan')

    axs[0, 1].legend(['Theoretical response', 'pn frequency scan'], loc='lower left')

    # Ajuste final
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()
    
def pn0_plot(fd0, Y_mag, Y_rad):
    
    # Definir límites del eje X (frecuencia)
    x_low_axis = fd0[0]
    x_up_axis = fd0[-1]
    
    # Configuración de estilos globales para las gráficas (incluyendo tamaño de textos)
    plt.rcParams.update({
        'axes.titlesize': 20,       # Tamaño de los títulos de los ejes
        'axes.labelsize': 20,       # Tamaño de las etiquetas de los ejes (X, Y)
        'xtick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje X
        'ytick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje Y
        'legend.fontsize': 20,      # Tamaño del texto de las leyendas
        'figure.titlesize': 20,     # Tamaño del título de la figura (si lo hubiera)
        'lines.linewidth': 2        # Ancho de las líneas de las curvas
    })
    
    # Crear subplots
    fig, axs = plt.subplots(4, 2, figsize=(15, 10))
    
    # Preparar listas para almacenar datos de magnitud y fase
    frequencies = sorted(Y_mag.keys())
    Ypp_mag, Ypn_mag, Ynp_mag, Ynn_mag = [], [], [], []
    Ypp_phase, Ypn_phase, Ynp_phase, Ynn_phase = [], [], [], []

    # Acumular valores para cada frecuencia
    for freq in frequencies:
        Y_pn0_mag = Y_mag[freq]
        Y_pn0_angle = Y_rad[freq]
        
        # Magnitudes
        Ypp_mag.append(20 * np.log10(Y_pn0_mag[0, 0]) if Y_pn0_mag[0, 0] > 0 else None)
        Ypn_mag.append(20 * np.log10(Y_pn0_mag[0, 1]) if Y_pn0_mag[0, 1] > 0 else None)
        Ynp_mag.append(20 * np.log10(Y_pn0_mag[1, 0]) if Y_pn0_mag[1, 0] > 0 else None)
        Ynn_mag.append(20 * np.log10(Y_pn0_mag[1, 1]) if Y_pn0_mag[1, 1] > 0 else None)

        # Fases
        Ypp_phase.append((180 / np.pi) * Y_pn0_angle[0, 0])
        Ypn_phase.append((180 / np.pi) * Y_pn0_angle[0, 1])
        Ynp_phase.append((180 / np.pi) * Y_pn0_angle[1, 0])
        Ynn_phase.append((180 / np.pi) * Y_pn0_angle[1, 1])

    # Graficar las magnitudes y fases acumuladas
    axs[0, 0].semilogx(frequencies, Ypp_mag, 'r-', label='pn frequency scan')
    axs[0, 0].set_title(r'$Y_{pp}(s)$')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].set_xlim([x_low_axis, x_up_axis])
    axs[0, 0].grid(True, which='major', axis='both')
    axs[0, 0].grid(True, which='minor', axis='x')

    axs[1, 0].semilogx(frequencies, Ypp_phase, 'r-', label='pn frequency scan')
    axs[1, 0].set_ylabel('Phase (Degrees)')
    axs[1, 0].set_xlim([x_low_axis, x_up_axis])
    axs[1, 0].grid(True, which='major', axis='both')
    axs[1, 0].grid(True, which='minor', axis='x')

    axs[0, 1].semilogx(frequencies, Ypn_mag, 'r-', label='pn frequency scan')
    axs[0, 1].set_title(r'$Y_{pn}(s)$')
    axs[0, 1].set_xlim([x_low_axis, x_up_axis])
    axs[0, 1].grid(True, which='major', axis='both')
    axs[0, 1].grid(True, which='minor', axis='x')

    axs[1, 1].semilogx(frequencies, Ypn_phase, 'r-', label='pn frequency scan')
    axs[1, 1].set_xlim([x_low_axis, x_up_axis])
    axs[1, 1].grid(True, which='major', axis='both')
    axs[1, 1].grid(True, which='minor', axis='x')

    axs[2, 0].semilogx(frequencies, Ynp_mag, 'r-', label='pn frequency scan')
    axs[2, 0].set_title(r'$Y_{np}(s)$')
    axs[2, 0].set_ylabel('Magnitude (dB)')
    axs[2, 0].set_xlim([x_low_axis, x_up_axis])
    axs[2, 0].grid(True, which='major', axis='both')
    axs[2, 0].grid(True, which='minor', axis='x')

    axs[3, 0].semilogx(frequencies, Ynp_phase, 'r-', label='pn frequency scan')
    axs[3, 0].set_ylabel('Phase (Degrees)')
    axs[3, 0].set_xlabel('Frequency (Hz)')
    axs[3, 0].set_xlim([x_low_axis, x_up_axis])
    axs[3, 0].grid(True, which='major', axis='both')
    axs[3, 0].grid(True, which='minor', axis='x')

    axs[2, 1].semilogx(frequencies, Ynn_mag, 'r-', label='pn frequency scan')
    axs[2, 1].set_title(r'$Y_{nn}(s)$')
    axs[2, 1].set_xlim([x_low_axis, x_up_axis])
    axs[2, 1].grid(True, which='major', axis='both')
    axs[2, 1].grid(True, which='minor', axis='x')

    axs[3, 1].semilogx(frequencies, Ynn_phase, 'r-', label='pn frequency scan')
    axs[3, 1].set_xlabel('Frequency (Hz)')
    axs[3, 1].set_xlim([x_low_axis, x_up_axis])
    axs[3, 1].grid(True, which='major', axis='both')
    axs[3, 1].grid(True, which='minor', axis='x')

    axs[0, 1].legend(['pn frequency scan'], loc='lower left')

    # Ajuste final
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()
    
    return fig
