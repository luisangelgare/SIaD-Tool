
# SIaD Tool

Welcome to the **Stability and Interactions assessment in the frequency-Domain (SIaD)** tool for modern power system applications. This free-access, multi-platform tool is implemented in **MATLAB/Simulink** and **Python/PSCAD** environments.

The **SIaD Tool** has been developed as part of the **MSCA-ADOreD** project, funded by the European Union’s Horizon Europe Research and Innovation Programme under the **Marie Skłodowska-Curie grant agreement No. 101073554**.

---

## 🚨 Important Notice
This tool is **experimental** and under **active development**. We encourage users to report any issues and contribute to its improvement.

---

## 📋 Requirements

### For MATLAB/Simulink Version:
- **MATLAB**: Version 2022b or newer.

### For Python/PSCAD Version:
- **Python**: Version 3.12.7.
- **PSCAD/EMTDC**: Version 5.0.2.

### Recommended Hardware (For Optimal Performance):
- **Processor**: Intel Core i5 or higher.
- **RAM**: 8 GB or more.

---

## ⚙️ Installation

### For Python/PSCAD:

1. Clone the repository:
   git clone https://github.com/luisangelgare/FD-Scanning-Tool.git

2. Copy and use the "Frequency Domain Scanner" component in your PSCAD project:
   - Open your PSCAD canvas.
   - Use `Paste Special > Paste Transfer` to insert the component into your main project.

3. Place the `FDScanningTool.py` file in the same directory as your PSCAD project and follow the parameter and initialization instructions provided in the file.

4. Refer to the examples included in the repository for detailed usage instructions.
5. For the moment, just voltage perturbation and single-tone strategy is available in PSCAD.

### For MATLAB/Simulink:

1. Clone the repository:
   git clone https://github.com/luisangelgare/FD-Scanning-Tool.git

2. Copy the "Frequency Domain Scanner" block to your Simulink workspace:
   - Add the `powergui` block to your Simulink model.
   - Configure the step time and simulation time in the **Configuration Parameters**.
   - Set the solver to **Ode1**.
   - Enable **Accelerator Mode** in the simulation settings.

3. Copy and follow the parameter settings of the `exe_FDScanningTool.m` program into your project directory.

4. Refer to the **Manual_V3_SIaD-Tool.pdf** for the details about the configuration and installation process.

### Examples will be available soon!

## 🚨 Important Advice
 Recommended to Use first the MATLAB version in order to understand the operating steps of the tool.
