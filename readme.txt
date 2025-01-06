Program Overview
The GSTA (Generalized Smoothed Trajectory Analysis) plugin for PLUMED provides a flexible framework for 
the convolution of molecular dynamics trajectories using a weighted moving average. Users can specify custom kernel 
functions to define the weights, enabling a tailored approach to trajectory analysis.

Features
Filters trajectory data using a kernel-defined weighted moving average.
Supports arbitrary kernel shapes provided in an external file (kernel.dat).
Outputs results in XYZ format for easy visualization.
Configurable options for input/output file names and atom selections.

Installation
Ensure PLUMED is installed on your system. For installation instructions, refer to the PLUMED website.
Place the GSTA plugin file (GSTA.cpp) in the appropriate PLUMED src/analysis directory.
Rebuild PLUMED. The code, GSTA.cpp was tested with the version of plumed2-2.9.2.

Configuration
Prepare the plumed.dat configuration file. Below is an example configuration:

GSTA ...
  ATOMS=@allatoms
  OUTPUT_STRIDE=10
  OUTPUT_FILE=filt_water_traj.xyz
  SYMBOL_FILE=water_traj.xyz
  KERNEL_FILE=kernel.dat
... GSTA

PRINT ARG=@0.dummy FILE=/dev/null

Explanation of Parameters:
ATOMS=@allatoms: Processes all atoms in the input file. Specify a subset, e.g., 1-100, to limit the selection.
OUTPUT_STRIDE=10: Writes smoothed data every 10 steps.
OUTPUT_FILE=filt_water_traj.xyz: Output file for the filtered trajectory.
SYMBOL_FILE=water_traj.xyz: Input file to read atomic symbols.
KERNEL_FILE=kernel.dat: File containing kernel weights for smoothing.

Running the Program
Use the following command to run the GSTA plugin via the PLUMED driver:
plumed driver --plumed plumed.dat --box 0.0,0.0,0.0 --ixyz water_traj.xyz
Explanation:
plumed driver: Executes PLUMED in standalone mode to process trajectory data.
--plumed plumed.dat: Specifies the configuration file containing the GSTA settings.
--box 0.0,0.0,0.0: Indicates no periodic boundary conditions. Set the box size if periodic boundaries are required.
--ixyz water_traj.xyz: Specifies the input trajectory file in XYZ format.

Output
The program generates the following files:
filt_water_traj.xyz: Filtered trajectory in XYZ format.
First line: Number of atoms.
Second line: Comment or metadata (e.g., time step).
Subsequent lines: Atomic symbols and averaged coordinates.

Contact
For questions or feedback, please contact: Dr. Ádám Madarász
Email: madarasz.adam@ttk.hu

