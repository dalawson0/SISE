# SISE - Simultaneous Input and State Estmation for MATLAB

**Author:** Daniel Lawson
**Version:** v0.0.0
**MATLAB Release:** R2025b
**License:** MIT

---

## Overview
**SISE** is a MATLAB toolbox providing algorithms for **Simultaneous Input and State Estimation (SISE)** in linear systems. The toolbox is designed to support research, simulation, and practical applications in state estimation and controls.

* **Core estimation function:** `ULISE.m` — implements the Unified Linear Input and State Estimator for discrete-time, time-invariant systems with one or zero timestep delay. 
* **Live demonstration:** `example_ULISE_fault_identification.mlx` — an interactive example showcasing the ULISE algorithm on simulated data. 
* **Supporting resources:** documentation, example data, and guides to help you get started with the toolbox. 

Future versions will expand the toolbox with additional estimation algorithms and potential Simulink models.

___

## Installation

The SISE Toolbox can be installed in MATLAB once the `.mltbx` release is available:

**Future Release (Recommended): Download MATLAB Toolbox File**
A `.mltbx` installer will be provided in the GitHub releases. Allowing for a one-click installation.

**Current Option: GitHub**
1. Download and extract the ZIP file or clone the repository:
	```bash
	git clone https://github.com/dalawson0/SISE.git
	```
2.Add the SISE folder (and all subfolders) to your MATLAB path:
	```matlab
	addpath(genpath("path/to/SISE_folder"));
	```
