# Two-Sample Testing for Persistence Intensity Functions

This repository contains the code used in the paper

**“A Two-Sample Test on Weighted Persistence Intensity Functions in Topological Data Analysis.”**

The proposed method performs two-sample testing for persistence
intensity functions using kernel-based test statistics and permutation
calibration. The repository also contains implementations of competing
methods and scripts for the simulation and real-data experiments presented
in the paper.

## Repository Structure

```text
.
├── source/                     # Core functions for test methods
├── simulations/
│   ├── Circles/               
│   ├── Instrument/         
│   ├── orbit5k/              
│   └── Torus/
├── Per_Image_testing/
│   ├── PI_orbit5k_data/
│   ├── PI_circle_data/
│   └── PI_Torus_data/
├── results/
│   ├── figures/
│   │   └── figure_file/
│   └── simulation_results/
│       ├── Circle_results/
│       ├── Torus_results/
│       └── Orbit5k_results/
├── data/
│   ├── clarinet_diagrams/                
│   ├── flute_diagrams/          
│   ├── json/ 
│   ├── orbit_data/
│   └── Rdata
├── Figures.ipynb               # Notebook for generating figures
└── README.md
```

## Folder Description

### `source/`

This directory contains Python modules implementing the proposed bandwidth-aggregated two-sample test and other competing methods. 
The implementation of the aggregation procedure is adapted from the MMDAgg code of Schrab et al. (2023), with modifications for persistence intensity functions.


### `simulations/`

Contains all scripts for the simulation studies.

- `Circles/` – Circles simulation experiment.
- `Instrument/` – Musical instrument experiment (flute vs. clarinet).
- `orbit5k/` – ORBIT5K simulation experiment.
- `Torus/` – Torus simulation experiment.

### `Per_Image_testing/`

Implementation of the persistence image two-sample test (Moon and Lazar, 2023).
The implementation of the Two-stage two-sample testing method is adapted from the code of Moon and Lazar. (2023).
This directory contains the implementation of the two-stage testing method based on persistence images.

### `results/`

Stores processed experimental outputs.

### `data/`

Data used in the experiments.
This directory contains source code for generating the data used in the paper rather than the generated data themselves

### `Figures.ipynb`

Jupyter notebook used to reproduce all figures in the paper.