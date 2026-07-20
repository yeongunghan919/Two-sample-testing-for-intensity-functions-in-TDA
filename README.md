# Two-Sample Testing for Persistence Diagram Intensity Functions

This repository contains the code used in the paper

**“A Two-Sample Test on Weighted Persistence Intensity Functions in Topological Data Analysis.”**

The proposed method performs two-sample testing for persistence diagram
intensity functions using kernel-based test statistics and permutation
calibration. The repository also contains implementations of competing
methods and scripts for the simulation and real-data experiments presented
in the paper.

## Repository Structure

```text
.
├── source/                  # Core functions for the proposed test
├── simulations/             # Simulation experiments
├── Per_Image_testing/       # Persistence image-based competing method
├── Tikz_intensity_testing/  # Code for figures and intensity visualizations
├── results/                 # Processed experimental results
├── Figures.ipynb            # Notebook for generating figures
└── README.md