# Project - RNA Folding Objective Function

## Project Overview

This project aims to create an objective function for the RNA folding problem, where the goal is to identify the native fold of a given ribonucleotide chain. The native fold corresponds to the structure with the lowest Gibbs free energy. To estimate this energy, we will develop an objective function based on interatomic distance distributions from experimentally determined RNA 3D structures.

## Implementation Overview

The project consists of three Python scripts:

1. **Training Script**: Computes interatomic distances from a dataset of PDB files and derives scoring functions based on base-pair distance distributions.
2. **Plotting Script**: Visualizes the scoring profiles, showing how the estimated Gibbs free energy varies with interatomic distances.
3. **Scoring Script**: Evaluates RNA structures from the RNA-Puzzles dataset using the trained objective function.

## Getting Started

### Prerequisites

Before running the scripts, ensure you have the following dependencies installed:

- Python 3.x
- NumPy
- Pandas
- Biopython
- Matplotlib (for visualization)
- R (for advanced plotting, if applicable)
- ggplot2 (if using R for plotting)

You can install the necessary Python packages using:

```bash
pip install numpy pandas biopython matplotlib
