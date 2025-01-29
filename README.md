# Project - RNA Folding Objective Function

## Project Overview
The goal of this project is to infer the Gibbs free energy of a pure RNA structure using an objective function trained on interatomic distance distributions. The function uses known experimentally determined 3D structures to compute interatomic distances of C3' atoms of nucleotides and evaluates the stability of RNA structures.

## Getting Started  
The dataset consists of 100 RNA structures downloaded from the [RCSB website](https://www.rcsb.org/stats/growth/growth-rna), filtered based on the following criteria:

- **Polymer Entity Type**: RNA
- **Experimental Method**: X-ray Diffraction
- **Polymer Entity Sequence Length**: >= 2000
- **Number of Chains**: < 3 (for greater number of intrachain pairs)

### RNA Pairings and Distance Binning
The objective function considers 10 different RNA pairings:  
- AA, AC, AG, AU, CC, CG, CU, GG, GU, UU

These pairings are classified into 20 distance bins, each 1 Å wide, ranging from 0 Å to 20 Å (excluded). The choice of 20 Å as the upper limit is based on the fact that interatomic interactions weaken significantly at greater distances and have minimal impact beyond this threshold.
Note: The maximum scoring value is arbitrarily set to 10 to prevent overly high penalties and maintain a balanced evaluation of RNA structure stability.

## Training Script
The first part of the project focuses on training the objective function using a subset of the downloaded RNA structures.

- **Testing Phase**: Initially, the script was tested with a random selection of 5 RNA structures (1ffk.pdb, 1j5a.pdb, 1jj2.pdb, 1jzx.pdb, 1jzy.pdb) from the 100 downloaded structures. Pseudo scores, plots, and Gibbs energy scores were obtained and saved in the `testing` folder.
- **Full Training**: After successfully testing with 5 structures, the script was applied to the full dataset of 100 RNA structures. The results were saved in the `main_results` folder.

## Folder Structure
- `testing/`: Contains the `test_data`, `RNA_Puzzles_data`, `test_output` (which are the pseudoenergy scores), and plot scoring profiles as `5_Interaction_Profiles.pdf`, and the Gibbs energy scores as `test_gibbs_score.tsv`.
- `main_results/`: Contains all the results of the training phase (`output`, `100_Interaction_Profiles.pdf`, `gibbs_score.tsv`) with 100 RNA structures.

## Scripts
### 1. Training Script (`training.py`)
This script processes experimentally determined RNA structures to extract interatomic distance distributions and compute a pseudo-energy function for RNA stability evaluation.

**Description:**  
- Extracts C3' atom coordinates from RNA structures and computes interatomic distances.  
- Constructs interatomic distance distributions for 10 RNA base pairings.  
- Computes pseudo-energy scores based on observed distance frequencies.  
- Saves the trained objective function as scoring files for further analysis.  

### 2. Plotting Script (`PLOT.rmd`)

**Description:**  
- Reads pseudo-energy score files for 10 RNA base pairings.  
- Extracts and processes interaction scores from text files.  
- Plots scoring profiles as a function of interatomic distance bins.  
- Saves the final visualization as a PDF (`5_Interaction_Profiles.pdf` and `100_Interaction_Profiles.pdf`).  

### 3. Evaluation Script (`scores.py`)
**Description:**  
This script evaluates the predicted RNA structures using the trained objective function with [RNA Puzzles dataset](https://github.com/RNA-Puzzles/raw_dataset_and_for_assessment).
- Computes interatomic distances for a given RNA structure using the same thresholds (20 Å and i, i+4).  
- Assigns a pseudo-energy score to each distance using linear interpolation.  
- Sums all scores to estimate the Gibbs free energy of the RNA conformation.  
- Outputs the final stability score of the evaluated RNA structure.  

## Installation and Requirements
For this project, I used both **Python** and **R** for different stages of analysis.  

R Dependencies for the script `PLOT.rmd`:
- `install.packages("tidyverse")`  # Includes ggplot2 and dplyr for data manipulation and visualization
- `library(ggplot2)`
- `library(tidyverse)`

To run the other two scripts, you will need the following python packages:
- `numpy`
- `matplotlib`
- `scipy`
- `biopython`
- `pandas`

You can also install them using pip:
```bash
pip install numpy matplotlib scipy biopython pandas
