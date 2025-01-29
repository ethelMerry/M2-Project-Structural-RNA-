# Project - RNA Folding Objective Function

## Project Overview
The goal of this project is to infer the Gibbs free energy of a pure RNA structure using an objective function trained on interatomic distance distributions. The function uses known experimentally determined 3D structures to compute interatomic distances of C3' atoms of nucleotides and evaluates the stability of RNA structures.

## Getting Started  
The dataset consists of 100 RNA structures-  (1FFK,1J5A,1JJ2,1JZX,1JZY,1JZZ,1K01,1K73,1K8A,1K9M,1KC8,1KD1,1KQS,1M1K,1M90,1N8R,1NJI,1NJM,1NJN,1NJO,1NJP,1NKW,1NWX,1NWY,1OND,1P9X,1Q7Y,1Q81,1Q82,1Q86,1QVF,1QVG,1S72,1SM1,1VQ4,1VQ5,1VQ6,1VQ7,1VQ8,1VQ9,1VQK,1VQL,1VQM,1VQN,1VQO,1VQP,1VVJ,1VY4,1VY5,1VY6,1VY7,1W2B,1XBP,1Y69,1YHQ,1YI2,1YIJ,1YIT,1YJ9,1YJN,1YJW,1Z58,2AAR,2D3O,2O43,2O44,2O45,2OGM,2OGN,2OGO,2OTJ,2OTL,2QA4,2QEX,2ZJP,2ZJQ,2ZJR,3CC2,3CC4,3CC7,3CCE,3CCJ,3CCL,3CCM,3CCQ,3CCR,3CCS,3CCU,3CCV,3CD6,3CF5,3CMA,3CME,3CPW,3CXC,3DLL,3FWO,3G4S,3G6E,3G71,3I55,3I56,3JQ4,3OW2,3PIO), downloaded from the [RCSB website](https://www.rcsb.org/stats/growth/growth-rna), filtered based on the following criteria:

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

- **Testing Phase**: Initially, the script was tested with a random selection of 5 RNA structures (1ffk.pdb, 1j5a.pdb, 1jj2.pdb, 1jzx.pdb, 1jzy.pdb) from the 100 downloaded structures. Pseudo scores, plots, and Gibbs energy scores were obtained and saved in the [testing](testing/) folder.
- **Full Training**: After successfully testing with 5 structures, the script was applied to the full dataset of 100 RNA structures. The results were saved in the [Main_Results](Main_Results/) folder.

## Folder Structure
### [testing_folder](testing/): 
- [test_data](testing/test_data): The 5 data initially used to test the script
- [RNA_Puzzles_data](testing/RNA_Puzzles_data): To evaluate the estimated Gibbs free energy of the evaluated RNA conformation
- [test_output](testing/test_output): Computed pseudoenergy scores of 5 RNA structure
- [5_Interaction_Profiles.pdf](testing/5_Interaction_Profiles.pdf) : Interaction plot (Distance vs Score)
- [test_gibbs_free_energy_results.tsv](testing/test_gibbs_free_energy_results.tsv): Estimated Gibbs energy scores

### [Main_Results](Main_Results/): 
- [output_pseudoscores](Main_Results/output_pseudoscores): Computed pseudoenergy scores of 100 RNA structures.
- [100_Interaction_Profiles.pdf](Main_Results/100_Interaction_Profiles.pdf): Interaction plot (Distance vs Score)
- [100_test_gibbs_free_energy_results.tsv](Main_Results/100_test_gibbs_free_energy_results.tsv): Estimated Gibbs energy scores

## Scripts
### 1. Training Script [training.py](scripts/training.py)
This script processes experimentally determined RNA structures to extract interatomic distance distributions and compute a pseudo-energy function for RNA stability evaluation.

**Description:**  
- Extracts C3' atom coordinates from RNA structures and computes interatomic distances.  
- Constructs interatomic distance distributions for 10 RNA base pairings.  
- Computes pseudo-energy scores based on observed distance frequencies.  
- Saves the trained objective function as scoring files for further analysis.  

### 2. Plotting Script [PLOT.Rmd](scripts/PLOT.Rmd)

**Description:**  
- Reads pseudo-energy score files for 10 RNA base pairings.  
- Extracts and processes interaction scores from text files.  
- Plots scoring profiles as a function of interatomic distance bins.  
- Saves the final visualization as a PDF (`5_Interaction_Profiles.pdf` and `100_Interaction_Profiles.pdf`).  

### 3. Evaluation Script [scores.py](scripts/scores.py)
**Description:**  
This script evaluates the predicted RNA structures using the trained objective function with [RNA Puzzles dataset](https://github.com/RNA-Puzzles/raw_dataset_and_for_assessment).
- Computes interatomic distances for a given RNA structure using the same thresholds (20 Å and i, i+4).  
- Assigns a pseudo-energy score to each distance using linear interpolation.  
- Sums all scores to estimate the Gibbs free energy of the RNA conformation.  
- Outputs the final stability score of the evaluated RNA structure.  

## Installation
To install and set up this project, follow the steps below:

1. Clone the repository to your local machine:
   ```bash
   git clone https://github.com/ethelMerry/M2-Project-Structural-RNA-
   cd M2-Project-Structural-RNA-
   
## Requirements
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
