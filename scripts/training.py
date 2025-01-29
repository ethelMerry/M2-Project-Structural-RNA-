import os
import numpy as np
from Bio.PDB import PDBParser
import math

# Step 1: Specify the input folder and output folder
input_folder = "Final/data"  # Replace with the folder containing your PDB files
output_folder = "Final/output"  # Replace with the folder where you want to save results

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Define the 10 base pairs to track
base_pairs = [('A', 'A'), ('A', 'U'), ('A', 'C'), ('A', 'G'),
              ('U', 'U'), ('U', 'C'), ('U', 'G'),
              ('C', 'C'), ('C', 'G'), ('G', 'G')]

# Define 20 distance bins from 0 to 20 Å
bins = np.linspace(0, 20, 21)

# Initialize dictionaries to aggregate distances
pair_distances = {pair: [] for pair in base_pairs}
xx_distances = []  # For the reference frequency (XX pair)

# Step 2: Process each PDB file and aggregate distances
for pdb_filename in os.listdir(input_folder):
    if pdb_filename.endswith(".pdb"):
        pdb_path = os.path.join(input_folder, pdb_filename)
        try:
            # Parse the PDB file
            parser = PDBParser()
            structure = parser.get_structure(pdb_filename, pdb_path)

            C3_atoms = []
            # Extract C3' atom coordinates and their corresponding residues
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.resname in ['A', 'U', 'C', 'G']:
                            for atom in residue:
                                if atom.get_name() == "C3'":
                                    C3_atoms.append((residue.resname, atom))

            # Compute pairwise distances (intrachain, separated by at least 3 positions)
            for i in range(len(C3_atoms)):
                res_i, atom_i = C3_atoms[i]
                for j in range(i + 4, len(C3_atoms)):  # Only consider i and i+4, i+5, etc.
                    res_j, atom_j = C3_atoms[j]
                    
                    # Calculate distance between C3' atoms
                    coord_i = atom_i.get_coord()
                    coord_j = atom_j.get_coord()
                    distance = np.linalg.norm(coord_i - coord_j)

                    # Only consider intrachain distances
                    if atom_i.get_parent().get_full_id()[2] == atom_j.get_parent().get_full_id()[2]:
                        # Check if the pair is one of the base pairs we are tracking
                        if (res_i, res_j) in pair_distances:
                            pair_distances[(res_i, res_j)].append(distance)
                        elif (res_j, res_i) in pair_distances:
                            pair_distances[(res_j, res_i)].append(distance)
                        
                        # Add to the reference (XX) pair regardless of the nucleotides
                        xx_distances.append(distance)
        except Exception as e:
            print(f"Error processing {pdb_filename}: {e}")

# Step 3: Compute observed and reference frequencies
observed_frequencies = {pair: np.zeros(len(bins) - 1) for pair in base_pairs}
reference_frequencies = np.zeros(len(bins) - 1)

# Compute frequencies for each base pair
for pair, distances in pair_distances.items():
    observed_frequencies[pair], _ = np.histogram(distances, bins=bins)

# Compute reference frequencies (XX pair)
reference_frequencies, _ = np.histogram(xx_distances, bins=bins)

# Normalize frequencies (observed and reference)
for pair in base_pairs:
    observed_frequencies[pair] = observed_frequencies[pair].astype(float)  # Convert to float
    observed_frequencies[pair] /= sum(observed_frequencies[pair])

# Normalize reference frequencies
reference_frequencies = reference_frequencies.astype(float)  # Convert to float
reference_frequencies /= sum(reference_frequencies)

# Step 4: Compute the log-ratio of the observed to reference frequencies
pseudo_energies = {pair: np.zeros(len(bins) - 1) for pair in base_pairs}

for pair in base_pairs:
    for r in range(len(bins) - 1):
        obs_freq = observed_frequencies[pair][r]
        ref_freq = reference_frequencies[r]

        # Avoid division by zero or log(0) by setting a small epsilon value
        if obs_freq > 0 and ref_freq > 0:
            pseudo_energies[pair][r] = -math.log(obs_freq / ref_freq)
        else:
            pseudo_energies[pair][r] = 10  # Cap the value at 10

# Step 5: Write the results (10 files, 1 per base pair)
for pair in base_pairs:
    pair_name = f"{pair[0]}{pair[1]}"
    output_file = os.path.join(output_folder, f"{pair_name}_scores.txt")
    with open(output_file, 'w') as file:
        for r, energy in enumerate(pseudo_energies[pair]):
            file.write(f"Distance bin {bins[r]}-{bins[r+1]} Å: {energy:.3f}\n")

print("Pseudo-energies have been computed and saved.")
