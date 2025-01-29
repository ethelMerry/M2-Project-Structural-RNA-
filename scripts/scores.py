import os
import math
import numpy as np
from Bio.PDB import PDBParser

# 1) LOAD PRE-TRAINED SCORES FROM .TXT FILES
def load_pairwise_scores(score_dir):
    """
    Reads the 10 text files (one per base-pair type) from 'score_dir' and
    returns a dictionary {('A','A'): [score_bin0, score_bin1, ... , score_bin19], ...}
    """
    # The 10 base pairs, matching your training code
    base_pairs = [
        ('A','A'), ('A','U'), ('A','C'), ('A','G'),
        ('U','U'), ('U','C'), ('U','G'),
        ('C','C'), ('C','G'),
        ('G','G')
    ]
    
    pair_scores = {}
    for bp in base_pairs:
        pair_name = f"{bp[0]}{bp[1]}"
        file_name = os.path.join(score_dir, f"{pair_name}_scores.txt")
        
        scores = []
        with open(file_name, 'r') as f:
            for line in f:
                # e.g. line = "Bin 0.0-1.0 Å: 10.000"
                # We can parse out the last token as the score
                parts = line.strip().split()
                # parts[-1] should be something like "10.000"
                score_val = float(parts[-1])
                scores.append(score_val)
        
        # We should have 20 scores if each file has 20 bins
        if len(scores) != 20:
            raise ValueError(f"File {file_name} does not have 20 lines of scores.")
        
        # Store them in the dictionary
        pair_scores[bp] = scores
    
    return pair_scores

# 2) DISTANCE-BASED SCORING VIA LINEAR INTERPOLATION
def interpolate_score(distance, scores):
    """
    Given a distance (0 <= distance < 20) and a list of 20 pseudo-energy scores
    (one for each integer bin [0,1), [1,2), ... [19,20)),
    do linear interpolation.
    
    For example, if distance = 2.3:
      - lower_idx = 2
      - upper_idx = 3
      - fraction = 0.3
      final_score = scores[2] + 0.3 * (scores[3] - scores[2])

    Edge cases:
      - If distance < 1, lower_idx=0, upper_idx=1
      - If distance ~ 19.7, lower_idx=19, upper_idx=19, fraction=0 => scores[19]
    """
    if distance < 0:
        return 0.0
    
    # Clamp to [0, 19.999...) for safety
    if distance >= 20:
        distance = 19.9999999
    
    lower_idx = int(distance)  # floor
    # For the upper index, if distance is 19.2 => upper_idx=20, but we only have scores[0..19]
    # So we clamp at 19
    if lower_idx == 19:
        # no room for interpolation above 19
        return scores[19]
    
    upper_idx = lower_idx + 1
    fraction = distance - lower_idx
    
    lower_score = scores[lower_idx]
    upper_score = scores[upper_idx]
    
    # Linear interpolation
    interpol_val = lower_score + fraction * (upper_score - lower_score)
    return interpol_val


# 3) PARSE A PDB FILE, COMPUTE THE SUM OF INTERPOLATED SCORES
def parse_c3_atoms(pdb_file):
    """
    Parses a single PDB file and returns a list of (res_name, atom_object)
    for all C3' atoms in standard RNA residues (A, U, C, G).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('RNA', pdb_file)
    
    c3_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname in ['A', 'U', 'C', 'G']:
                    for atom in residue:
                        if atom.get_name() == "C3'":
                            c3_list.append((residue.resname, atom))
    return c3_list

def compute_pseudoenergy_for_structure(pdb_file, pair_scores):
    """
    For a single PDB file, compute the total 'pseudoenergy':
      1. Extract C3' atoms for each A/U/C/G residue
      2. For each pair i, j where j >= i+4 and in the same chain, distance < 20 Å:
         - Identify base pair (A,A), (A,U), (A,C), etc.
         - Look up the 20-bin scores from `pair_scores`
         - Linearly interpolate the score for that specific distance
         - Add it to the running total
    Returns the total sum.
    """
    # Step 1: Get all relevant atoms
    c3_atoms = parse_c3_atoms(pdb_file)
    
    total_score = 0.0
    
    # Step 2: Loop over all pairs (i+4)
    n = len(c3_atoms)
    for i in range(n):
        res_i, atom_i = c3_atoms[i]
        coord_i = atom_i.get_coord()
        chain_i = atom_i.get_parent().get_parent().id
        
        for j in range(i+4, n):
            res_j, atom_j = c3_atoms[j]
            chain_j = atom_j.get_parent().get_parent().id
            
            # Same chain?
            if chain_i == chain_j:
                coord_j = atom_j.get_coord()
                dist = np.linalg.norm(coord_i - coord_j)
                
                if dist < 20.0:
                    # Determine pair type. If (res_i, res_j) is not in pair_scores,
                    # we try (res_j, res_i).
                    # This step is needed because we only stored symmetrical
                    # pairs e.g. (A,U) but not necessarily (U,A) in the dictionary
                    pair = (res_i, res_j)
                    if pair not in pair_scores:
                        pair_rev = (res_j, res_i)
                        if pair_rev in pair_scores:
                            pair = pair_rev
                        else:
                            continue  # skip if it's truly not in the dictionary
                    # Interpolate
                    score_contribution = interpolate_score(dist, pair_scores[pair])
                    total_score += score_contribution
    
    return total_score


# 4) MAIN: LOAD SCORES, SCORE EACH PDB, PRINT OR SAVE
if __name__ == "__main__":
    # 1. Folder with the 10 score files
    score_folder = "Final/output_new_pdbs"  # Path to the precomputed scores (adjust to your path)
    bp_scores = load_pairwise_scores(score_folder)  # Load the base-pair score tables

    # 2. Main folder with PDB files downloaded from RNA PUZZLES (here I used PZ9)
    main_pdb_dir = "Final/RnaPuzzles_pdbs/PZ9"  

    # 3. Output folder for results
    output_dir = "Final/Res_gibbsnew_pdb" 
    os.makedirs(output_dir, exist_ok=True)  # Ensure the output folder exists

    # Collect all PDB files in the main directory (no subfolders)
    pdb_files = [f for f in os.listdir(main_pdb_dir) if f.endswith(".pdb")]
    if not pdb_files:
        print(f"No PDB files found in {main_pdb_dir}")
    
    # Sort PDB files for consistent output
    pdb_files.sort()

    # Output file for results
    results_file = os.path.join(output_dir, "new_gibbs_free_energy_results.tsv")
    with open(results_file, 'w') as out:
        out.write("PDB_File\tPseudoenergy\n")

        # Process each PDB file
        for pdb_file in pdb_files:
            pdb_path = os.path.join(main_pdb_dir, pdb_file)
            pseudoenergy_val = compute_pseudoenergy_for_structure(pdb_path, bp_scores)

            # Log the result
            print(f"{pdb_file} => {pseudoenergy_val:.4f}")
            out.write(f"{pdb_file}\t{pseudoenergy_val:.4f}\n")

    print(f"Scoring complete. Results saved to {results_file}")


