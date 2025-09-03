from itertools import groupby
import matplotlib.pyplot as plt
from natsort import natsorted
import numpy as np
import pandas as pd
from pyDOE import lhs
from sklearn.decomposition import PCA

# Define lists of row and column names
n_384_rows, n_384_cols = 16, 24
n_96_rows, n_96_cols = 8, 12
rows_384 = list(map(chr, range(65, 65+n_384_rows)))
cols_384 = list(map(str, range(1, 25)))
rows_96 = list(map(chr, range(65, 65+n_96_rows)))
cols_96 = list(map(str, range(1, 13)))

# Converts list of 96-well positions to list of 384-well positions
# Shifts column of first entry by number of columns specified in offset variable
# NOTE: row offset not yet functional
# inputs: reps (number of sample replicates), offset_384 (column at which 384-well list begins)

# mode 'alphabetical': 1D list for (DLS), 
# mode 'per_well': list of {number of repetitions}-element lists (which 384 wells map to which 96 wells)
# mode 'per_step': list of {number of repetitions} lists (e.g. Hamilton step)

def map_96_to_384(positions_96, reps, mode, rows_96, cols_96, offset_384):
  

  full_cols = [sum([int(i[1:])==int(n) for i in positions_96])for n in cols_96]
  col_groups = [[key, len(list(group))] for key, group in groupby(full_cols)]
  if reps == 4:
    multiplier = 2
  elif reps == 8:
    multiplier = 4
  elif reps == 2:
    multiplier = 1

  ic = 0
  positions_384 = []
  for g in col_groups:
    for i in range(offset_384+ic*multiplier, offset_384+ic*multiplier+g[1]):
      for j in range(0, g[0]*2,2):
        indices = []
        for k in range(0,2):
          for l in range(0,multiplier):
            indices += [rows_384[j+k] + cols_384[i+l*g[1]]]
        positions_384.append(indices)
    ic = ic+g[1]

  organized_positions = []
  if mode=='per_well':
    return positions_384
  elif mode=='per_step':
    for i in range(len(indices)):
      organized_positions.append([x[i] for x in positions_384])
  elif mode=='alphabetical':
    for i in range(len(positions_384)):
      for j in range(len(indices)):
        organized_positions.append(positions_384[i][j])
    organized_positions = natsorted(organized_positions)

  return organized_positions

# Converts list of 384-well positions to list of 96-well positions
# Shifts column of first entry by number of columns specified in offset variable
# NOTE: row offset not yet functional
# inputs: reps (number of sample replicates), offset_96 (column at which 96-well list begins)

def map_384_to_96(positions_384, reps, rows_384, cols_384, offset_96):

  full_cols = [sum([int(i[1:])==int(n) for i,_ in positions_384])for n in cols_384]
  col_groups = [[key, len(list(group))] for key, group in groupby(full_cols)]

  if reps == 4:
    multiplier = 2
  elif reps == 8:
    multiplier = 4
  elif reps == 2:
    multiplier = 1

  ic = 0
  positions_96 = []
  for g in col_groups:
    for i in range(offset_96+ic//multiplier,offset_96+ic+g[1]//multiplier):
      for j in range(0,g[0]//2):
        id = rows_384[j]+cols_384[i]
        positions_96.append(id)
    ic = ic+g[1]

  return positions_96

# Sample points from experiment design space
# DP and HEA concentration sampled using LHS
# Remaining monomer feed ratios sampled using Dirichlet sampling
# input n_samples: number of desired samples
# input bounds: upper and lower limits of variables chosen using LHS

def lhs_with_constraints(n_samples, bounds):

    # Sample DP and HEA using LHS
    # Constraint: HEA can only occur between 75% and 100% in steps of 5%
    lhs_inst = lhs(2, samples=n_samples)
    scaled_samples = np.floor(lhs_inst * (1+(bounds[:,1] - bounds[:,0]) / bounds[:,2]))
    backbone_vals = scaled_samples[:,0]*bounds[0,2]+75
    dp_vals = scaled_samples[:,1]*bounds[1,2]+bounds[1,0]

    samples = []

    # Apply constraints to monomer selection:
    # Constraint: All feed ratios must add to 100
    # Constraint: Up to 3 out of 4 possible remaining polymers can be used
    # Constraint: The sum of all monomers other than HEA must be less than the feed ratio of HEA

    for bb, dp in zip(backbone_vals, dp_vals):
        s = 100 - bb
        step = int(s/bounds[0,2])

        # Choose monomer to drop from the random selection
        zero_index = np.random.choice(4)

        # Split remaining non-HEA polymer content among remaining 3 copolymers
        raw = np.random.dirichlet(np.ones(3))
        dirichlet_inst = np.floor(raw*step).astype(int)
        diff = step - dirichlet_inst.sum()

        # Rounding procedure to convert random Dirichlet-selected floats into valid integers
        if diff > 0:
            residuals = raw*step - dirichlet_inst
            dirichlet_inst[np.argsort(-residuals)[:diff]] += 1
        elif diff < 0:
            residuals = raw*step - dirichlet_inst
            dirichlet_inst[np.argsort(residuals)[:abs(diff)]] -= 1

        # Combine dropped monomer index with rounded concentratinos for other comonomers
        full_dirichlet = np.insert(dirichlet_inst, zero_index, 0) * bounds[0,2]

        # Combine DP / HEA data with other co-monomer data
        sample = np.concatenate(([bb], full_dirichlet, [dp]))
        samples.append(sample)

    return np.array(samples)


# Format randomly selected polymers into Excel table to be used during preprocessing
def create_polymer_table(samples, names, exp_name):
    
    # Separate DP and feed ratios
    sample_mmers = pd.DataFrame(samples[:,:5],columns=names)
    sample_dps = pd.DataFrame(samples[:,5],columns=['Degree of Polymerization'])

    # Convert "one column per monomer" format into format with
    # 4 available monomers per polymer, with the name and percentage specified
    setup_rows = []
    max_mon = 4
    for _, row in sample_mmers.iterrows():
        mon_entries = list(row.items())
        mon_entries = [(mon, perc) for mon, perc in mon_entries if perc != 0 and not pd.isna(perc)]

        # Pad row with empty values if less than 4 monomers are used per polymer
        mon_entries += [('', '')] * (max_mon - len(mon_entries))
        mon_entries = mon_entries[:max_mon]

        row_dict = {}
        for i, (mon, perc) in enumerate(mon_entries, start=1):
            row_dict[f'Mon {i}'] = mon
            row_dict[f'% Mon {i}'] = perc

        setup_rows.append(row_dict)

    new_setup_data = pd.concat([sample_dps,pd.DataFrame(setup_rows)], axis=1)

    print(new_setup_data)

    # Add generic values to table
    # NOTE: Concentration (mM) of photoinitiator (ZnTPP), CTA and monomer are
    # set constant during automatic polymer table generation.
    # However, preproc.py can handle a diverse range of concentrations for
    # manually generated polymer libraries, including different concentrations
    # for each substance at each row.
    
    polymers = pd.DataFrame()
    polymers["Polymer ID"] = [f"{exp_name}_sample_{i+1}" for i in range(len(new_setup_data))]
    polymers["Experiment Code"] = None
    polymers["Degree of Polymerization"] = new_setup_data["Degree of Polymerization"]
    polymers["Volume"] = 200
    polymers["Solvent"] = "DMSO"
    polymers["Photoinitiator"] = "ZnTPP"
    polymers["Photoinitiator mM"] = 1
    polymers["CTA"] = 1
    polymers["CTA mM"] = 25

    for i in range(1, 5):
        mon_col = f"Mon {i}"
        pct_col = f"% Mon {i}"
        polymers[mon_col] = new_setup_data.get(mon_col)
        polymers[f"{mon_col} mM"] = np.where(polymers[mon_col] != '', 2000, np.nan  )
        polymers[f"% Mon {i}"] = new_setup_data.get(pct_col)

    return polymers

# Display polymer distribution as 2D PCA decomposition
def plot_polymer_pca(samples):
    pca = PCA(n_components=2)
    samples_pca = pca.fit_transform(samples)
    plt.figure(figsize=(8, 6))
    plt.scatter(samples_pca[:, 0], samples_pca[:, 1], alpha=0.7, edgecolors='k')
    plt.title("PCA of LHS Samples (projected to 2D)")
    plt.xlabel("Principal Component 1")
    plt.ylabel("Principal Component 2")
    plt.grid(True)
    plt.show()

def random_polymer_sampling(n_reagents, reagent_names, exp_name, plot_distribution=False):

    # bounds defined as: minimum, maximum, step
    dp_bounds = [150,800,25]
    mmer_bounds = [0,25,5]

    bounds = np.array([mmer_bounds, dp_bounds])

    # Sample a random selection of copolymers consisting of 
    # 5 possible monomers at different DPs, up to 4 monomers per polymer.
    lhs_samples = lhs_with_constraints(n_reagents, bounds)
    print(lhs_samples)

    if plot_distribution:
        plot_polymer_pca(lhs_samples)

    if len(reagent_names) != 5:
        print('Wrong number of reagent names provided!')
        return

    # Format selected polymers into Excel table
    polymers = create_polymer_table(lhs_samples, reagent_names, exp_name)
    return polymers