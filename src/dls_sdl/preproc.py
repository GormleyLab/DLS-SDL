from collections import Counter
from itertools import cycle, islice
import math
import numpy as np
import os
import pandas as pd

from sampling import *

# Edits path title to avoid overwriting existing file (adds a (1) or (n) if the file is found in the directory)
# input: file
# output: new file title

def check_dupe(path):
    base, ext = os.path.splitext(path)
    path_n = 1
    while os.path.exists(path):
        path = f"{base} ({path_n}){ext}"
        path_n += 1
    return path

# Insert step into existing template
# input xp: experiment DataFrame
# input steps: step DataFramses
# input positions: numerical positions of where to add steps
# output: new experiment DataFrame

def add_step(xp, steps, positions):
  if isinstance(steps, pd.DataFrame):
    steps = [steps]
  if isinstance(positions, int):
    positions = [positions]

  for i, step in enumerate(steps):
    erase_step_label_row = step.loc[(step[0]=='Step')].index
    step.iloc[erase_step_label_row,1] = 0

    split_label = xp.loc[(xp[0]=='Step') & (xp[1]==int(positions[i]))].index[0] - 1
    xp = pd.concat([xp.loc[:split_label], step, xp.loc[split_label:]]).reset_index(drop=True)

  return xp

# Change all step numbers to ascending order starting from 1
def update_step_numbers(xp):
  for i, row in enumerate(xp.loc[xp[0]=='Step'].index):
    xp.iloc[row,1]=i+1
  return xp

# Add necessary variables to liquid transfer SmartStep
def fill_transfer_step(xp, step, volumes, n_samples, aliquots, positions, pickup=0, dropoff=0):

    step_label = xp.loc[(xp[0]=='Step') & (xp[1]==int(step))].index
    if xp.iloc[step_label+1,1].values == 'Transfer':

        xp.iloc[step_label+2,1] = ",".join(map(str, volumes)).strip()
        xp.iloc[step_label+3,1] = n_samples
        xp.iloc[step_label+4,1] = ",".join(map(str, aliquots)).strip()
        xp.iloc[step_label+5,1] = ",".join(map(str, positions)).strip()
        if pickup: xp.iloc[step_label+6,1] = pickup
        if dropoff: xp.iloc[step_label+7,1] = dropoff

    return xp

# Change order of numerical positions into as many possible
# concatenated ascending lists as possible
# e.g. sort [1 1 1 1 2 2 3 3 3] into [1 2 3 1 2 3 1 3 1]
# Used to generate more efficient aspiration pattern
# input seq: list of numbers to be sorted
# input split: sort two halves of the list separately, with split as the separation value

def cyclic_sort_with_split(seq, split=0):

  counts = Counter(seq)
  sorted_seq = sorted(counts.keys())
  range_1 = [n for n in sorted_seq if n < split]
  range_2 = [n for n in sorted_seq if n >= split]

  out_seq = []

  def cycle_range(range_list):

    nonlocal out_seq
    while any(counts[n] > 0 for n in range_list):
      for n in range_list:
        if counts[n] > 0:
          out_seq.append(n)
          counts[n] -= 1

  cycle_range(range_1)
  cycle_range(range_2)

  return out_seq

# For several data_lists of identical length, perform cyclic sort on one list, seq,
# and adjust the positions of items according to the new positions of the
# cyclically sorted list
# input seq: list of numbers to be sorted, *data_lists = any number of additional lists
# input split: sort two halves of the list separately, with split as the separation value

def sort_with_partners(seq,split=0, *data_lists):

    def seq_to_dict(seq, data):
        seq_to_data = {}
        for n, d in zip(seq, data):
            if n not in seq_to_data:
                seq_to_data[n] = []
            seq_to_data[n].append(d)
        return seq_to_data

    seq_to_data_dicts = [seq_to_dict(seq, data) for data in data_lists]
    sorted_seq = cyclic_sort_with_split(seq, split)

    out_data_lists = [[] for _ in range(len(data_lists))]
    for num in sorted_seq:
        for i, seq_to_data in enumerate(seq_to_data_dicts):
            out_data_lists[i].append(seq_to_data[num].pop(0))

    return (sorted_seq, *out_data_lists)

# Accept spreadsheet of polymer samples and organize it into
# 4 additional spreadsheets containing commands to be executed
# by instruments involved in the workflow

def preprocess_dls(exp_folder, exp_name, script_dir, data_folder, n_reagents, reagent_names, n_384_reps = 4, 
                   sample_polymers_from_seed=False, create_hamilton_file=True, create_lightbox_file=True,
                   create_dls_prep_file=True, create_dls_run_file=True):

    template_subfolder = 'PolyCraft Templates'
    template_folder = os.path.join(script_dir, template_subfolder)

    dmso_step_id = 3676

    # Constants selected for PET-RAFT
    final_monomer_conc = 1000
    cta_to_photoinitiator = 50

    # Artificial limit to prevent Hamilton deck layout from defaulting mid-run
    tip_limit = 96*4

    # Maximum volume to be eluted from 1.5 mL Eppendorf tube containing ~1 mL of liquid
    volume_threshold = 800 #microliters

    # Number of glass tubes containing solvent
    n_solvent_tubes = 4

    print_results = True

    # Use cyclic sort for small tube positions or not
    efficient_aspiration = True

    # Aspirate solvent using fresh tips each time or not (usually False)
    replace_solvent_tips = False

    if sample_polymers_from_seed:
        # LHS + Dirichlet sampling
        setup_data = random_polymer_sampling(n_reagents, reagent_names, exp_name)
        save_path = os.path.join(exp_folder,f'samples_{exp_name}.xlsx')
        save_path = check_dupe(save_path)
        setup_data.to_excel(save_path, index=False)
    else:
        # Sampling from previous generation of active learning
        setup_file = os.path.join(exp_folder,f'samples_{exp_name}.xlsx')
        setup_data = pd.read_excel(setup_file, index_col=None)

    n_setup_samples = len(setup_data)

    # Categorize each unique combination of a monomer and its concentration as a different monomer
    monomers_w_conc_set = [setup_data.groupby([f'Mon {i}', f'Mon {i} mM'], sort=False).size().reset_index()
                    .rename(columns={f'Mon {i}': 'Reag', f'Mon {i} mM': 'Reag mM', 0:'size'}) for i in range(1, 5)]
    monomers_w_conc = pd.concat(monomers_w_conc_set, ignore_index=True).drop_duplicates(subset=['Reag', 'Reag mM']).reset_index(drop=True)

    print(monomers_w_conc)

    # Isolate DP column
    dp_list = setup_data['Degree of Polymerization']

    # Categorize each unique combination of a CTA type and its concentration as a unique CTA
    ctas = setup_data.groupby(['CTA','CTA mM'], sort=False).size().reset_index().rename(columns={'CTA': 'Reag', 'CTA mM': 'Reag mM', 0:'size'})
    ctas['Reag'] = [str(x).replace(str(x), 'CTA '+str(x)) for x in ctas['Reag']]

    # Categorize each unique combination of a photoinitiator (ZnTPP) and its concentration as a photoinitiator
    photoinits = setup_data.groupby(['Photoinitiator','Photoinitiator mM'], sort=False).size().reset_index()
    photoinits = photoinits.rename(columns={'Photoinitiator': 'Reag', 'Photoinitiator mM': 'Reag mM', 0:'size'})

    # Combine all reagents together into single set
    reagents = pd.concat([monomers_w_conc,ctas,photoinits], axis=0, ignore_index=True)
    reagents = reagents.replace(np.nan,'')
    reagent_titles = (reagents['Reag']+' '+reagents['Reag mM'].astype(float).astype(str)).str.strip()

    # Set up solvent names in the same format as other reagents
    solvents = setup_data.groupby(['Solvent'], sort=False).size().reset_index().rename(columns={'Solvent': 'Reag', 0:'size'})
    solvent_titles = (solvents['Reag'].astype(str)).str.strip()

    # Define 96-well positions for polymer reagents in terms of how many samples are being polymerized
    offset_96 = 0
    positions_96 = [rows_96[x%8]+cols_96[offset_96+x//8] for x in range(n_setup_samples)]

    # Set up empty dictionaries to be used as spreadsheet inputs later on

    aliquots = {col: [] for col in reagent_titles}
    used_reagent = {col: [] for col in reagent_titles}
    feed_ratios = pd.DataFrame(columns=[row['Reag']+' '+str(row['Reag mM']) for _,row in monomers_w_conc.iterrows()])

    reagent_ids = {col: [] for col in reagent_titles}
    reagent_300_ids = {col: [] for col in reagent_titles}
    reagent_50_ids = {col: [] for col in reagent_titles}
    solvent_300_ids = {col: [] for col in solvent_titles}
    solvent_50_ids = {col: [] for col in solvent_titles}

    reagent_volumes = {col: [] for col in reagent_titles}
    reagent_300_volumes = {col: [] for col in reagent_titles}
    reagent_50_volumes = {col: [] for col in reagent_titles}
    solvent_300_volumes = {col: [] for col in solvent_titles}
    solvent_50_volumes = {col: [] for col in solvent_titles}

    reagent_positions = {col: [] for col in reagent_titles}
    reagent_300_positions = {col: [] for col in reagent_titles}
    reagent_50_positions = {col: [] for col in reagent_titles}
    solvent_300_positions = {col: [] for col in solvent_titles}
    solvent_50_positions = {col: [] for col in solvent_titles}

    reagent_aliquots = {col: [] for col in reagent_titles}
    reagent_300_aliquots = {col: [] for col in reagent_titles}
    reagent_50_aliquots = {col: [] for col in reagent_titles}
    solvent_300_aliquots = {col: [] for col in solvent_titles}
    solvent_50_aliquots = {col: [] for col in solvent_titles}

    # Generic method
    # Add reagent name (ID) to reagent list
    # Add reagent volume to volume list
    # Add reagent 96-well position to 96-well position list

    def fill_reagent_dicts(ids, volumes, positions, id, volume, position):
        ids[id] += [id]
        volumes[id] += [volume]
        positions[id] += [position]
        return ids, volumes, positions
    
    # Iterate through each polymer to be synthesized and create separate aspiration data lists for monomers, CTA, photoinitiators and solvents

    for i, row in setup_data.iterrows():

        # Set DP and total polymer volume as constants for volume calculation
        dp = row['Degree of Polymerization']
        sample_volume = row['Volume']
        total_mon_volume = 0
        feed_ratio = {col: 0 for col in [row['Reag']+' '+str(row['Reag mM']) for _,row in monomers_w_conc.iterrows()]}

        for j in range(1, 5):

            mon = row[f'Mon {j}']
            mon_conc = float(row[f'Mon {j} mM'])
            mon_perc = row[f'% Mon {j}']

            # Check if monomer and concentration combo actually exists
            if sum(reagents.drop(columns='size').isin([mon, mon_conc]).all(axis=1)):

                # Calculate monomer volume depending on feed ratio (mon_perc), concentration and sample volume
                mon_volume = final_monomer_conc/mon_conc*sample_volume * mon_perc/100

                # Add monomer information to generic reagent list, with 96-well position corresponding to the order of the polymer entry in the list / loop
                fill_reagent_dicts(reagent_ids, reagent_volumes, reagent_positions, mon+' '+str(mon_conc), mon_volume, positions_96[i])
                total_mon_volume += mon_volume

                # Sort into either 50 uL tip list or 300 uL tip list
                if mon_volume > 50:
                    fill_reagent_dicts(reagent_300_ids, reagent_300_volumes, reagent_300_positions, mon+' '+str(mon_conc), mon_volume, positions_96[i])
                elif mon_volume <=300:
                    fill_reagent_dicts(reagent_50_ids, reagent_50_volumes, reagent_50_positions, mon+' '+str(mon_conc), mon_volume, positions_96[i])
                else:
                    print(f'No tips available to aliquot {mon_volume} uL of {mon} ({str(mon_conc)})!!!')

                feed_ratio[mon+' '+str(mon_conc)] = mon_perc

        feed_ratios = pd.concat([feed_ratios, pd.DataFrame.from_dict(feed_ratio, orient='index').T], ignore_index=True)

        cta = 'CTA '+str(row['CTA'])
        cta_conc = float(row['CTA mM'])

        # Check if CTA and concentration combo actually exists
        if sum(reagents.drop(columns='size').isin([cta, cta_conc]).all(axis=1)):

            # Calculate CTA volume depending on concentration, DP and total voume of all monomers
            cta_volume = final_monomer_conc/dp * sample_volume/cta_conc

            # Add CTA information to generic reagent list, with 96-well position corresponding to the order of the polymer entry in the list / loop
            fill_reagent_dicts(reagent_ids, reagent_volumes, reagent_positions, cta+' '+str(cta_conc), cta_volume, positions_96[i])

            # Sort into either 50 uL tip list or 300 uL tip list
            if cta_volume > 50:
                fill_reagent_dicts(reagent_300_ids, reagent_300_volumes, reagent_300_positions, cta+' '+str(cta_conc), cta_volume, positions_96[i])
            elif cta_volume <=300:
                fill_reagent_dicts(reagent_50_ids, reagent_50_volumes, reagent_50_positions, cta+' '+str(cta_conc), cta_volume, positions_96[i])

        phin = row['Photoinitiator']
        phin_conc = float(row['Photoinitiator mM'])

        # Check if photoinitiator and concentration combo actually exists
        if sum(reagents.drop(columns='size').isin([phin, phin_conc]).all(axis=1)):

            # Calculate photoinitiator volume as a measure of CTA concentration and volume
            phin_volume = cta_conc / phin_conc / cta_to_photoinitiator * cta_volume
            
            # Add photoinitiator information to generic reagent list, with 96-well position corresponding to the order of the polymer entry in the list / loop
            fill_reagent_dicts(reagent_ids, reagent_volumes, reagent_positions, phin+' '+str(phin_conc), phin_volume, positions_96[i])

            # Sort into either 50 uL tip list or 300 uL tip list
            if phin_volume > 50:
                fill_reagent_dicts(reagent_300_ids, reagent_300_volumes, reagent_300_positions, phin+' '+str(phin_conc), phin_volume, positions_96[i])
            elif phin_volume <= 300:
                fill_reagent_dicts(reagent_50_ids, reagent_50_volumes, reagent_50_positions, phin+' '+str(phin_conc), phin_volume, positions_96[i])

        solvent = row['Solvent']

        # Check if solvent actually exists
        if solvent in solvents['Reag'].values.tolist():

            # Calculate solvent volume as remainder of total sample volume not yet filled by other reagents
            solvent_volume = float(np.round(sample_volume-total_mon_volume-cta_volume-phin_volume, 4))

            # Sort into either 50 uL tip list or 300 uL tip list
            if solvent_volume > 50:
                fill_reagent_dicts(solvent_300_ids, solvent_300_volumes, solvent_300_positions, solvent, solvent_volume, positions_96[i])
            elif solvent_volume <= 300:
                fill_reagent_dicts(solvent_50_ids, solvent_50_volumes, solvent_50_positions, solvent, solvent_volume, positions_96[i])

    # Assign Hamilton robot positions to Eppendorf tubes containing polymer reagents
    global_aliquot = 0

    for r, vols in reagent_volumes.items():

        aliquots[r] = list(range(global_aliquot+1,global_aliquot+1+math.ceil(sum(vols)/volume_threshold)))
        used_reagent[r] = sum(vols)

        # Determine total volume of each reagent that is needed, assign aliquot as the loop iterates
        total_reagent_volume = 0
        for v in vols:
            total_reagent_volume += v

            # Split large reagent quantities into multiple aliquots to volume threshold
            local_aliquot = math.ceil(global_aliquot + total_reagent_volume / volume_threshold)
            reagent_aliquots[r] += [local_aliquot]

            # Sort into either 50 uL tip list or 300 uL tip list
            if v > 50:
                reagent_300_aliquots[r] += [local_aliquot]
            else:
                reagent_50_aliquots[r] += [local_aliquot]

        global_aliquot = local_aliquot

    # Convert dictionaries into lists

    reagent_300_ids = [id for reag in reagent_300_ids.values() for id in reag]
    reagent_300_volumes = [vol for reag in reagent_300_volumes.values() for vol in reag]
    reagent_300_positions = [pos for reag in reagent_300_positions.values() for pos in reag]
    reagent_300_aliquots = [alq for reag in reagent_300_aliquots.values() for alq in reag]

    reagent_50_ids = [id for reag in reagent_50_ids.values() for id in reag]
    reagent_50_volumes = [vol for reag in reagent_50_volumes.values() for vol in reag]
    reagent_50_positions = [pos for reag in reagent_50_positions.values() for pos in reag]
    reagent_50_aliquots = [alq for reag in reagent_50_aliquots.values() for alq in reag]

    # Account for cyclic sort when adjusting reagent information
    if efficient_aspiration:

        first_cta_id = aliquots[ctas['Reag'][0]+' '+str(float(ctas['Reag mM'][0]))][0]
        reagent_50_aliquots, reagent_50_ids, reagent_50_volumes, reagent_50_positions = sort_with_partners(
            reagent_50_aliquots, first_cta_id, reagent_50_ids, reagent_50_volumes, reagent_50_positions)

        reagent_300_aliquots, reagent_300_ids, reagent_300_volumes, reagent_300_positions = sort_with_partners(
            reagent_300_aliquots, first_cta_id, reagent_300_ids, reagent_300_volumes, reagent_300_positions)

    solvent_300_ids = [id for solv in solvent_300_ids.values() for id in solv]
    solvent_300_volumes = [vol for solv in solvent_300_volumes.values() for vol in solv]
    solvent_300_positions = [pos for solv in solvent_300_positions.values() for pos in solv]
    solvent_300_aliquots = list(islice(cycle(range(1,1+n_solvent_tubes)), len(solvent_300_volumes)))

    solvent_50_ids = [id for solv in solvent_50_ids.values() for id in solv]
    solvent_50_volumes = [vol for solv in solvent_50_volumes.values() for vol in solv]
    solvent_50_positions = [pos for solv in solvent_50_positions.values() for pos in solv]
    solvent_50_aliquots = list(islice(cycle(range(1,1+n_solvent_tubes)), len(solvent_50_volumes)))

    # Calculate the number of each tip type required for  polymer preparation

    if replace_solvent_tips:
        n_solvent_50_tips = len(solvent_50_volumes)
        n_solvent_300_tips = len(solvent_300_volumes)
    else:
        n_solvent_50_tips = min(len(solvent_50_volumes),8)
        n_solvent_300_tips = min(len(solvent_300_volumes),8)

    n_reagent_50_tips = len(reagent_50_volumes)
    n_reagent_300_tips = len(reagent_300_volumes)
    total_50_tips = n_reagent_50_tips+n_solvent_50_tips
    total_300_tips = n_reagent_300_tips+n_solvent_300_tips+n_setup_samples

    if total_300_tips > tip_limit:
        print(f'WARNING: {total_300_tips - tip_limit} more 300 uL tips required')

    if total_50_tips > tip_limit:
        print(f'WARNING: {total_50_tips - tip_limit} more 50 uL tips required')

    # Console output

    if print_results:

        print('Solvent IDs (50 uL)')
        print(solvent_50_ids)
        print('Solvent Volumes (50 uL)')
        print(solvent_50_volumes)
        print('Solvent Positions (50 uL)')
        print(solvent_50_positions)
        print('Solvent Aliquots (50 uL)')
        print(solvent_50_aliquots)
        print('--------------------------------------------\n')

        print('Solvent IDs (300 uL)')
        print(solvent_300_ids)
        print('Solvent Volumes (300 uL)')
        print(solvent_300_volumes)
        print('Solvent Positions (300 uL)')
        print(solvent_300_positions)
        print('Solvent Aliquots (300 uL)')
        print(solvent_300_aliquots)
        print('--------------------------------------------\n')

        print('Reagent IDs (50 uL)')
        print(reagent_50_ids)
        print('Reagent Volumes (50 uL)')
        print(reagent_50_volumes)
        print('Reagent Volumes (50 uL)')
        print(reagent_50_positions)
        print('Reagent Aliquots (50 uL)')
        print(reagent_50_aliquots)
        print('--------------------------------------------\n')

        print('Reagent IDs (300 uL)')
        print(reagent_300_ids)
        print('Reagent Volumes (300 uL)')
        print(reagent_300_volumes)
        print('Reagent Positions (300 uL)')
        print(reagent_300_positions)
        print('Reagent Aliquots (300 uL)')
        print(reagent_300_aliquots)
        print('--------------------------------------------\n')

        print(f'Total 50 uL Volumes: {total_50_tips}')
        print(f'Total 300 uL Volumes: {total_300_tips}')
        print(f'Total solvent used: {sum(solvent_300_volumes)+sum(solvent_50_volumes)}')

        print('\nAliquot Positions')
        for key, values in aliquots.items():
            print(key, values)

        print('\nTotal Used Reagent Volume')
        for key, values in used_reagent.items():
            print(key, values)
    
    # Create sample positions spreadsheet summarizing Hamilton aliquot positions to the user
    summary_data = {
        'Metric': ['Total 50 uL Volumes', 'Total 300 uL Volumes', 'Total solvent used'],
        'Value': [
            total_50_tips,
            total_300_tips,
            sum(solvent_300_volumes) + sum(solvent_50_volumes)
        ]
    }

    df_summary = pd.DataFrame(summary_data)

    df_aliquots = pd.DataFrame.from_dict(aliquots, orient='index')
    df_aliquots.index.name = 'Aliquot'
    df_aliquots.reset_index(inplace=True)

    df_reagents = pd.DataFrame.from_dict(used_reagent, orient='index')
    df_reagents.index.name = 'Reagent'
    df_reagents.reset_index(inplace=True)

    sample_pos_path = os.path.join(exp_folder,f'sample_positions_{exp_name}.xlsx')
    sample_pos_path = check_dupe(sample_pos_path)
    with pd.ExcelWriter(sample_pos_path, engine='openpyxl') as writer:
        df_summary.to_excel(writer, sheet_name='Summary', index=False)
        df_aliquots.to_excel(writer, sheet_name='Aliquot Positions', index=False, header=False)
        df_reagents.to_excel(writer, sheet_name='Used Reagents', index=False, header=False)

    # Prepare spreadsheet for Checkpoint 2: polymer preparation using Hamilton
    def hamilton_steps(exp_name, save=False):
        xp = pd.read_excel(os.path.join(template_folder,'Polymer Synthesis Template.xlsx'), header=None)

        fill_transfer_step(xp, 2, solvent_300_volumes, len(solvent_300_volumes), solvent_300_aliquots, solvent_300_positions)
        fill_transfer_step(xp, 3, solvent_50_volumes, len(solvent_50_volumes), solvent_50_aliquots, solvent_50_positions)
        fill_transfer_step(xp, 4, reagent_300_volumes, len(reagent_300_volumes), reagent_300_aliquots, reagent_300_positions)
        fill_transfer_step(xp, 5, reagent_50_volumes, len(reagent_50_volumes), reagent_50_aliquots, reagent_50_positions)
        fill_transfer_step(xp, 6, [150.0]*n_setup_samples, n_setup_samples, positions_96, positions_96)

        if len(solvent_300_volumes)==0:
            step_label = xp.loc[(xp[0]=='Step') & (xp[1]==2)].index[0]-1
            next_step_label = xp.loc[(xp[0]=='Step') & (xp[1]==3)].index[0]
            xp = pd.concat([xp.loc[:step_label],xp.loc[next_step_label:]]).reset_index(drop=True)
        
        if len(solvent_50_volumes)==0:
            step_label = xp.loc[(xp[0]=='Step') & (xp[1]==3)].index[0]-1
            next_step_label = xp.loc[(xp[0]=='Step') & (xp[1]==4)].index[0]
            xp = pd.concat([xp.loc[:step_label],xp.loc[next_step_label:]]).reset_index(drop=True)

        save_path = os.path.join(exp_folder,f'hamilton_{exp_name}.xlsx')
        if save:
            save_path = check_dupe(save_path)
            update_step_numbers(xp).to_excel(save_path, index=False, header=False)
        return save_path

    # Prepare spreadsheet for Checkpoint 3: lightbox polymerization
    def lightbox_steps(exp_name, hours, save=False):
        xp = pd.read_excel(os.path.join(template_folder,'Lightbox Template.xlsx'), header=None)

        light_on_step = xp.loc[(xp[0]=='Step') & (xp[1]==4)].index
        xp.iloc[light_on_step+2,1] = 3600*hours
        xp.iloc[light_on_step+4,1] = ",".join(map(str, positions_96)).strip()

        save_path = os.path.join(exp_folder,f'lightbox_{exp_name}.xlsx')
        if save:
            save_path = check_dupe(save_path)
            update_step_numbers(xp).to_excel(save_path, index=False, header=False)
        return save_path

    # Prepare spreadsheet for Checkpoint 4: plate preparation for DLS
    def dls_prep_steps(exp_name, n_384_reps, save=False, add_spare_solvent = False):
        xp = pd.read_excel(os.path.join(template_folder,'DLS Preparation Template.xlsx'), header=None)

        dls_sample_96_volume = 250
        dls_sample_384_volume = 70
        dilution_constant = 0.04
        dls_polymer_volume = dls_sample_96_volume * dilution_constant * 2
        dls_solvent_volume = dls_sample_96_volume - dls_polymer_volume
        dls_sample_96_halfvol = dls_sample_96_volume/2

        solvent_pickup = 'SmallTubes1'
        spare_solvent_position = 32
        n_dls_96_samples = n_setup_samples
        full_positions_96 = positions_96

        if add_spare_solvent==True:

            spare_step = pd.read_excel(os.path.join(template_folder, 'DMSO Step.xlsx'), header=None)
            spare_step = spare_step.dropna(how='all')
            spare_step.loc[len(spare_step)] = np.nan
            id = dmso_step_id

            if spare_step.iloc[spare_step.loc[(spare_step[0]=='Step')].index,1].values[0] == id:
                n_spare_samples = 8-n_setup_samples%8
                n_dls_96_samples = n_setup_samples + n_spare_samples
                spare_volumes = [dls_polymer_volume]*n_spare_samples
                spare_aliquots = [spare_solvent_position]*n_spare_samples
                spare_positions = [rows_96[x%8]+cols_96[offset_96+x//8] for x in range(n_setup_samples, n_dls_96_samples)]
                full_positions_96 = positions_96+spare_positions
                spare_step = fill_transfer_step(spare_step, id, spare_volumes, n_spare_samples, spare_aliquots, spare_positions, pickup=solvent_pickup)

                xp = add_step(xp, spare_step, 5)

        n_solvent_tubes = 4

        offset_96 = 0
        pos_96_1 = [rows_96[x%8]+cols_96[offset_96+x//8] for x in range(n_setup_samples)]
        full_pos_96_1 = [rows_96[x%8]+cols_96[offset_96+x//8] for x in range(n_dls_96_samples)]
        offset_96 = 0
        pos_96_2 = [rows_96[x%8]+cols_96[offset_96+x//8] for x in range(n_setup_samples)]
        full_pos_96_2 = [rows_96[x%8]+cols_96[offset_96+x//8] for x in range(n_dls_96_samples)]

        offset_384 = 0
        dls_384_positions = map_96_to_384(full_positions_96, n_384_reps, 'per_step', rows_96, cols_96, offset_384)

        if n_384_reps == 4:
            solvent_aliquots = list(islice(cycle(range(1,1+n_solvent_tubes)), n_dls_96_samples))
            fill_transfer_step(xp, 2, [100]*n_setup_samples, n_setup_samples, positions_96, positions_96)
            fill_transfer_step(xp, 3, [dls_solvent_volume]*n_dls_96_samples, n_dls_96_samples, solvent_aliquots, full_pos_96_1)
            fill_transfer_step(xp, 4, [dls_sample_96_halfvol]*n_dls_96_samples, n_dls_96_samples, solvent_aliquots, full_pos_96_1)
            fill_transfer_step(xp, 5, [dls_polymer_volume]*n_setup_samples, n_setup_samples, positions_96, full_pos_96_1)
            fill_transfer_step(xp, 6, [dls_sample_96_halfvol]*n_dls_96_samples, n_dls_96_samples, full_pos_96_1, full_pos_96_1)
            fill_transfer_step(xp, 7, [dls_sample_96_halfvol]*n_dls_96_samples, n_dls_96_samples, full_pos_96_1, full_pos_96_1)
            fill_transfer_step(xp, 8, [dls_sample_96_halfvol]*n_dls_96_samples, n_dls_96_samples, solvent_aliquots, full_pos_96_1)
            fill_transfer_step(xp, 9, [dls_sample_384_volume]*n_dls_96_samples, n_dls_96_samples, full_pos_96_1, dls_384_positions[0])
            fill_transfer_step(xp, 10, [dls_sample_384_volume]*n_dls_96_samples, n_dls_96_samples, full_pos_96_1, dls_384_positions[2])
            fill_transfer_step(xp, 11, [dls_sample_384_volume]*n_dls_96_samples, n_dls_96_samples, full_pos_96_1, dls_384_positions[1])
            fill_transfer_step(xp, 12, [dls_sample_384_volume]*n_dls_96_samples, n_dls_96_samples, full_pos_96_1, dls_384_positions[3])

        if n_384_reps == 8:
            solvent_aliquots = list(islice(cycle(range(1,1+n_solvent_tubes)), n_dls_96_samples*2))
            fill_transfer_step(xp, 2, [100]*n_setup_samples, n_setup_samples, positions_96, positions_96)
            fill_transfer_step(xp, 3, [dls_solvent_volume]*n_dls_96_samples*2, n_dls_96_samples*2, solvent_aliquots, full_pos_96_1+full_pos_96_2)
            fill_transfer_step(xp, 4, [dls_sample_96_halfvol]*n_dls_96_samples*2, n_dls_96_samples*2, solvent_aliquots, full_pos_96_1+full_pos_96_2)
            fill_transfer_step(xp, 5, [dls_polymer_volume]*n_setup_samples*2, n_setup_samples*2, positions_96+positions_96, pos_96_1+pos_96_2)
            fill_transfer_step(xp, 6, [dls_sample_96_halfvol]*n_dls_96_samples*2, n_dls_96_samples*2, full_pos_96_1+full_pos_96_2, full_pos_96_1+full_pos_96_2)
            fill_transfer_step(xp, 7, [dls_sample_96_halfvol]*n_dls_96_samples*2, n_dls_96_samples*2, full_pos_96_1+full_pos_96_2, full_pos_96_1+full_pos_96_2)
            fill_transfer_step(xp, 8, [dls_sample_96_halfvol]*n_dls_96_samples*2, n_dls_96_samples*2, solvent_aliquots, full_pos_96_1+full_pos_96_2)
            fill_transfer_step(xp, 9, [dls_sample_384_volume]*n_dls_96_samples*2, n_dls_96_samples*2, full_pos_96_1+full_pos_96_1, dls_384_positions[0]+dls_384_positions[1])
            fill_transfer_step(xp, 10, [dls_sample_384_volume]*n_dls_96_samples*2, n_dls_96_samples*2, full_pos_96_2+full_pos_96_2, dls_384_positions[4]+dls_384_positions[5])
            fill_transfer_step(xp, 11, [dls_sample_384_volume]*n_dls_96_samples*2, n_dls_96_samples*2, full_pos_96_1+full_pos_96_1, dls_384_positions[2]+dls_384_positions[3])
            fill_transfer_step(xp, 12, [dls_sample_384_volume]*n_dls_96_samples*2, n_dls_96_samples*2, full_pos_96_2+full_pos_96_2, dls_384_positions[6]+dls_384_positions[7])

        save_path = os.path.join(exp_folder,f'dls_prep_{exp_name}.xlsx')
        if save:
            save_path = check_dupe(save_path)
            update_step_numbers(xp).to_excel(save_path, index=False, header=False)

        dls_wells = map_96_to_384(full_positions_96, n_384_reps, 'alphabetical', rows_96, cols_96, offset_384)
        return save_path, dls_wells
    
    # Prepare spreadsheet for Checkpoint 5: DLS characterization
    def dls_run_steps(exp_name, dls_wells, save=False):
        xp = pd.read_excel(os.path.join(template_folder,'DLS Run Template.xlsx'), header=None)

        action_label = xp.loc[xp[1]=='Run DLS'].index

        xp.iloc[action_label+1,1] = script_dir
        xp.iloc[action_label+2,1] = data_folder
        xp.iloc[action_label+3,1] = exp_name
        xp.iloc[action_label+4,1] = ",".join(map(str, dls_wells)).strip()
        xp.iloc[action_label+5,1] = 384
        xp.iloc[action_label+6,1] = 5
        xp.iloc[action_label+7,1] = 8
        xp.iloc[action_label+8,1] = 25
        xp.iloc[action_label+9,1] = 20
        xp.iloc[action_label+10,1] = True
        xp.iloc[action_label+11,1] = 'D'
        xp.iloc[action_label+12,1] = 'M'
        xp.iloc[action_label+13,1] = 'R'
        xp.iloc[action_label+14,1] = False
        xp.iloc[action_label+15,1] = False
        
        save_path = os.path.join(exp_folder,f'dls_run_{exp_name}.xlsx')
        if save:
            save_path = check_dupe(save_path)
            update_step_numbers(xp).to_excel(save_path, index=False, header=False)
        return save_path

    # Prepare all spreadsheets and return paths
    if create_hamilton_file:
        hamilton_path = hamilton_steps(exp_name, save=True)
    if create_lightbox_file:
        lightbox_path = lightbox_steps(exp_name, hours=3, save=True)
    if create_dls_prep_file:
        dls_prep_path, dls_wells = dls_prep_steps(exp_name, n_384_reps, save=True)
    if create_dls_run_file:
        dls_run_path = dls_run_steps(exp_name, dls_wells, save=True)
    
    return hamilton_path, lightbox_path, dls_prep_path, dls_run_path