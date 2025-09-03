import time

from sampling import *
from preproc import *
from postproc import *

loop_actions = ['Active Learning with DLS']

#Function that creates the next file for a continuous multifile cycle to detect
def update_cycle_number(cycle, n_cycles, action, template, save_new_file):
    if cycle < n_cycles:

        # Selects appropriate file name
        base, ext = os.path.splitext(template)
        if cycle == 1:
            in_path = f'{base}{ext}'
        else:
            if save_new_file:
                in_path = f'{base}_{cycle}{ext}'
            else:
                in_path = f'{base}{ext}'

        # Locates and updates cycle number
        xp = pd.read_excel(in_path, header = None)
        action_label = xp.loc[xp[1]==action].index
        xp.iloc[action_label+1,1] +=1
        cycle = xp.iloc[action_label+1,1].values[0]

        start_chkp_label = xp.loc[xp[0]=='start checkpoint'].index
        xp.iloc[start_chkp_label,1] = 0
        start_chkp_label = xp.loc[xp[0]=='end checkpoint'].index
        xp.iloc[start_chkp_label,1] = 6

        # Saves experimental file for next iteration
        if save_new_file:
            out_path = f'{base}_{cycle}{ext}'
        else:
            out_path = f'{base}{ext}'
        print(f"Next file saved at: {out_path}")
        xp.to_excel(out_path, index=False, header=False)

#Function that runs a single generation of DBTL
def active_learning_with_dls(cycle, n_cycles, exp_name, data_subfolder, seed_size, next_gen_size, reagent_names, 
                            n_reps_per_384_plate, acq_select, acq_param, start_checkp=0, end_checkp=6):

    # Locate or create appropriate folders for file organization
    script_dir = os.getcwd()
    data_folder = os.path.join(script_dir, data_subfolder)
    if not os.path.exists(data_folder):
        os.makedirs(data_folder, exist_ok = True)
    exp_folder = os.path.join(data_folder, exp_name)
    if not os.path.exists(exp_folder):
        os.makedirs(exp_folder, exist_ok = True)

    # Configure the generation ID used to generate / locate experimental files
    if cycle == 1:
        exp_name_with_cycle = exp_name
    else:
        exp_name_with_cycle = f'{exp_name}_{cycle-1}'
    
    exp_name_next_cycle = f'{exp_name}_{cycle}'

    # Check variable types
    seed_size = int(seed_size)
    next_gen_size = int(next_gen_size)
    reagent_names = reagent_names.split(',')
    n_reps_per_384_plate = int(n_reps_per_384_plate)
    acq_param = float(acq_param)

    # Identify how many sub-experiment checkpoints are being ran
    # For a complete generation, 0 through 6 should be used
    controlled_steps = list(range(start_checkp, end_checkp))
    print(f'Checkpoints: {[x for x in controlled_steps]}')

    # Checkpoint 1: sampling and instrument command spreadsheet processing
    # Instruments used: none
    if 0 in controlled_steps:

        # Generate polymers randomly for seed generation
        if cycle==1:
            hamilton_path, lightbox_path, dls_prep_path, dls_run_path = preprocess_dls(
                exp_folder, exp_name_with_cycle, script_dir, data_folder, seed_size, reagent_names, n_reps_per_384_plate, sample_polymers_from_seed=True)
        
        # Select polymers from active learning results for subsequent generations
        else:
            hamilton_path, lightbox_path, dls_prep_path, dls_run_path = preprocess_dls(
                exp_folder, exp_name_with_cycle, script_dir, data_folder, seed_size, reagent_names, n_reps_per_384_plate, sample_polymers_from_seed=False)
    
        print(f'Hamilton steps saved to {hamilton_path}')
        print(f'Lightbox steps saved to {lightbox_path}')
        print(f'Hamilton steps saved to {dls_prep_path}')
        print(f'DLS steps saved to {dls_run_path}')

        time.sleep(5)

    # Checkpoint 2: preparation of polymer reagents for polymerization
    # Instruments used: Hamilton
    if 1 in controlled_steps:
        # Backup path attribute if spreadsheet preparation checkpoint isn't used
        if 0 not in controlled_steps:
            hamilton_path = os.path.join(exp_folder,f'hamilton_{exp_name_with_cycle}.xlsx')
        # Run spreadsheet
        hamilton_exp_seq = Sequence_Importer(hamilton_path.encode('unicode_escape').decode())
        Experimental_Runner(hamilton_exp_seq)
        reset_instruments()

    # Checkpoint 3: plate transfer to lightbox, polymerization, transfer back to Hamilton
    # Instruments used: xArm, lightbox
    if 2 in controlled_steps:
        # Backup path attribute if spreadsheet preparation checkpoint isn't used
        if 0 not in controlled_steps:
            lightbox_path = os.path.join(exp_folder,f'lightbox_{exp_name_with_cycle}.xlsx')
        # Run spreadsheet
        lightbox_exp_seq = Sequence_Importer(lightbox_path.encode('unicode_escape').decode())
        Experimental_Runner(lightbox_exp_seq)
        reset_instruments()

    # Checkpoint 4: dilution of polymers in preparation for DLS
    # Instruments used: Hamilton
    if 3 in controlled_steps:
        # Backup path attribute if spreadsheet preparation checkpoint isn't used
        if 0 not in controlled_steps:
            dls_prep_path = os.path.join(exp_folder,f'dls_prep_{exp_name_with_cycle}.xlsx')
        # Run spreadsheet
        dls_prep_exp_seq = Sequence_Importer(dls_prep_path.encode('unicode_escape').decode())
        Experimental_Runner(dls_prep_exp_seq)
        reset_instruments()

    # Checkpoint 5: DLS characterization
    # Instruments used: xArm, DLS
    if 4 in controlled_steps:
        # Backup path attribute if spreadsheet preparation checkpoint isn't used
        if 0 not in controlled_steps:
            dls_run_path = os.path.join(exp_folder,f'dls_run_{exp_name_with_cycle}.xlsx')
        # Run spreadsheet
        dls_run_exp_seq = Sequence_Importer(dls_run_path.encode('unicode_escape').decode())
        Experimental_Runner(dls_run_exp_seq)
        reset_instruments()

    # Checkpoint 6: Data processing, ML, active liearning
    # Instruments used: none
    if 5 in controlled_steps:
        cutoffs = [85, 0.0035, 0.025, 10]
        new_data = postprocess_dls(exp_folder, exp_name, cycle, n_reps_per_384_plate, next_gen_size, reagent_names, 
                                        acq_select, acq_param, cutoffs)

        # Save the next polymer sheet
        new_data_path = os.path.join(exp_folder, f'samples_{exp_name_next_cycle}.xlsx')
        new_data.to_excel(new_data_path, index=False, header=True)

        # Save the next loop checkpoint file
        exp_clone_template_path = os.path.join(script_dir, data_folder, f'{exp_name}.xlsx')
        update_cycle_number(cycle, n_cycles, 'Active Learning with DLS', exp_clone_template_path, save_new_file=False)
    