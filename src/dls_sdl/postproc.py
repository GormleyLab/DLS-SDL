from itertools import product
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import stats
import shap

from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF,ConstantKernel,DotProduct,WhiteKernel
from sklearn.metrics import mean_squared_error,mean_absolute_error,r2_score
# from sklearn.metrics import root_mean_squared_error
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_score, train_test_split, RandomizedSearchCV, GridSearchCV
from sklearn.preprocessing import StandardScaler

from sampling import rows_96, cols_96, map_96_to_384, create_polymer_table

import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings('ignore', category=ConvergenceWarning)

# Convert "4 monomer column" format into format where
# each monomer feed ratio is stored in a separate column
def get_monomers(data):
    
    monomer_names = [x for x in pd.unique(data[[f'Mon {i}' for i in range(1, 5)]].values.ravel()).tolist() if str(x)!='nan']
    monomer_table = pd.DataFrame(columns=monomer_names)

    for _, row in data.iterrows():
        temp = {col: 0 for col in monomer_names}
        for i in range(1, 5):
            mon = row[f'Mon {i}']
            perc = row[f'% Mon {i}']
            if mon in monomer_names:
                temp[mon] = perc
        monomer_table = pd.concat([monomer_table, pd.DataFrame.from_dict(temp, orient='index').T], ignore_index=True)

    return monomer_table

# Generate all possible polymer formulations for the design space that the polymers
# Represents all possible outputs of lhs_with_constraints

def get_design_space():

    options = [float(x) for x in list(np.arange(0,30,5))+list(np.arange(75,105,5))]
    x_range = []
    for c in product(options, repeat = 5):
        if sum(c) == 100 and all(x<=25 for x in c[1:]) and c[0]>=75 and sum(1 for x in c[1:] if x != 0) <= 3:
            x_range.append(c)
    
    x_range = [list(x) for x in x_range]

    full_x_range = np.array([])

    dp_options = [float(x) for x in list(np.arange(150,825,25))]

    for opt in dp_options:
        new_x_range = np.hstack((x_range,np.full((len(x_range), 1), opt)))
        if full_x_range.size == 0:
            full_x_range = new_x_range
        else:
            full_x_range = np.vstack((full_x_range,new_x_range))

    return full_x_range

# Cut dataset remainder to allow for identical number of datapoints per stratified K-fold split
# Samples with highest variance removed during training, retained as part of dataset for future active learning generations (but might still get removed)
# Feed ratio and DP parameters used as features, Rh used as labels, variance used for noise regularization 

def prep_data_for_gpr(dataset, n_cv, rd_state):
    
    cut_data = dataset.sample(frac=1, random_state=rd_state)
    print('Removed samples due to high standard deviation:')
    print(cut_data.nlargest(len(cut_data)%n_cv, 'alpha Rh'))
    cut_data.drop(cut_data.nlargest(len(cut_data)%n_cv, 'alpha Rh').index,inplace=True)
    
    X = cut_data.drop(['Rh', 'alpha Rh'], axis=1)
    y = cut_data['Rh']
    print(f'{len(y)} samples remaining in filtered data:')
    print(X)
    z = cut_data['alpha Rh']
    
    return X, y, z

# Process data, train model and obtain new polymer samples using machine learnning

def postprocess_dls(exp_folder, exp_name, cycle, n_384_reps, batch_size, reagent_names, acq_select, acq_param, cutoffs=None, plot_figures = True):
    
    # Select variable to be tested during machine learning (datalog_var) from DLS data
    # Select "comparison variable" if one DLS column can be used to filter the other

    #datalog_var = 'Mw-R (kDa)'
    #datalog_var = 'Range1 Radius (M) (nm)'
    datalog_var = 'Range1 Radius (I) (nm)'
    comp_var = 'Range1 %Mass (M)'
    #datalog_var = 'Radius (nm)'
    #datalog_var = 'Diffusion Coefficient (cm²/s)'
    #datalog_var = 'Polydispersity (nm)'

    # Configure the generation ID used to generate / locate experimental files
    if cycle == 1:
        exp_name_with_cycle = exp_name
    else:
        exp_name_with_cycle = f'{exp_name}_{cycle-1}'

    exp_name_next_cycle = f'{exp_name}_{cycle}'

    # Use DLS datalog or another DLS heuristic (true until said heuristic is developed)
    source_radius_from_datalog = True
    # Check if API was used to generate ACF file (set to False if manually downloaded from DYNAMICS)
    source_acf_from_api = True
    # Check if API was used to generate regularization peaks (hist) file (set to False if manually downloaded from DYNAMICS)
    source_hist_from_api = True

    # default filter values (comparison variable from datalog, ACF baseline value, ACF amplitude value, ACF SOS value)
    default_cutoffs = [75, 0.01, 0, 100]
    if cutoffs is None:
        cutoffs = default_cutoffs

    comp_var_cutoff = cutoffs[0]
    baseline_cutoff = [1-cutoffs[1],1+cutoffs[1]]
    amp_cutoff = cutoffs[2]
    sos_cutoff = cutoffs[3]

    comp_var_filter = False
    baseline_filter = True
    amp_filter = True
    sos_filter = False

    rd_state = 42
    nan_threshold = 3
    pca = PCA(n_components=2)

    # Configure which plots should be generated 
    if plot_figures:
        if cycle == 1:
            plot_gen_label = 'Seed Library'
        else:
            plot_gen_label = f'Active Learning Round {cycle-1}'
        plot_sample_distro = True
        plot_regression = True
        plot_shap = True
        plot_future_distro = True
    else:
        plot_sample_distro = False
        plot_regression = False
        plot_shap = False
        plot_future_distro = False
    
    # Load DLS datalog, ACF file and regularization peak (hist) file

    dls_table = os.path.join(exp_folder,f'{exp_name_with_cycle}_Datalog.csv')
    dls_table_data = pd.read_csv(dls_table, encoding= 'unicode_escape',index_col=0).dropna(axis=1, how='all').dropna(axis=0, how='all')
    dls_table_data.index = dls_table_data.index.str.strip()

    if source_acf_from_api:

        acf_file = os.path.join(exp_folder,f'{exp_name_with_cycle}_Raw.csv')
        acf_delay_file = os.path.join(exp_folder,f'{exp_name_with_cycle}_DelayTimes.csv')
        acf_file_data = pd.read_csv(acf_file, encoding= 'unicode_escape',index_col=None).dropna(axis=1, how='all')
        acf_file_data.columns = dls_table_data.index.tolist()
        acf_file_data.index = pd.read_csv(acf_delay_file, encoding= 'unicode_escape',index_col=0).index[:len(acf_file_data)]

    else:

        #acf_file_data = pd.read_excel(acf_file, index_col=None)
        acf_file = os.path.join(exp_folder,f'{exp_name_with_cycle}_Raw.csv')
        acf_file_data = pd.read_csv(acf_file, encoding= 'unicode_escape',index_col=0).dropna(axis=1, how='all')
        acf_file_data.columns = acf_file_data.columns.str.strip()

    if source_hist_from_api:

        hist_x_file = os.path.join(exp_folder,f'{exp_name_with_cycle}_xValues MR.csv')
        hist_y_file = os.path.join(exp_folder,f'{exp_name_with_cycle}_yValues MR.csv')
        hist_x_data = pd.read_csv(hist_x_file, encoding= 'unicode_escape', index_col=None).dropna(axis=1, how='all')
        hist_y_data = pd.read_csv(hist_y_file, encoding= 'unicode_escape', index_col=None).dropna(axis=1, how='all')
        hist_x_data.columns, hist_y_data.columns = dls_table_data.index.tolist(), dls_table_data.index.tolist()

        hist_file_data = pd.DataFrame()
        
        for x_col, y_col in zip(hist_x_data.columns, hist_y_data.columns):
            temp_dict = dict(zip(hist_x_data[x_col], hist_y_data[y_col]))
            temp_df = pd.DataFrame({x_col: temp_dict.values()}, index=temp_dict.keys())
            temp_df.index.name = 'Radius (nm)'
            hist_file_data = pd.concat([hist_file_data, temp_df], axis=1)

        hist_file_data = hist_file_data.sort_index()

    else:

        hist_file = os.path.join(exp_folder,f'{exp_name_with_cycle}_Peaks.csv')
        hist_file_data = pd.read_csv(hist_file, encoding= 'unicode_escape',index_col=0).dropna(axis=1, how='all')
        #hist_file_data = hist_file_data.dropna(how='all')
        hist_file_data.columns = hist_file_data.columns.str.strip()

    # Create list of empty lists for ACF and regularization peak data
    # Set index equal to names of wells collected during DLS
    acf_list = pd.Series([None]*len(dls_table_data))
    hist_list = pd.Series([None]*len(dls_table_data))
    acf_list.index = dls_table_data.index
    hist_list.index = dls_table_data.index

    for i in range(len(dls_table_data)):
        acf_list.iloc[i] = acf_file_data.iloc[:, i].tolist()
        hist_list.iloc[i] = hist_file_data.iloc[:, i].tolist()

    # Select variable column selected from DLS datalog to represent Rh
    if source_radius_from_datalog == True:
        var_label = datalog_var
        var_series = pd.to_numeric(dls_table_data[datalog_var].replace('--',np.nan))

    # Attempt to set alternate DLS heuristic (unfinished)
    elif source_radius_from_datalog == False:

        var_label = 'Hist Radius (nm)'
        rh_range = [0,16]
        hist_cutoff = 30
        var_series = pd.Series([None]*len(dls_table_data))

        for i,hist in enumerate(hist_list):
            abridged_hist = [hist[id] for id in [i for i,x in enumerate(hist_file_data[dls_table_data.index[i]].index) if x>rh_range[0] and x<=rh_range[1]]]
            if np.nanmax(abridged_hist)>hist_cutoff:
                #var_series.loc[i] = hist_file_data.index[abridged_hist.index(np.nanmax(abridged_hist))]
                var_series.loc[i] = np.nanmean(abridged_hist)

    # Set index of Rh values equal to names of wells collected during DLS
    var_series = pd.to_numeric(var_series)
    var_series.index = dls_table_data.index

    offset_96 = 0
    offset_384 = 0
    positions_96 = [rows_96[x%8]+cols_96[offset_96+x//8] for x in range(len(dls_table_data)//n_384_reps)]
    positions_384_by_well = map_96_to_384(positions_96, n_384_reps, 'per_well', rows_96, cols_96, offset_384)

    # Perform data filtering based on DLS datalog comparison variable
    if comp_var_filter==True:
        comp_var_series = pd.to_numeric(dls_table_data[comp_var].replace('--',np.nan))
        var_series = var_series.replace([x for x,y in zip(var_series, comp_var_series) if float(y)<comp_var_cutoff],np.nan)

    # Perform data filtering based on ACF baseline thresholds
    if baseline_filter==True:
        base_var_series = pd.to_numeric(dls_table_data['Baseline'].replace('--',np.nan))
        var_series = var_series.replace([x for x,y in zip(var_series, base_var_series) if (float(y)<baseline_cutoff[0] or float(y)>baseline_cutoff[1])],np.nan)

    # Perform data filtering based on ACF amplitude thresholds
    if amp_filter==True:
        amp_var_series = pd.to_numeric(dls_table_data['Amplitude'].replace('--',np.nan))
        var_series = var_series.replace([x for x,y in zip(var_series, amp_var_series) if float(y)<amp_cutoff],np.nan)

    # Remove ACFs and regularization peaks for polymers that were digitally filtered out
    for i in range(len(var_series)):
        if pd.isna(var_series.iloc[i]):
            acf_list.iloc[i] = []
            hist_list.iloc[i] = []

    # Compute metrics across sample replicates for each polymer: mean, SD, SEM, and select ACFs and regularization peaks for good acquisitions only
    means, stds, sems = pd.Series(dtype=object, index=positions_96), pd.Series(dtype=object, index=positions_96), pd.Series(dtype=object, index=positions_96)
    acfs, hists = pd.Series(dtype=object, index=positions_96), pd.Series(dtype=object, index=positions_96)

    for i, pos in enumerate(positions_96):

        indices = positions_384_by_well[i]
        non_nan_indices = sum([not math.isnan(x) for x in var_series[indices]])

        if non_nan_indices >= nan_threshold:
            means.loc[pos] = np.nanmean(var_series[indices])
            stds.loc[pos] = np.nanstd(var_series[indices])
            sems.loc[pos] = np.nanstd(var_series[indices])/np.sqrt(non_nan_indices)
            peak_choice=np.where(~np.isnan(var_series[indices]))[0]
            acfs.loc[pos] = acf_list[[indices[x] for x in peak_choice]]
            hists.loc[pos] = hist_list[[indices[x] for x in peak_choice]]


    dp_assay = False
    copolymer_assay = False
    active_learning_assay = True
    pd.options.mode.chained_assignment = None

    # Load polymer sample spreadsheet, set up dataset for ML
    setup_file = os.path.join(exp_folder,f'samples_{exp_name_with_cycle}.xlsx')
    setup_data = pd.read_excel(setup_file, index_col=None)

    # If DP variation only (discrete DP series)
    if dp_assay:

        dp_list = setup_data['Degree of Polymerization']
        design_space_dp_list = list(range(100,725,25))
        sel = slice(0,len(dp_list))
        dataset = pd.DataFrame({'DP': dp_list, 'Rh': means[sel], 'alpha Rh': stds[sel]},index=means.index[sel])
        # Represent design space as entire DP array
        x_range = np.array(design_space_dp_list).reshape(-1, 1)

    # If monomer variation only (constant DP, not included as feature)
    if copolymer_assay:

        monomer_table = get_monomers(setup_data)
        dataset = monomer_table.copy()
        dataset['Rh'] = means.values
        dataset['alpha Rh'] = stds.values ** 2

        # Used during current experiment to exclude DP series from the dataset
        # dataset.drop([40,41,42,43,44,45,46,47],inplace=True)

        # Design space generation for monomer variation only (note 2.5% thresholds as well as lack of DP component)
        options = [float(x) for x in list(np.arange(0,27.5,2.5))+list(np.arange(75,102.5,2.5))]
        x_range = []
        for c in product(options, repeat = 5):
            if sum(c) == 100 and all(x<=25 for x in c[1:]) and c[0]>=75 and sum(1 for x in c[1:] if x != 0) <= 3:
            x_range.append(c)
        x_range = [list(x) for x in x_range]

    # If active learning (desired use case of program)
    if active_learning_assay:

        dp_list = setup_data['Degree of Polymerization']
        monomer_table = get_monomers(setup_data)
        dataset = monomer_table.copy()
        dataset['DP'] = dp_list
        dataset['Rh'] = means.values
        # Variance array to be used by GPR as noise regularization parameter
        dataset['alpha Rh'] = stds.values ** 2
        
        # Plot 2D PCA decomposition of polymer features
        if plot_sample_distro:
            
            lhs_samples_reconstructed = np.hstack((monomer_table, np.array(dp_list).reshape(-1, 1)))
            samples_pca = pca.fit_transform(lhs_samples_reconstructed)
            fig = plt.figure(figsize=(8, 6)); plt.grid(True)
            plt.scatter(samples_pca[:, 0], samples_pca[:, 1], alpha=0.7, edgecolors='k')
            plt.title("PCA of LHS Samples (projected to 2D)")
            plt.xlabel("Principal Component 1"); plt.ylabel("Principal Component 2")
            plt.savefig(os.path.join(exp_folder,f'fig_sample_distro_{exp_name_with_cycle}.png'), dpi=fig.dpi)

        # Get design space to be used during active learning
        x_range = get_design_space()

    # Make sure no filtered polymers end up in the dataset
    dataset.reset_index(drop=True, inplace=True)
    dataset = dataset[dataset['Rh'].notna()]
    dataset['Rh'] = pd.to_numeric(dataset['Rh'], errors='coerce')
    dataset['alpha Rh'] = pd.to_numeric(dataset['alpha Rh'], errors='coerce')
    reordered_cols = reagent_names+['DP','Rh','alpha Rh']

    # Add blank column if dataset is completely missing any previously used monomer (i.e. later in active learning)
    for name in reagent_names:
        if name not in dataset.columns:
            dataset[name]=0

    # Configure file name to save dataset path depending on cycle number
    if cycle == 1:

        dataset_path = os.path.join(exp_folder,f'{exp_name}_dataset.xlsx')
        dataset[reordered_cols].to_excel(dataset_path, index=False)
        
    else:

        if cycle == 2:
            prior_data_path = os.path.join(exp_folder,f'{exp_name}_dataset.xlsx')
        else:
            prior_data_path = os.path.join(exp_folder,f'{exp_name}_{cycle-2}_dataset.xlsx')

        # For all cycles after seed library, append new samples to existing samples when creating dataset
        prior_data = pd.read_excel(prior_data_path, index_col=None).dropna(axis=1, how='all')
        dataset = pd.concat([prior_data,dataset], ignore_index=True).reset_index(drop=True)

        dataset_path = os.path.join(exp_folder,f'{exp_name_with_cycle}_dataset.xlsx')
        dataset[reordered_cols].to_excel(dataset_path, index=False)
        
    # Ensure that column order matches the order of the design space
    dataset = dataset[reordered_cols].fillna(0)

    # Machine learning method based on GPR
    def run_gpr(X, y, z, n_cv, rd_state, plot_regression):

        # Bounds for kernel hyper-parameter tuning
        signal_variances = [0.1, 1.0, 10.0]
        length_scales = [0.1, 1.0, 10.0]
        sigma_zeros = [0.1, 1.0, 10.0]
        noise_levels = [1e-5, 1e-2, 0.1]
        
        print(f'CV: {n_cv}')

        kernel_candidates = [
            ConstantKernel(variance)*RBF(lscale) + WhiteKernel(noise)
            for variance in signal_variances
            for lscale in length_scales
            for noise in noise_levels
        ]

        param_dist = {"kernel": kernel_candidates}
        #scores = ['r2','neg_mean_squared_error','neg_mean_absolute_error']
        scores = ['r2']

        # Kernel hyper-parameter tuning setup
        for score in scores:
            #print("# Tuning hyper-parameters for %s \n" % score)
            clf = RandomizedSearchCV(
                estimator=GaussianProcessRegressor(normalize_y=True),
                param_distributions=param_dist,
                scoring='%s' % score,
                n_iter=10,
                cv=n_cv,
                random_state=rd_state
            )

        # Split dataset into n_cv folds for cross-validation (CV)
        kf = KFold(n_splits=n_cv, shuffle=True, random_state=rd_state)
        y_tests, y_means, y_stds = np.array([]), np.array([]), np.array([])#, np.array([]), np.array([])

        # Normalize data
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Initialize figure to plot regression results
        if plot_regression:
            fig, axs = plt.subplots(1,n_cv+1,figsize=(n_cv*4,2))
            # fig.suptitle(f'{datalog_var}: Cross-validation', y=1.1)

        # Iterate through each KFold and perform a train-test split
        for i, (train_index, test_index) in enumerate(kf.split(X,pd.qcut(y, q=n_cv, labels=False))):
            X_train, X_test = X_scaled[train_index], X_scaled[test_index]
            y_train, y_test = y.iloc[train_index], y.iloc[test_index]
            alpha_train = np.array(z.iloc[train_index])

            # Select best hyper-parameter models and train the model
            clf.fit(X_train, y_train)
            gpr = clf.best_estimator_
            gpr.alpha = alpha_train
            gpr.fit(X_train, y_train)

            # Test the model and evaluate R^2 accuracy
            y_mean, y_std = gpr.predict(X_test, return_std=True)
            y_tests = np.concatenate((y_tests,y_test))
            y_means = np.concatenate((y_means,y_mean))
            y_stds = np.concatenate((y_stds,y_std))
            fold_r2 = r2_score(y_test, y_mean)

            # Make fold-specific regression subplot
            if plot_regression:
                ax = axs[i]
                ax.scatter(y_test, y_mean, s=10)
                ax.set(xlabel=f'Actual', ylabel=f'Predicted', title=f'R2: {fold_r2:.3f}')
                ax.plot(np.linspace(0,10,100),np.linspace(0,10,100),color='k',alpha=0.3)
                ax.set_xlim(left = 0)
                ax.set_ylim(bottom=0)
                ax.grid(color='gainsboro')
                ax.set_facecolor('w')
            
            print(f'R2 for fold {i+1}: {fold_r2}')

        # Fit to entire dataset, train and test
        clf.fit(X_scaled, y)
        gpr = clf.best_estimator_
        gpr.alpha = np.array(z)
        gpr.fit(X_scaled, y)
        r2 = r2_score(y_tests, y_means)

        if plot_regression:
            ax = axs[-1]
            ax.scatter(y_tests, y_means, s=10)
            ax.set(xlabel=f'Actual', ylabel=f'Predicted', title=f'Total R2: {r2:.3f}')
            ax.plot(np.linspace(0,10,100),np.linspace(0,10,100),color='k',alpha=0.3)
            ax.set_xlim(left = 0)
            ax.set_ylim(bottom=0)
            ax.grid(color='gainsboro')
            ax.set_facecolor('w')
            plt.savefig(os.path.join(exp_folder,f'fig_r2_reg_{exp_name_with_cycle}.png'), dpi=fig.dpi)

        return gpr, scaler, r2

    # Call ML method using n_cv cross-validation folds
    n_cv = 4
    X, y, z = prep_data_for_gpr(dataset, n_cv, rd_state)
    gpr, scaler, r2 = run_gpr(X, y, z, n_cv, rd_state, plot_regression)
    print(f'R2 after {n_cv} folds: {r2}')

    # Define EI activation function
    # input eta = exploration parameter
    def expected_improvement(y_pred, y_std, best_y, eta=1.0):
        z = (y_pred - best_y + eta) / y_std
        return (y_pred - best_y + eta) * stats.norm.cdf(z) + y_std * stats.norm.pdf(z)

    # Define LCB activation function
    # input lambda = exploration parameter
    def lcb(y_pred, y_std, lam=1.0):
        return y_pred - lam * y_std

    # Method that uses acquisition function (acq_select) with select exploit-explore parameter (acq_param)
    # input x_range = design space to select polymers from
    # input batch_size = number of new polymers to get selected

    def choose_batch(model, scaler, x_range, X, y, z, y_pred, y_std, acq_select, acq_param, batch_size):

        new_ids = []
        new_forms = np.empty((0,x_range.shape[1]), float)

        X_array = X.to_numpy().astype(float)

        # Mask that helps ensure no polymers already in the dataset are selected again
        unique_mask = ~np.any(np.all(np.isclose(x_range[:, None, :], X_array[None, :, :], atol=1e-8), axis=2), axis=1)

        # Prior to each polymer, use Kriging Believer strategy to find which formulation minimizes the acquisition function
        for i in range(batch_size):

            # Compute acquisition function depending on which function is selected
            if acq_select == 'ei' or acq_select == 'EI':
                acq = expected_improvement(y_pred, y_std, best_y=np.min(y), eta=acq_param)
                print(f'EI performed using eta = {acq_param}')
            elif acq_select == 'lcb' or acq_select == 'LCB':
                acq = lcb(y_pred, y_std, lam=acq_param)
                print(f'LCB performed using lambda = {acq_param}')

            # Mark design space formulations used previously in the dataset,
            # as well as polymer formulations already selected during current round of batch selection,
            # as used IDs not to be repeatedly selected
             
            biased_acq = acq.copy()
            used_ids = np.append(np.where(~unique_mask)[0],new_ids).astype(int)
            #used_ids = new_ids
            for id in used_ids:
                biased_acq[id] = np.inf

            # Identify scaled formulation ID that minimizes acquisition function 
            new_id = np.argmin(biased_acq)
            print(f'Sample {i+1} acq func value: {np.min(biased_acq)}')
            new_ids.append(new_id)

            # Use scaled formulation ID to select true formulation
            new_form = x_range[new_id]
            new_forms = np.append(new_forms, np.array([new_form]), axis=0)

            # Use GPR trained on entire design space to predict the model uncertainty for the polymer in advance
            scaled_new_id =  scaler.transform(x_range)[new_id]
            _, new_std = model.predict(scaled_new_id.reshape(1,-1), return_std = True)
            
            # Add new formulation, as well as hypothetical Rh predicted by design-space trained GPR model,
            # as additional point to the model dataset
         
            X = np.vstack((X, new_form.reshape(1,-1)))
            y = np.hstack((y, y_pred[new_id].reshape(-1)))

            # Use predicted uncertainty value squared as variance to be added to GPR for noise regularization 
            z = np.append(z, (new_std**2).reshape(-1))

            # Add new variance array as GPR parameter
            model = GaussianProcessRegressor(
                kernel=model.kernel,
                alpha=z,
                random_state=rd_state, 
                normalize_y=True
            )

            # Re-train the model using new formulation,
            # as though the newest batch selected polymer's 
            # in-silico predicted Rh was already an experimental datapoint.
            # This will be used to update the acquisition function 
            # once it is recalculated to obtain the next batch selected polymer
            model.fit(scaler.transform(X), y)

            # Use the new model to predict Rh values across the entire design space
            y_pred, y_std = model.predict(scaler.transform(x_range), return_std=True)

        return new_forms

    best_idx = np.argmin(y)
    best_x = X.iloc[best_idx]
    best_y = y.iloc[best_idx]
    print(best_idx)
    print(best_x)
    print(best_y)

    # Predict 'mean' and SD of Rh across entire design space using model trained on experimental data
    x_range_scaled = scaler.fit_transform(x_range)
    range_pred, range_std = gpr.predict(x_range_scaled, return_std=True)

    # Using batch selection, select the several samples for the next active learning generation
    new_formulations = choose_batch(gpr, scaler, x_range, X, y, z, range_pred, range_std, acq_select, acq_param, batch_size)

    # Use current model to generate in-silico predictions for the batch selected samples
    acq_pred, acq_std = gpr.predict(scaler.transform(new_formulations), return_std=True)
    acq_pred_table = pd.DataFrame(np.hstack((new_formulations.astype(int), acq_pred.reshape(-1,1), acq_std.reshape(-1,1))),columns=reordered_cols)
    acq_pred_path = os.path.join(exp_folder,f'{exp_name_next_cycle}_predicted_dataset.xlsx')
    acq_pred_table.to_excel(acq_pred_path, index=False)

    # Plot PCA decomposition to represent spread of batch selected polymers in current design space
    if plot_future_distro:

        pca_X = pca.fit_transform(X)
        pca_new_forms = pca.transform(new_formulations)
        fig = plt.figure(figsize=(8,6)); plt.grid(True)
        plt.scatter(pca_X[:, 0], pca_X[:, 1], alpha=0.7, edgecolors='k')
        plt.scatter(pca_new_forms[:, 0], pca_new_forms[:, 1], edgecolors='k', color='red')
        plt.title(f'{plot_gen_label}: Sample distribution\n(projected to 2D)')
        plt.xlabel('Principal Component 1'); plt.ylabel('Principal Component 2')
        plt.savefig(os.path.join(exp_folder,f'fig_predicted_distro_{exp_name_next_cycle}.png'), dpi=fig.dpi)

    # Plot SHAP analysis of design space predictions
    if plot_shap:
        fig = plt.figure(figsize=(4,3))
        explainer = shap.Explainer(gpr.predict, x_range_scaled, feature_names = X.columns)
        shap_vals = explainer(x_range_scaled)
        print("SHAP shape:", shap_vals.values.shape)
        print("SHAP min/max:", shap_vals.values.min(), shap_vals.values.max())
        print("SHAP sum (for sanity):", shap_vals.values.sum())
        shap.plots.beeswarm(shap_vals, show=False)
        plt.title(f'SHAP Analysis: {plot_gen_label}')
        plt.savefig(os.path.join(exp_folder,f'fig_shap_{exp_name_with_cycle}.png'), dpi=fig.dpi)

    # Call sampling function to format batch selection polymers into table
    # that can be used in the next active learning generation
    return create_polymer_table(new_formulations, reagent_names, exp_name_next_cycle)