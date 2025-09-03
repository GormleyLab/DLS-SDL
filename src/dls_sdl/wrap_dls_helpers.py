import os, time, re
from alive_progress import alive_bar
from collections import Counter
from itertools import count
import pandas as pd

from wrap_dls_connect import *
from wrap_dls_dicts import *

#Converts label to DYNAMICS enum value
#inputs: 'x', label string, and 'dict', the dictionary representing the enum
#outputs: 'correct_key', DYNAMICS enum class item corresponding to the label string

def val2key(x, dict):
    d2l_keys, d2l_vals = list(dict.keys()), list(dict.values())
    correct_key = d2l_keys[d2l_vals.index(x)]
    return(correct_key)

#Truncates DYNAMICS enum class into dictionary for ease of access / printing
#inputs: 'x', the enum class, and 'n', the number of items to extract from the class
#outputs: 'info' dictionary containing the proper attributes

def unravel_obj(x,n):

    # Returns last n keys of a DYNAMICS enum
    keys = sorted((dir(x))[-n:])
    
    info=dict.fromkeys(keys)
    for k in keys:
        info[k]=getattr(x,k)
    return(info)

#----------DLS COMMAND HELPERS----------

#Prints label followed by array elements if the array is not empty
def print_array(x, label):
    x=[i for i in x]
    if len(x)>0: 
        print(label)
        print(x)

#Converts well names from 'A01' to 'A1', removing leading zero (MoveToWell)
def remove_leading_zero(well):
    if well==None:
        return None
    well_alpha = re.findall(r'[^\W\d]+',well)
    well_num = re.findall(r'\d+',well)
    well = well_alpha[0]+str(int(well_num[0]))
    return well


#Displays a progress bar while the door is being opened / closed (ControlDoor)
#inputs: 'n' number of seconds that the door is in motion
#outputs: visual display

def dls_door_bar(n):
    res=1000
    with alive_bar(res) as bar:
        for i in range(res):
            time.sleep(n/res)
            bar()


#Converts array of solvent information into a SolventInfo Enum object
#inputs: args list formatted as [str name, str model, float refractive_index, float viscosity, float viscosity_temp]
#outputs: 'solv', SolventInfo Enum with model string converted to a SolventTemperatureModel Enum attribute

def arr2solv(args):
    if len(args)!=5:
        pass #results in error
    
    solv = get_dynamics_lib().SolventInfo()
    solv.name = args[0]
    solv.model = val2key(args[1], dict_tempmdl)

    if not isinstance(args[2], float):
        print('Refractive index must be a number')
        solv.refractiveIndexAt589nm20C = float(input('Refractive index: '))
    else:
        solv.refractiveIndexAt589nm20C = args[2]

    if not isinstance(args[3], float):
        print('Viscosity must be a number')
        solv.viscosity = float(input('Viscosity: '))
    else:
        solv.viscosity = args[3]

    if not isinstance(args[4], float):
        print('Viscosity Temperature must be a number')
        solv.viscosityTemp = float(input('Viscosity temperature: '))
    else:
        solv.viscosityTemp = args[4]
    
    return solv

#Edits path title to avoid overwriting existing file (adds a (1) or (n) if the file is found in the directory)
#inputs: file, directory
#outputs: new file title
def check_path_dupe(path):
    base, ext = os.path.splitext(path)
    path_n = 1
    while os.path.exists(path):
        path = f"{base} ({path_n}){ext}"
        path_n += 1
    return path


#----------MAIN HELPERS--------
#Set up output folder for exported data
def export_path_setup(exp_folder, export_title):
    return check_path_dupe(os.path.join(exp_folder, export_title))

#Set up output folders and locates correct preset file
def dls_path_setup(script_dir, data_dir, exp_name, plate_type, acqTime):
    preset_folder = 'DLS Presets'
    preset_title = 'preset_'+str(plate_type)+'_'+str(acqTime)+'s.dpst'
    plate_title = 'plate_setup.xlsx'

    exp_title = exp_name+'.dexp'
    export_title = exp_name+'_Datalog.csv'

    preset_path = os.path.join(script_dir, preset_folder, preset_title)
    plate_path = os.path.join(script_dir, plate_title)
    if not os.path.exists(preset_path):
        print('Please create a DYNAMICS preset file containing the correct acquisition time!')
    
    # Checks if experiment name contains a cycle post-script (lets method run without cycle arg)
    if exp_name.rfind('_') != -1:
        if exp_name[exp_name.rfind('_') + 1:].isdigit():
            base_expname = '_'.join(exp_name.split('_')[:-1])
        else:
            base_expname = exp_name
    else:
        base_expname = exp_name

    #exp_folder = os.path.join(script_dir, data_subfolder, '_'.join(exp_name.split('_')[:-1]))
    exp_folder = os.path.join(data_dir, base_expname)
    if not os.path.exists(exp_folder):
        os.makedirs(exp_folder, exist_ok = True)

    save_path = check_path_dupe(os.path.join(exp_folder, exp_title))
    export_path = export_path_setup(exp_folder, export_title)
    
    return preset_path, plate_path, exp_folder, save_path, export_path

#Convert ACF / Regu data into preliminary dataframes
#df: output for GetACFData / GetReguData converted to Pandas DataFrame
#m: current measurement
#out: output file for ACF or Regu data
#outputs: DataFrame with corrected column labels

def add2frame(df, m, out):
    new_cols  = ['Meas {m} '.format(m=m+1)+col for col in df.columns]
    df = df.set_axis(new_cols, axis=1)
    if out.empty:
        return df
    else:
        return pd.concat([out,df], axis=1)

#Saves .xlsx Excel spreadsheets organized by collected ACF / Regu data category
#inputs: 'data', ACF / Regu output for all measurements, 'cats' list of category strings, 'title' experiment name string
#outputs: saved .xlsx file(s) sorted by data category

def sort_write(data, cats, exp_folder, title):
    for _, substr in enumerate(cats):
        df = data[[c for c in data.columns if substr in c]].copy()
        df.dropna(axis = 0, how = 'all', inplace = True)
        if not df.empty:
            path = export_path_setup(exp_folder, '_'.join((title,substr))+'.csv')
            df.to_csv(path, index=False)

#Returns list of .xlsx cells that contain samples
#inputs: 'file' path of .xlsx template file
#outputs: 'list' of sample names formatted as Excel cell names

def xls2measurement(filename, unique=False):

    df = pd.read_excel(filename, 'Sheet1', nrows=16, usecols = 'C:Z')
    df.index = [chr(ord('a') + x).upper() for x in df.index]
    df.dropna(axis = 0, how = 'all', inplace = True)
    df.dropna(axis = 1, how = 'all', inplace = True)
    cells = df.where(df != 'Blank').stack()
    names, values = cells.index.tolist(), cells.values.tolist()
    names = [x[0]+str(x[1]) for x in names]
    
    if unique==False:
        values = names
    if unique==True:
        iters = {k: count(1) for k, v in Counter(values).items() if v > 1}
        values = [x+f' ({str(next(iters[x]))})' if x in iters else x for x in values]
    
    return names, values