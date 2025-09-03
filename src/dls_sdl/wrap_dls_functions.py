import inspect
import numpy as np

from wrap_dls_connect import *
from wrap_dls_helpers import *
from wrap_dls_dicts import *

#Prints ResultCode only if an error occurs
def print_dls_result(x):
    if x != resultcode.WRESULT_OK:
        print(x)

#----------BASIC RUN----------

#RETURN
def load_dls_preset(preset):
    dynapro = get_dynapro(get_dynamics_lib())
    _, handle, result = dynapro.LoadPreset(preset)
    #handle, result = dynapro.LoadPreset(preset)
    print_dls_result(result)
    return handle #UInt32()

#RETURN
def load_dls_exp(path):
    dynapro = get_dynapro(get_dynamics_lib())
    handle, messages, result = dynapro.LoadExperiment(path)
    print_dls_result(result)
    print_array(messages, inspect.stack()[0][3]+' messages: ')
    return handle #UInt32()

def save_dls_exp(handle, path):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.SaveExperiment(handle, path)
    print_dls_result(result)
    print(f'Experiment saved to {path}')

def close_dls_exp(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.CloseExperiment(handle)
    print_dls_result(result)

def load_dls_plate_templates(handle, path): # earlier version only had handle input
    dynapro = get_dynapro(get_dynamics_lib())
    messages, result = dynapro.LoadPlateTemplates(handle, path)
    print_dls_result(result)
    print_array(messages, inspect.stack()[0][3]+' messages: ')

def connect_to_dls(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    result, errorMsg = dynapro.ConnectToHardware(handle)
    if result == resultcode.WRESULT_EXP_ALEADY_CONNECTED:
        print("Connection errors: ", result)
        return result
    elif errorMsg: 
        print("Connection errors: ", errorMsg)
        return errorMsg

def disconnect_from_dls(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    result, errorMsg = dynapro.DisconnectFromHardware(handle)
    print_dls_result(result)
    if errorMsg!=None: print("Disconnection errors: ", errorMsg)

def set_dls_laser(handle, power):
    dynapro = get_dynapro(get_dynamics_lib())
    if not (isinstance(power, int) and power>=0 and power<=100):
        print('Laser power setting must be number between 0 and 100')
        power=int(input('Laser power: '))
    result = dynapro.SetLaserPower(handle, power)
    print_dls_result(result)

#RETURN
def get_dls_laser(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    power, result = dynapro.GetLaserPower(handle)
    print_dls_result(result)
    print("Power: "+str(power))
    return power

def start_dls_exp(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.StartRecordingData(handle)
    print_dls_result(result)
    print('start')

def stop_dls_exp(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.StopRecordingData(handle)
    print_dls_result(result)
    print('stop')

def pause_dls_exp(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.SuspendRecordingData(handle) 
    print_dls_result(result)
    print('pause')

def resume_dls_exp(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.ResumeRecordingData(handle)
    print_dls_result(result)
    print('resume')



#----------FORMATTING----------

def set_dls_exp_title(handle, string):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.SetExperimentTitle(handle, string)
    print_dls_result(result)

#RETURN
def get_dls_exp_title(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    string, result = dynapro.GetExperimentTitle(handle)
    print_dls_result(result)
    print(string)
    return string

def set_dls_exp_comment(handle, string):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.SetExperimentComment(handle, string)
    print_dls_result(result)

#RETURN
def get_dls_exp_comment(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    string, result = dynapro.GetExperimentComment(handle)
    print_dls_result(result)
    print(string)
    return string

def set_dls_meas_solv(handle, meas_index, solvent_name): # earlier version only had handle input
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.SetMeasSolvent(handle, meas_index, solvent_name)
    print_dls_result(result)

def set_dls_meas_name(handle, meas_index, meas_name):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.SetMeasName(handle, meas_index, meas_name)
    print_dls_result(result)

def set_dls_meas_val(handle, meas_index, name, value):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.SetMeasUserDefinedValue(handle, meas_index, name, value)
    print_dls_result(result)



#----------SOLVENTS AND INFO----------

#NO RESULT
def get_dls_software_version():
    dynapro = get_dynapro(get_dynamics_lib())
    vers = dynapro.GetSoftwareInfo()
    print(unravel_obj(vers,2))

def add_dls_solvent(handle, solvent):
    dynamicsLib = get_dynamics_lib() 
    dynapro = get_dynapro(dynamicsLib)
    if isinstance(solvent, dynamicsLib.SolventInfo):
        pass
    elif isinstance(solvent, list):
        solvent = arr2solv(solvent)
    _, result = dynapro.AddSolvent(handle, solvent)
    print_dls_result(result)

def modify_dls_solvent(handle, solvent):
    dynapro = get_dynapro(get_dynamics_lib())
    _, result = dynapro.ModifySolvent(handle, solvent)
    print_dls_result(result)

def delete_dls_solvent(handle, solvent):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.DeleteSolvent(handle, solvent)
    print_dls_result(result)

def get_dls_hardware_info(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    hw, result = dynapro.GetHardwareInfo(handle)
    print_dls_result(result)
    print(unravel_obj(hw,12))

#RETURN; NO RESULT
def get_dls_solvent_names(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    solvent_names = dynapro.GetSolventNames(handle)
    # Converts a newline-delimited string to a list
    return solvent_names.rstrip().split('\n')

#RETURN
def get_dls_solvent(handle, solvent_name):
    dynapro = get_dynapro(get_dynamics_lib())
    solv, result = dynapro.GetSolvent(handle, solvent_name)
    print_dls_result(result)
    info = unravel_obj(solv,5)
    info['model']=handle_enum(info['model'], tempmdl)
    print(info)
    return(solv)


#----------DATA COLLECTION----------

#RETURN
def get_dls_datalog_val(handle, meas_index, valueName):
    dynapro = get_dynapro(get_dynamics_lib())
    v, result = dynapro.GetDatalogValue(handle, meas_index, valueName)
    print(valueName+": "+str(v))
    print_dls_result(result)

#RETURN
def export_dls_datalog(handle, meas_index, acqIndex, exportPath):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.ExportDatalogTable(handle, meas_index, acqIndex, exportPath) #exportData
    print_dls_result(result) #if meas_index=-1, acqIndex doesn't matter

#RETURN
def get_dls_acf_data(handle, string, meas_index, acqIndex):
    dynamicsLib = get_dynamics_lib()
    dynapro = get_dynapro(dynamicsLib)
    choices = dynamicsLib.ACFChoices()
    if 'D' in string:
        choices.raw=1
    if 'C' in string:
        choices.cumulants=1
    if 'R' in string:
        choices.regularization=1

    #np.set_printoptions(formatter={'float_kind':"{:.3f}".format})
    pre_acf = dynapro.GetACFData(handle, meas_index, acqIndex, choices)
    print_dls_result(pre_acf[-1])
    pre_acf=list(pre_acf)
    pre_acf.pop(-1)

    acf_dict_names=['DelayTimes', 'Raw', 'Cumulants', 'Regu']
    acf_arrays = dict.fromkeys(acf_dict_names)
    for i, arr in enumerate(pre_acf):
        if arr != None:
            arr=np.array(arr)
            acf_arrays[acf_dict_names[i]]=arr

    return acf_arrays

#RETURN
def get_dls_regu_data(handle, dist_str, xstr, meas_index, acq_index):
    dynapro = get_dynapro(get_dynamics_lib())
    dist=None #distribution 
    rdist_letters=['I','M','N']
    if any(x in dist_str for x in rdist_letters):
        dist = handle_enum(dist_str, regudist)
    else: 
        print("Distribution type must be one of the following:\n'I' (intensity)\n'M' (mass)\n'N' (number)")
        dist = input('Enter the distribution type: ')
    
    xtype=None #x-axis setting
    xtype_letters=['C','D','R','T']
    if any(x in xstr for x in xtype_letters):
        xtype = handle_enum(xstr, reguxval)
    else:
        print("X-axis type must be one of the following:\n'C' (diffusion coefficient)\n'D' (diameter)\n'R' (radius)\n'T' (decay time)")
        xstr=input('Enter the X-axis type: ')

    #np.set_printoptions(formatter={'float_kind':"{:.3f}".format})
    x_vals, y_vals, result = dynapro.GetReguData(handle, meas_index, acq_index, dist, xtype)
    print_dls_result(result)
    
    regu_dict = {'xValues '+dist_str+xstr: np.array(x_vals), 'yValues '+dist_str+xstr: np.array(y_vals)}
    return regu_dict

def collect_dls_meas(handle, num_acqs, capture_images, force):
    dynapro = get_dynapro(get_dynamics_lib())
    if isinstance(num_acqs, int) == False or num_acqs<1:
        print("Number of acquisitions must be a positive number")
        num_acqs=int(input("How many acquisitions? "))
    result = dynapro.CollectAcquisitions(handle, num_acqs, capture_images, force)
    print_dls_result(result)

#----------TEMP / RATE----------

def set_dls_temp(handle,temp):
    dynapro = get_dynapro(get_dynamics_lib())
    if not isinstance(temp, float):
        print('Temperature must be a number')
        temp = float(input("Set the instrument's target temperature: "))
    result = dynapro.SetTemperature(handle, temp)
    print_dls_result(result)

def set_dls_ramp(handle,rate):
    dynapro = get_dynapro(get_dynamics_lib())
    if not isinstance(rate, float):
        print('Ramp rate must be a number')
        rate = float(input("Set the instrument's temperature ramp rate: "))
    result = dynapro.SetRampRate(handle, rate)
    print_dls_result(result)

def set_dls_temp_and_ramp(handle,temp,rate):
    dynapro = get_dynapro(get_dynamics_lib())
    if not isinstance(temp, float):
        print('Temperature must be a number')
        temp = float(input("Set the instrument's target temperature: "))
    if not isinstance(rate, float):
        print('Ramp rate must be a number')
        temp = float(input("Set the instrument's temperature ramp rate: "))
    result = dynapro.SetTemperatureAndRampRate(handle, temp, rate)
    print_dls_result(result)

#RETURN
def get_dls_temp(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    temp, rate, result = dynapro.GetTemperature(handle)
    print_dls_result(result)
    print(str(temp)+' °C, '+str(rate)+' °C/min')
    return temp, rate

#RETURN
def get_dls_sample_temp(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    temp, result = dynapro.GetSampleTemperature(handle)
    print_dls_result(result)
    print(str(temp)+' °C')
    return temp

#RETURN
def get_is_dls_temp_locked(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    locked, temp, result = dynapro.CheckTemperatureLock(handle)
    print_dls_result(result)
    print(str(locked)+' at: '+str(temp)+' °C')
    return locked, temp



#----------ATTENUATION----------

def enable_dls_auto_atten(handle, enable):
    dynapro = get_dynapro(get_dynamics_lib())
    if enable in ['y','Y',1,True]: enable=True
    elif enable in ['n','N',0,False]: enable=False
    else:
        print("Select 'y' or 1 or True for 'enable'. Select 'n' or 0 or False for 'disable'.")
        enable=input('Enable auto-attenuation? ')
    result = dynapro.EnableAutoAttenuation(handle, enable)
    print_dls_result(result)

def get_is_dls_auto_atten_enabled(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    enable, result = dynapro.AutoAttenuationEnabled(handle)
    print_dls_result(result)
    print(enable)

def get_is_dls_auto_atten(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    attenuating, result = dynapro.IsAutoAttenuating(handle)
    print_dls_result(result)
    print(attenuating)

#RETURN
def get_dls_atten_level(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    level, result = dynapro.GetAttenuationLevel(handle)
    print_dls_result(result)
    print(level)
    return(level)

def gst_dls_atten_level(handle, level):
    dynapro = get_dynapro(get_dynamics_lib())
    if not (isinstance(level, int) and level>=0 and level<=100):
        print('Attenuation level must be a number between 0 and 100')
        level=int(input('Laser power: '))
    result = dynapro.SetAttenuationLevel(handle, level)
    print_dls_result(result)



#----------WELLS / DOOR----------

#RETURN
def pos_dls(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    wellName, result = dynapro.GetCurrentWell(handle)
    print_dls_result(result)
    return wellName

# Used as part of async_dls_move
def move_dls(handle, well):
    dynapro = get_dynapro(get_dynamics_lib())
    dynapro = get_dynapro(get_dynamics_lib())
    well = remove_leading_zero(well)
    result = dynapro.MoveToWell(handle, well)
    print_dls_result(result)

def open_dls_door(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.ControlDoor(handle, 1) #Int32(1)
    print_dls_result(result)
    #lag = ~24.25 #seconds to open the door
    #doorBar(lag)

def close_dls_door(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.ControlDoor(handle, 0) #Int32(0)
    print_dls_result(result)
    #lag = ~16.25 #seconds to close the door
    #doorBar(lag)

#RETURN
def get_dls_door_state(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    door, result = dynapro.GetDoorState(handle)
    print_dls_result(result)
    #print(door)
    return(door)

#RETURN
"""
def getDoorFullyOpen(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    door, result = dynapro.GetDoorFullyOpen(handle)
    print_dls_result(result)
    print(door)
    return(door)
"""

#----------SAMPLES----------

#RETURN
def get_dls_sample_name(handle):
    dynapro = get_dynapro(get_dynamics_lib())
    sample_names, result = dynapro.GetSampleNames(handle)
    print_dls_result(result)
    # Converts a newline-delimited string to a list
    return sample_names.rstrip().split('\n')

#RETURN
def get_dls_sample(handle, sample_name):
    dynapro = get_dynapro(get_dynamics_lib())
    samp, result = dynapro.GetSample(handle, sample_name)
    print_dls_result(result)
    info = unravel_obj(samp,8)
    info['mwrModel']=handle_enum(info['mwrModel'], mwrmdl)
    info['rgModel']=handle_enum(info['rgModel'], rgmdl)
    print(info)
    return samp

def set_dls_sample(handle, meas_index, sample_name):
    dynapro = get_dynapro(get_dynamics_lib())
    result = dynapro.SetMeasSample(handle, meas_index, sample_name)
    print_dls_result(result)
    print("Measurement "+str(meas_index)+" set to "+sample_name)