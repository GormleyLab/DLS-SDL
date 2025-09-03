from wrap_dls_connect import *
from wrap_dls_dicts import *
from wrap_dls_helpers import *
from wrap_dls_functions import *
from wrap_dls_events import *

# Python wrapper for DYNAMICS SDK 3.0.1 (DLS)
# Allows for pythonization of instrument run control

class DynamicsDLS:

    # Initialize control server tracking boolean and connection state
    def __init__(self):
        self.dynapro = None
        self.connected = False

    # Connect to DYNAMICS SDK dynamicsLib, then to control server dynapro
    # Bind events to control server to enable asynchronous event tracking
    def initialize(self):
        dynamicsLib = get_dynamics_lib()
        self.dynapro = get_dynapro(dynamicsLib)
        bind_events(self.dynapro, events)
        self.connected = True

    # Upon closing an experimental handle, DYNAMICS will have already closed
    # Unbind events from control server to disable asynchronous event tracking
    def shutdown(self):
        if self.dynapro:
            unbind_events(self.dynapro, events)
            self.dynapro = None
            self.connected = False

    # Open DLS door
    def opendoor(self, step):

        print('Running step %d' % step)

        #Set up output folders and locates correct preset file
        script_dir = os.getcwd()
        preset_folder = 'DLS Presets'
        # preset_title = 'preset_'+str(plate_type)+'_'+str(acqTime)+'s.dpst'
        open_title = 'Door Control Dummy.dexp'
        open_path = os.path.join(script_dir, preset_folder, open_title)
        print(f'Temporary path: {open_path}')

        # Open door by loading preset file and connecting to instrument
        handle = load_dls_exp(open_path)
        errorMsg = connect_to_dls(handle)
        if errorMsg: return

        open_dls_door(handle)
        # time.sleep plus cool visual
        dls_door_bar(24.25)
        disconnect_from_dls(handle)
        close_dls_exp(handle)


    # CLose DLS door
    def closedoor(self, step):
        
        print('Running step %d' % step)

        #Set up output folders and locates correct preset file
        script_dir = os.getcwd()
        preset_folder = 'DLS Presets'
        # preset_title = 'preset_'+str(plate_type)+'_'+str(acqTime)+'s.dpst'
        close_title = 'Door Control Dummy.dexp'
        close_path = os.path.join(script_dir, preset_folder, close_title)
        print(f'Temporary path: {close_path}')

        # Open door by loading preset file and connecting to instrument
        handle = load_dls_exp(close_path)
        errorMsg = connect_to_dls(handle)
        if errorMsg: return

        close_dls_door(handle)
        # time.sleep plus cool visual
        dls_door_bar(16.25)
        disconnect_from_dls(handle)
        close_dls_exp(handle)

    async def run(self, step, script_dir, data_dir, exp_name, dls_wells, plate_type = 384,
    acq_time = 5, num_acqs = 8, start_temp = 25.0, laser_power = 15, atten = True, acf_choices = 'D',
    dist_type = 'M', regXVal= 'R', capture_images = 0, save_per_acq = True):
        
        print('Running step %d' % step)

        preset_path, plate_path, exp_folder, save_path, export_all_meas = dls_path_setup(
                                script_dir, data_dir, exp_name, plate_type, acq_time)

        # Converts string of measurement wells into list
        measList = dls_wells.split(",")
        #measList, sampleNames = xls2measurement(plate_path, False)
        sampleNames = measList
        print(measList)

        # Convert regularization settings to list if multiple peak types (e.g. intensity and mass) are desired
        if len(dist_type) > 1:
            dist_type = dist_type.split(",")
        
        # Convert regularization settings to list if multiple data types (e.g. radius and diff. coeff.) are desired
        if len(regXVal) > 1:
            regXVal = regXVal.split(",")

        # Setting that truncates printed decimals to 3 digits
        np.set_printoptions(formatter={'float_kind':"{:.3f}".format})

        # Start experiment by loading correct preset file and connecting to instrument
        handle = load_dls_preset(preset_path)
        errorMsg = connect_to_dls(handle)
        if errorMsg: return

        # Position the plate reader so that it's calibrated at position A1,
        # then start laser, laser auto-attenuation and temperature stabilization

        await asyncio.sleep(20)
        await async_dls_move(handle, 'A1')
        set_dls_laser(handle, int(laser_power))
        if atten == True:
            enable_dls_auto_atten(handle, 1)
        else:
            enable_dls_auto_atten(handle, 0)
        set_dls_temp(handle, float(start_temp))

        # DO NOT DECREASE the stabilization time. Otherwise the DLS won't move wells
        stable_time = 30
        print(f'Waiting {stable_time} seconds to stabilize')
        await asyncio.sleep(stable_time)
        await wait_for_stable_dls_temperature(handle, target_temp=start_temp)

        # Begin data collection. 'acf_output' and 'regu_output' store the auto-correlation
        # function and regularization function data for each measured well

        print(f'Begin experiment {exp_name}')
        acf_all_meas, regu_all_meas = pd.DataFrame(), pd.DataFrame()
        old_acf_cols, old_regu_cols = [], []

        # For each measurement, adjust the plate so that the correct well is above the laser.
        # Collect the right number of acquisitions at the correct acuiqistion time,
        # then save the data and generate a Pandas dataframe for each ACF / regularization output
        for m, meas in enumerate(measList):

            current_temp, current_rate = get_dls_temp(handle)
            print(f'Current Temp: {float(current_temp):.2f} °C; Current Rate: {float(current_rate):.3f} °C/s')
            await asyncio.sleep(1)
            
            print(f'Moving to well {meas} for measurement {m+1}')                 
            await async_dls_move(handle, meas)

            tracker.reset(expected=num_acqs)
            collect_dls_meas(handle, num_acqs, capture_images, force=0)
            print(f'Waiting for {num_acqs} acquisitions to complete')
            await wait_for_dls_acquisitions(timeout=acq_time*(num_acqs+5))
            
            # Change measurement names to well names
            set_dls_meas_name(handle, m, sampleNames[m])

            # Saves data for one measurement to data file
            save_dls_exp(handle, save_path)
        
            # ACF function takes up to three consecutive character keys as acceptable settings
            acf_per_meas = pd.DataFrame({key:pd.Series(val,dtype='object') for key, val in get_dls_acf_data(handle, acf_choices, m, -1).items()})
            acf_all_meas = add2frame(acf_per_meas, m, acf_all_meas)

            # For multiple peak regularization settings, a comma-separated list must be provided
            # FUNCTIONALITY NOT DEVELOPED
            regu_per_meas = pd.DataFrame()
            # for i, dt in enumerate(dist_type):
            #     for j, xv in enumerate(regXVal):
            #        if i==0 and j==0:
            if m==0:
                # regu_per_meas = pd.DataFrame(get_dls_regu_data(handle, dt, xv, m, -1))
                regu_per_meas = pd.DataFrame(get_dls_regu_data(handle, dist_type, regXVal, m, -1))
            else:
                # regu_per_meas = pd.concat([regu_per_meas, pd.DataFrame(get_dls_regu_data(handle, dt, xv, m, -1))], ignore_index=True)
                regu_per_meas = pd.concat([regu_per_meas, pd.DataFrame(get_dls_regu_data(handle, dist_type, regXVal, m, -1))], ignore_index=True)
            regu_all_meas = add2frame(regu_per_meas, m, regu_all_meas)

            if m==len(measList)-1:
                old_acf_cols = acf_per_meas.columns
                old_regu_cols = regu_per_meas.columns
            
            if save_per_acq:
                sort_write(acf_per_meas, old_acf_cols, exp_folder, f'{exp_name}_{meas}')
                sort_write(regu_per_meas, old_regu_cols, exp_folder, f'{exp_name}_{meas}')
                
                # Export datalog per measurement to PC
                # NOTE: the Regudist_type saved to the datalog will be (I) even if mass% or number% peaks are collected
                export_dls_datalog(handle, m, -1, export_path_setup(exp_folder, f'{exp_name}_{meas}_Datalog.csv'))
                print(f'Data for measurement {m+1} exported')

        # Saves ACF and peak data
        os.chdir(exp_folder)
        sort_write(acf_all_meas, old_acf_cols, exp_folder, exp_name)
        # for i, dt in enumerate(dist_type):
            # for j, xv in enumerate(regXVal):
        sort_write(regu_all_meas, old_regu_cols, exp_folder, exp_name)
        os.chdir(script_dir)

        # Export overall datalog to PC
        # NOTE: the Regudist_type saved to the datalog will be (I) even if mass% or number% peaks are collected
        export_dls_datalog(handle, -1, -1, export_all_meas)

        # Shut down instrument and experiment
        enable_dls_auto_atten(handle, 0)
        set_dls_laser(handle, 0)

        # DO NOT DECREASE the shutdown time. Otherwise the DLS won't move wells
        stable_time = 30
        print(f'Waiting {stable_time} seconds to shut down')
        await asyncio.sleep(stable_time)

        disconnect_from_dls(handle)
            
        # Closing the DLS experiment will automatically close DYNAMICS software when run via API
        close_dls_exp(handle)