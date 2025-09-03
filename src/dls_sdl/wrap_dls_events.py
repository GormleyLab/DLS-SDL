import time, asyncio
import threading

from wrap_dls_connect import *
from wrap_dls_helpers import *
from wrap_dls_functions import pos_dls, move_dls, get_dls_temp

#---------ASYNCHRONOUS EVENT TRACKING CLASSES------------

class DLSMeasurementTracker:
    def __init__(self):
        self.event = threading.Event()
        self.lock = threading.Lock()
        self.count = 0
        self.expected = 0

    def reset(self, expected):
        with self.lock:
            self.count = 0
            self.expected = expected
            self.event.clear()

    def increment(self):
        with self.lock:
            self.count += 1
            print(f'Acquisition {self.count}/{self.expected}')
            if self.count >= self.expected:
                self.event.set()

tracker = DLSMeasurementTracker()

class DLSMoveTracker:
    def __init__(self):
        self.event = threading.Event()
        self.lock = threading.Lock()
        self.target_well = None

    def start(self, well):
        with self.lock:
            self.target_well = remove_leading_zero(well)
            self.event.clear()

    def signal_complete(self, well):
        well = remove_leading_zero(well)
        with self.lock:
            if self.target_well == well:
                print(f'Move to {well} completed')
                self.event.set()

    def wait_for_completion(self, well, timeout=30):
        #with self.lock:
            self.event.wait(timeout)            

move_tracker = DLSMoveTracker()

#---------ASYNCHRONOUS EVENT HANDLING------------

# Function that detects whether or not an acquisition has occurred
async def wait_for_dls_acquisitions(timeout=60):
    start = time.time()
    while True:
        with tracker.lock:
            if tracker.count >= tracker.expected:
                print(f'Received {tracker.count}/{tracker.expected} acquisitions')
                return
        if time.time() - start > timeout:
           raise TimeoutError(f'Only {tracker.count}/{tracker.expected} acquisitions received')
        await asyncio.to_thread(tracker.event.wait)

def wait_for_dls_move(handle, well, timeout=20):
    well = remove_leading_zero(well)
    startMove = time.time()
    move_tracker.wait_for_completion(well, timeout)
    if remove_leading_zero(pos_dls(handle)) != remove_leading_zero(well):
        return False
    elapsed = time.time()-startMove
    print('Time elapsed: %.3f seconds' % elapsed)
    return True

async def async_dls_move(handle, well, timeout=20):
    global move_tracker
    move_tracker.start(well)
    move_dls(handle, well)
    successful_move = wait_for_dls_move(handle, well)
    if successful_move:
        await asyncio.sleep(0.5)
        move_tracker.event.clear()
    if not successful_move:
        #raise TimeoutError(f'Failed to move to {well}')
        print(f'Restart move to {well}')
        move_tracker.event.clear()
        move_tracker.start(well)
        move_dls(handle, well)
        wait_for_dls_move(handle, well, timeout)

# Called at the beginning of run to identify whether the fixed temperature has stopped ramping
async def wait_for_stable_dls_temperature(handle, target_temp, tolerance=0.05, timeout=300):
    start_time = time.time()
    while True:
        current_temp, current_rate = get_dls_temp(handle)  
        print(f'Current Temp: {float(current_temp):.2f} °C; Current Rate: {float(current_rate):.3f} °C/s; Target: {float(target_temp):.2f} °C')
        if abs(current_temp - target_temp) <= tolerance:
            print('Temperature is within acceptable range')
            break
        if time.time() - start_time > timeout:
            raise TimeoutError(f'Temperature did not stabilize within {timeout} seconds')
        await asyncio.sleep(5)

#---------DYNAMICS EVENT HANDLING------------

def on_dls_error(handle, msg):
    print(msg)

def on_dls_scheduler_done(handle):
    print('Scheduler done')

def on_dls_event_scheduler_msg(handle):
    print('Event scheduler message received')

def on_dls_acq_done(handle):
    tracker.increment()

def on_dls_meas_done(handle):
    print('Measurement complete')

def on_dls_temp_lock(handle, locked, temp):
    if locked==1:
        print('Temperature locked')
    else:
        print('Temperature unlocked')

def on_dls_auto_atten(handle, start):
    if start==1:
        print('Start enabling')
    else:
        print('Finish enabling')

def on_dls_move_done(handle, wellName):
    move_tracker.signal_complete(wellName)

def on_dls_door_state_change(handle, open):
    if open==1:
        print('Door open')
    else:
        print('Door closed')

# Unused feature from  DYNAMICS 8
def on_dls_door_fully_open(handle, open):
    if open==1:
        print('Door fully opened')
    else:
        print("Door not fully open")

#----------EVENT BINDING-------------
events = [on_dls_error, on_dls_scheduler_done, on_dls_event_scheduler_msg, 
    on_dls_acq_done, on_dls_meas_done, on_dls_temp_lock, 
    on_dls_auto_atten, on_dls_move_done, on_dls_door_state_change]

def bind_events(dynapro, events):
    on_error, on_scheduler_done, on_event_scheduler_msg, on_acq_done, on_meas_done, on_temp_lock, on_auto_atten, on_move_done, on_door_state_change = events
    dynapro.OnDataCollectionError += on_error
    dynapro.OnSchedulerDone += on_scheduler_done
    dynapro.OnEventSchdulerMessage += on_event_scheduler_msg
    dynapro.OnAcquisitionDone += on_acq_done
    dynapro.OnMeasurementDone += on_meas_done
    dynapro.OnTemperatureLock += on_temp_lock
    dynapro.OnAutoAttenuation += on_auto_atten
    dynapro.OnMoveDone += on_move_done
    dynapro.OnDoorStateChange += on_door_state_change
    #dynapro.OnDoorFullyOpen += on_door_fully_open

def unbind_events(dynapro, events):
    on_error, on_scheduler_done, on_event_scheduler_msg, on_acq_done, on_meas_done, on_temp_lock, on_auto_atten, on_move_done, on_door_state_change = events
    dynapro.OnDataCollectionError -= on_error
    dynapro.OnSchedulerDone -= on_scheduler_done
    dynapro.OnEventSchdulerMessage -= on_event_scheduler_msg
    dynapro.OnAcquisitionDone -= on_acq_done
    dynapro.OnMeasurementDone -= on_meas_done
    dynapro.OnTemperatureLock -= on_temp_lock
    dynapro.OnAutoAttenuation -= on_auto_atten
    dynapro.OnMoveDone -= on_move_done
    dynapro.OnDoorStateChange -= on_door_state_change
    #dynapro.OnDoorFullyOpen += on_door_fully_open
