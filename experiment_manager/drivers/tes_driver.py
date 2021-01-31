# TES is an open source integration package for modelling exoplanet evolution.
# Copyright (C) <2021>  <Peter Bartram, Alexander Wittig>

# TES is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>. 

import numpy as np
import json
import subprocess
import os
import config_file_generation as config_file_gen
    
def validate_config_input(config):
    required_keys = np.array(['dQ max', 'dP max', 'rectifications per orbit', 'steps initial per orbit', 'time out'])
    keys_present = np.array([x in config.keys() for x in required_keys])
    
    if (keys_present == True).all():        
        return True
    else:
        keys_missing = required_keys[keys_present == False]
        [print("Missing key from config structure: {}".format(x)) for x in keys_missing]    
        raise ValueError()
        
def run(config):
    #  Unused internally but not yet removed.
    atol=1E-50 
    config['output file'] = config['results dir'] + "output.txt"
    config['h0'] = config['period'] / config['steps initial per orbit']
    
    config_file_gen.create_tes_config_file(config)
    exe_path = './../../drivers/tes/'
    
    if 'exe file' in config.keys():
        exe_file = config['exe file']
    else:
        exe_file = 'output'
    res = subprocess.call([exe_path+exe_file, 'input.json'])
    
    return extract_data(config['output file'])
   
def determine_version(output_file):
    with open(output_file) as f:
        data = f.readlines()
        sample = data[-2]
        footer = data[-1]
        footer = footer.split(" ")
        sample_len = len(sample.split(" "))
        n = int(footer[2])
   
        if sample_len == 10+15*n:
            return "v1"
        else:
            return "dev"
    
def extract_data(output_file):        
   # Extract the header to this file.
    output = {}
    with open(output_file) as f:
        footer = f.readlines()[-1]
        params = footer.split(" ")
        output["t0"] = np.float64(params[0])
        output["tEnd"] = np.float64(params[1])
        output["n"] = int(params[2])
        output["orbits"] = np.float64(params[3])
        output["aTol"] = np.float64(params[4])
        output["rTol"] = np.float64(params[5])
        output["period"] = np.float64(params[6])
        output["output interval"] = np.float64(params[7])
        output["step size initial"] = np.float64(params[8])
        output["rectifications per orbit"] = np.float64(params[9])
        output["dQ cutoff"] = np.float64(params[10])
        output["dP cutoff"] = np.float64(params[11])
        output["git revision"] = params[12]
        output["run time"] = np.float64(params[13])
        output["steps taken"] = np.float64(params[14])
    
    period = output["period"]
    n = output["n"]
    orbits = output["orbits"]
    
    data = np.genfromtxt(output_file, delimiter=" ", skip_footer=1, dtype=np.float64)
    
    time = data[:, 0]
    timeInOrbits = data[:, 0]/period
    samples = len(time)
    
    H = data[:, 1]
    dHsign = np.diff((H-H[0])/H[0])
    dHsign[dHsign > 0] = 1
    dHsign[dHsign < 0] = -1
    dH = np.abs((H-H[0])/H[0])
    dH[dH == 0.0] = 1E-16
    dHhp = data[:, 2]
    fCalls = data[:, 3]
    fCallsPerOrbit = fCalls[-1]/orbits
    h = data[:, 4]
    rectifications = data[:, 5]
    recti = data[:, 6]
    iterations = data[:, 7]
    
    state_idx = 8
    # Extract the integration state data.
    Q = data[:, state_idx:state_idx+3*n]
    P = data[:, state_idx+3*n:state_idx+6*n]
    dQ = data[:, state_idx+6*n:state_idx+9*n]
    dP = data[:, state_idx+9*n:state_idx+12*n]    

    # Extract step size control data
    b6Max = data[:, -2]
    accMax = data[:, -1]
    
    Q = np.reshape(Q, (samples, n, 3))
    P = np.reshape(P, (samples, n, 3))
    dQ = np.reshape(dQ, (samples, n, 3))
    dP = np.reshape(dP, (samples, n, 3))

    output["time"] = time 
    output["time in orbits"] = timeInOrbits 
    output["Q"] = Q
    output["P"] = P
    output["dQ"] = dQ
    output["dP"] = dP
    output["Qosc"] = Q-dQ
    output["Posc"] = P-dP
    output["samples"] = samples 
    output["H"] = H 
    output["dH"] = dH 
    output["dHsign"] = dHsign 
    output["fCalls"] = fCalls 
    output["fCallsPerOrbit"] = fCallsPerOrbit 
    output["h"] = h 
    output["rectifications"] = rectifications 
    output["rectified"] = recti
    output["iterations"] = iterations 
    output["b6Max"] = b6Max
    output["accMax"] =  accMax
    
    return output
