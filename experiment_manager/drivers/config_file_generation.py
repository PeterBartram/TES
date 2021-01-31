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


def create_tes_config_file(config):
    config_file = {}
    n = len(config['mass0'])
    config_file['n'] = n
    config_file['t0'] = 0
    config_file['period'] = config['period']
    config_file['orbits'] = config['orbits']
    config_file['outputPeriod'] = 0 # Remove me!
    config_file['output samples'] = config['output samples']
    config_file['rTol'] = config['rtol']
    config_file['aTol'] = 1E-50                 # Not yet removed from interface.
    config_file['outputFile'] = config['output file']
    config_file['G_AU_KG_DY'] = 1.4881806877180788e-34
    config_file['beta'] = 0
    config_file['timeOut'] = config['time out']
    config_file['hInitial'] = config['h0']
    
    if config['output spacing'].lower() == 'linear':
        config_file['output spacing'] = 0
    elif config['output spacing'].lower() == 'log':
        config_file['output spacing'] = 1
    else:
        raise ValueError('Unknown spacing type')
        
    if 'rectifications per orbit' in config.keys():
        config_file['rectisPerOrbit'] = config['rectifications per orbit']
        config_file['dQcutoff'] = config['dQ max']
        config_file['dPcutoff'] = config['dP max']
        
    
    config_file['mass'] = {}
    for i in range(n):
        config_file['mass']["{}".format(i)] = config['mass0'][i]   
        
        config_file['particles'] = {}
        for i in range(n):
            config_file['particles']["{}".format(i)] = {}
            config_file['particles']["{}".format(i)]['x'] = config['Q0'][i, 0]
            config_file['particles']["{}".format(i)]['y'] = config['Q0'][i, 1]
            config_file['particles']["{}".format(i)]['z'] = config['Q0'][i, 2]
            
            config_file['particles']["{}".format(i)]['vx'] = config['V0'][i, 0]
            config_file['particles']["{}".format(i)]['vy'] = config['V0'][i, 1]
            config_file['particles']["{}".format(i)]['vz'] = config['V0'][i, 2]    
            
    dat = json.dumps(config_file)
    with open('input.json', "w") as f:
        f.write(dat)
    