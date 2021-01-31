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

import sys
import pathlib
sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../drivers")
sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../tools")
sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../results")
sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../temp")
sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")

try:
    import cPickle
except:
    import pickle

import os
import glob
import numpy as np
import initial_conditions as init
import tes_driver
import time
import json
import MPIServer as mpi_interface
import orbit_tools as orbit

integrators = {}
integrators['tes'] = tes_driver

class Experiment():
    def __init__(self, config=None, problem=None, path=None):
        if path is None:
            if config is None:
                self.config = {}
                print("Experiment is missing config input.")
            else:
                if type(config) == dict:
                    # Check that we have the bare minimum number of configuration fields required for an integration.
                    required_keys = np.array(['orbits', 'rtol', 'name', 'integrator', 'output samples'])
                    keys_present = np.array([x in config.keys() for x in required_keys])
    
                    if (keys_present == True).all():                
                        if config['integrator'] in integrators:
                            config['driver'] = integrators[config['integrator']]
                            
                            # Now call into the driver to check that all integrator specific fields are present.
                            if config['driver'].validate_config_input(config):
                                self.config = config        
                        else:
                            raise ValueError("Unknown integrator selected. Acceptable options are: {}".format(integrators.keys()))
                    else:
                        keys_missing = required_keys[keys_present == False]
                        [print("Missing key from config structure: {}".format(x)) for x in keys_missing]
                        raise ValueError()
                else:
                    raise TypeError("Config input to Experiment class needs to be a dictionary.")
                    
            if problem is None:
                raise ValueError("A valid problem is required for an experiment.")
            else:
                self.problem = problem
                
                # Create output directories
                self.config['here'] = os.path.dirname(os.path.abspath(__file__)) + "/"
                file_name = str(np.random.randint(0, 1E18))
                self.config['file name'] = file_name + ".json"
                self.config['results dir'] = self.config['here'] + '../temp/' + file_name + "/"
                os.makedirs(self.config['results dir'], exist_ok=True)              
                             
                if callable(problem):
                    os.chdir(self.config['results dir'])   
                    self.config['Q0'], self.config['V0'], self.config['mass0'], \
                    self.config['period'], self.config['coords'] = self.problem()   
                    os.chdir(self.config['here'])                      
                elif type(problem) is str:
                    self.config = config
                    self.load_experiment_from_file(config, problem)
                  
        else:
            self.load(path)
        
    def load(self, file):
        # Extract fields from the param file
        with open(file, "r") as f:
            jsobj_str = f.read()

        experiment_data = json.loads(jsobj_str)            
        
        self.config = experiment_data['config']
        self.integration = experiment_data['integration']
           
    def run(self):
        os.chdir(self.config['results dir'])
        t0 = time.time()
        self.integration = self.config['driver'].run(self.config)        
        t1 = time.time()
        self.integration['run time'] = t1-t0
        os.chdir(self.config['here'])             
    
    def save(self, file_name=None):
        data = {}
        data['config'] = self.config.copy()
        del(data['config']['driver'])
        data['integration'] = self.integration.copy()
        
        if file_name == None:
            file_name = self.config['file name']
        file = self.config['final results dir'] + file_name

        # numpy objects need to be converted to lists in order to serialise.
        for field_key in data:
            field = data[field_key]
            for key in field.keys():
                data_type = type(field[key])
                
                if data_type == np.ndarray:
                    field[key] = field[key].tolist()
                    
                if data_type == np.float64:
                    field[key] = float(field[key])
                    
                
        with open(file, "w") as write_file:
            json.dump(data, write_file, skipkeys=True)

    def print_experiment(self):
        print("{} {} {} {}".format(self.config['name'], self.config['integrator'], self.config['rtol'], self.config['orbits']))
        
    
    def load_experiment_from_file(self, config, problem_file):
        f = open(problem_file)
        data = f.read()
        data = json.loads(data)
        
        n = len(data['particles'])
        
        particles = data['particles']
        mass = np.zeros(n)
        pos = np.zeros([n, 3])
        vel = np.zeros([n, 3])
        
        for i in range(n):
            p = particles[str(i)]
            pos[i, 0] = p['x']
            pos[i, 1] = p['y']
            pos[i, 2] = p['z']
            vel[i, 0] = p['vx']
            vel[i, 1] = p['vy']
            vel[i, 2] = p['vz']
            mass[i] = data['mass'][str(i)]        
        
        elems = orbit.calculate_orbital_elements(pos, vel, mass)
        self.config['period'] = elems['minimum period']
        
        self.config['Q0'] = pos
        self.config['V0'] = vel
        self.config['mass0'] = mass
        self.config['coords'] = False
        
class ExperimentManager():
    def __init__(self, name):
        seed = np.int64(time.time())
        # seed = 0
        # Need to seed the same for all MPI threads to ensure output dirs are the same.
        seed = mpi_interface.distribute_seed(seed)
        np.random.seed(seed)
        
        self.single_core = mpi_interface.running_single_core()
        self.experiments = []   
        self.experiment_name = name
        self.run_time = 0

        self.here = os.path.dirname(os.path.abspath(__file__)) + "/"
        file_path = str(np.random.randint(0, 1E18))
        self.results_dir = self.here + '../results/' + self.experiment_name + "_" + file_path + '/'
        
    def add_experiment(self, config, problem):
        # Final location to save outputs to.
        config['final results dir'] = self.results_dir
        self.experiments.append(Experiment(config, problem))
        
    def run_experiments(self):
        # dont create results dir unless we are going to run.
        os.makedirs(self.results_dir, exist_ok=True)
        # Store the temp directories incase we need to find the data later on.
        with open(self.results_dir + 'temp_dirs.txt', 'w') as file:
            for exp in self.experiments:
                file.write("{}\n".format(exp.config['results dir']))
        
        t0 = time.time()
        if self.single_core == True:
            for experiment in self.experiments:
                experiment.run()
        else:
            self.experiments = mpi_interface.run_experiments(self.experiments)
        t1 = time.time()
        self.run_time = t1-t0
            
    def save_experiments(self):        
        if self.single_core == True:
            if not os.path.exists(self.results_dir):
                os.makedirs(self.results_dir)              
            
            for i, experiment in enumerate(self.experiments):
                experiment.save((6-len(str(i)))*"0"+str(i)+'.json')
    
    def load_experiments(self, experiment_dir):
        load_dir = self.here + '../results/' + experiment_dir +'/'
        self.results_dir = load_dir
        
        files = sorted(glob.glob(load_dir + "*.json"))
        for file in files:            
            if not 'input.json' in file:
                self.experiments.append(Experiment(path=file))
            elif(len(files) == 1):
                print('No integration file found')

    def load_experiment_ordered(self, experiment_dir):
        ''' This function is the same as the generic load but instead looks
            at the temp_dirs.txt file as a means of ordering the experiments.'''
        dirs_file = experiment_dir + "/temp_dirs.txt"
        
        with open(dirs_file) as f:
            data = np.array(f.readlines())
        
        files = [experiment_dir + dat.split('/')[-2] +'.json' for dat in data]
        
        for file in files:            
            if not 'input.json' in file:
                self.experiments.append(Experiment(path=file))
            elif(len(files) == 1):
                print('No integration file found')
                        
    def print_experiments(self):
        print("***************************************************")
        print("* Experiments currently loaded:")
        print("***************************************************")
        for experiment in self.experiments:
            experiment.print_experiment()