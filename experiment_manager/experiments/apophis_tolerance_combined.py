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

# This file can be run using multiple cores simply by using this command from
# the command line in the directory of this file: mpiexec -n 8 python apophis_tolerance_combined.py 
# Replace 8 with the number of cores you so desire.
import os
import sys
import pathlib
sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../")
import numpy as np
from matplotlib import pyplot as plt
import initial_conditions as init
import experiment_manager as Exp

# Initial conditions are at t=1st Jan 1979 and we will integrate for 100 years.
orbits = 100

# Create an experiment manager for running experiments and handling data.
manager = Exp.ExperimentManager("apophis_baseline") 

# Get initial conditions function
problem = init.GetApophis1979

# Number of outputs over the span of our integration (careful as too many points will eat file storage)
output_samples=2500

# We want to vary the tolerance of TES from 10^-1 to 10^-6 and run an integration 
# at each value of tolerance chosen i.e. "samples" integrations
samples = 21
tols_tes = np.logspace(-2, -6, samples)

for tol in tols_tes:
    # Create a config dictionary for parametrising an "experiment".
    config = {}
    
    # specify the number of orbits we want to integrate for (the initial smallest orbital period is 
    # returned by the problem function above and used internally by the driver)
    config['orbits'] = orbits
    
    # configure the tolerance to use in the integration (this is a dimensionless value and 
    # setting it to 10^-6 will not result in an error of 10^-6 but something far smaller)
    config['rtol'] = tol
    
    # name the experiment - free choice here but useful for later identifying experiments.
    config['name'] = 'tes'
    
    # this must be set to 'tes'
    config['integrator'] = 'tes'
    
    # how many outputs do you want over the course of the integration?
    # this is not dense output and will output at the closest step end time.
    config['output samples'] = output_samples
    
    # choices here are 'linear' i.e. output linearly spaced in time
    # or 'log' to output on a base 10 log spacing.
    config['output spacing'] = 'linear'   
    
    # determines how large the deltas can grow before we rectify (just leave these be for now and ensure you
    # specify the orbital period correctly in the problem function)
    config['dQ max'] = 1E-3
    config['dP max'] = 1
    config['rectifications per orbit'] = 1.61803398875
    
    # how small the initial step size should be (1E4 is tiny but it doesnt harm to have it so small.)
    config['steps initial per orbit'] = 1E4
    
    # Set a time out period in seconds - not sure if this is enable atm
    config['time out'] = 2*86400
    
    # the executable file that will be called.
    config['exe file'] = 'output'
    
    # Finally, add this experiment configuration to the experiment manager.
    manager.add_experiment(config, problem)
   
# Run all experiments given to the manager.
manager.run_experiments()

# if this isnt called then results remain in the 'temp' directory, this will 
# transfer them to 'results' dir.
manager.save_experiments()
