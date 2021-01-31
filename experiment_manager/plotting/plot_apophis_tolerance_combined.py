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

import glob
import os
import sys
import pathlib
sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../results")
sys.path.append(str(pathlib.Path(__file__).parent.absolute()) + "/../experiments")
import experiment_manager as manager
import numpy as np
from matplotlib import pyplot as plt
import initial_conditions as init

# Get latest experiment from results dir
list_of_dirs = glob.glob('../results/*') 
path = max(list_of_dirs, key=os.path.getctime).split('/')[-1]

# Create an experiment manager
manager = manager.ExperimentManager("Apophis") 

# load previously run experiments
manager.load_experiments(path)

# print loaded experiments
manager.print_experiments()

# Within manager there is an array of 'experiment' objects called 'experiments'
# that contains all of the useful data.

# within each experiment object there are two dictionaries 'config' and 'integration'.
# config contains all initial condition data used to create the experiment whereas 
# integration contains all of the output data from the experiment.

# to see what is in each of these for experiment zero print the following:
# manager.experiments[0].config.keys() and manager.experiments[0].integration.keys()
# Important fields here are:
#
#   dH: energy conservation
#   Q: position vector (output in democratic heliocentric coordinates)
#   P: momentum vector (output in democratic heliocentric coordinates)
#   time in orbits: the time in units of orbits (relative to initial orbital period)
#   
# For more information on democratic heliocentric coordinates see the following: DOI 10.1086/300541

# Get the number of orbits run for
orbits = manager.experiments[0].config['orbits']
samples = len(manager.experiments)

tes_de = []
tes_rt = []
tes_fc = []
tes_st = []
tes_pos_apophis = []

for exp in manager.experiments:
    if exp.config['integrator'] == 'tes':
        tes_de.append(np.array(exp.integration['dH'][-1]))
        tes_rt.append(np.array(exp.integration['run time']))
        tes_fc.append(np.array(exp.integration['fCallsPerOrbit']))
        tes_st.append(np.array(exp.integration['steps taken'] / orbits))
        tes_pos_apophis.append(np.array(exp.integration['Q'])[-1, 2, :])

# samples-1 is the highest precision TES run.
Q = np.array(manager.experiments[samples-1].integration['Q'])
t = np.array(manager.experiments[samples-1].integration['time in orbits'])
encounter_idx = np.argmin(np.abs(t-50))

plt.figure(figsize=(4,4))
orbital_plot_marker_size=2
plt.plot(Q[:, 0, 0], Q[:, 0, 1], '.', label='Sun', markersize=orbital_plot_marker_size)
plt.plot(Q[:, 1, 0], Q[:, 1, 1], label='Earth', markersize=orbital_plot_marker_size)
plt.plot(Q[:, 2, 0], Q[:, 2, 1], label='Apophis', markersize=orbital_plot_marker_size)
plt.plot(Q[encounter_idx, 2, 0], Q[encounter_idx, 2, 1], 'X', label='2029 flyby', markersize=10)
plt.legend(loc='upper right')
plt.axis('equal')
plt.xlabel("x-position (AU)")
plt.ylabel("y-position (AU)")
plt.tight_layout()
plt.savefig(manager.results_dir + "apophis_orbital_trace.pdf", markersize=orbital_plot_marker_size)


sep = np.linalg.norm(Q[:, 2, :] - Q[:, 1, :], axis=1)
plt.figure()
plt.semilogy(t+1979, sep)
plt.xlabel('time (years)')
plt.ylabel('separation between Earth and Apophis (AU)')
plt.xlim(1979, 2080)
plt.grid()
plt.savefig(manager.results_dir + "apophis_separation.pdf")

marker_size=5


# Get solution function for comparison against (obtained with quadruple precision integration)
solution = init.GetApophis2079
pos_exact, vel_exact, _,_,_ = solution()

# Convert to heliocentric coordinates
pos_exact_helio = pos_exact - pos_exact[0, :]

tes_pos_err = np.linalg.norm(tes_pos_apophis - pos_exact_helio[-1], axis=1)

plt.figure()
plt.semilogy(tes_st, tes_pos_err, '*', label='TES', markersize=marker_size)
plt.grid()
plt.legend()
plt.xlabel("steps per orbit")
plt.ylabel(r"$\delta r$")
plt.savefig(manager.results_dir + "steps_vs_pos_err.pdf")
