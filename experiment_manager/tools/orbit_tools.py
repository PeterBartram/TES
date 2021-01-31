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
from matplotlib import pyplot as plt


def calculate_orbital_elements(r_in, v_in, m):
    period_min = 1E90
    n = len(m)
    
    elems = {}
    elems['a'] = np.zeros(n)
    elems['e'] = np.zeros(n)
    elems['i'] = np.zeros(n)
    elems['mean motion'] = np.zeros(n)
    elems['period'] = np.zeros(n)

    for i in range(1, n):
        r = r_in[i]
        v = v_in[i]
        # Assume that the masses have already been multiplied by G.
        mu = m[0] + m[i]
        r_norm = np.linalg.norm(r)
        v_norm = np.linalg.norm(v)
        v_circ = np.sqrt(mu/r_norm)
    
        a = -mu / (v_norm**2-2.0*v_circ**2)
        a_norm = np.linalg.norm(a)
    
        h = np.cross(r, v)
        h_norm = np.linalg.norm(h)
    
        e = (1.0/mu) * ((v_norm**2-(mu/r_norm))*r - np.dot(r, v)*v)
        e_norm = np.linalg.norm(e)
    
        nu = np.sqrt(mu/a_norm**3)
        period = 2.0*np.pi/nu
    
        inc = np.arccos(h[2] / h_norm)    
    
        elems['a'][i] = a_norm
        elems['e'][i] = e_norm
        elems['i'][i] = inc
        elems['mean motion'][i] = nu
        elems['period'][i] = period

        if(period < period_min):
            period_min = period
    
    elems['minimum period'] = period_min
    
    return elems

