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

class OrbitalElements:
    def __init__(self, r, v, m, t0=0, verbose=False):
        self.FLOATING_POINT_CUT_OFF = 1E-14
        mu = np.sum(m)
        nu = 0
        r_norm = np.linalg.norm(r)
        v_norm = np.linalg.norm(v)
        v_circ = np.sqrt(mu/r_norm)

        a = -mu / (v_norm**2-2.0*v_circ**2)
        a_norm = np.linalg.norm(a)

        h = np.cross(r, v)
        h_norm = np.linalg.norm(h)

        e = (1.0/mu) * ((v_norm**2-(mu/r_norm))*r - np.dot(r, v)*v)
        e_norm = np.linalg.norm(e)

        meanMotion = np.sqrt(mu/a_norm**3)
        period = 2.0*np.pi/meanMotion

        # This currently runs from 0-2*pi, should run from -pi to pi I think.
        i = np.arccos(h[2] / h_norm)

        if(np.isnan(i)):
            raise ValueError('Inclination is NaN.')

        # Unit z vector
        K = np.array([0, 0, 1])

        n = np.cross(K, h)
        n_norm = np.linalg.norm(n)

        energy = 0.5*v_norm**2 - mu / r_norm

        if(i > self.FLOATING_POINT_CUT_OFF):
            Omega = np.arccos(n[0] / n_norm)
        else:
            Omega = 0

        if n[1] < 0:
            Omega = 2*np.pi-Omega

        # Circular equitorial
        if(e_norm < self.FLOATING_POINT_CUT_OFF) and ((i < self.FLOATING_POINT_CUT_OFF) or ((i % np.pi) < self.FLOATING_POINT_CUT_OFF)):
            nu = np.arccos(r[0]/r_norm)
            if(r[1] < 0.0):
                nu = 2*np.pi-nu
            E = nu
            M = nu
            raise ValueError('Circular orbits dont exist.')
            self.mode = 0
        # Circular inclined
        elif(e_norm < self.FLOATING_POINT_CUT_OFF):
            nu = np.arccos(np.dot(n, r) / (n_norm*r_norm))
            if(r[2] < 0):
                nu = 2*np.pi - nu
            E = nu
            raise ValueError('Circular orbits dont exist.')
            self.mode = 1
        elif(e_norm < 1.0):
            nu = np.arccos(np.dot(e, r) / (e_norm*r_norm))
            if(np.dot(r, v) < 0):
                nu = np.pi*2 - nu
            E = np.arccos((e_norm+np.cos(nu))/(1+e_norm*np.cos(nu)))
            if(np.dot(r, v) <= 0):
                E = np.pi*2 - E
            self.mode = 2
        else:
            print("Orbit is either parabolic or hyperbolic")


        self.tempStore = E
        self.Eout = E

        M = E - e_norm*np.sin(E)

        p = h_norm**2/mu

        w = 0
        self.mu = mu
        self.a = a_norm
        self.e = e_norm
        self.i = i
        self.p = p
        self.h = h_norm
        self.r = r_norm
        self.period = period
        self.n = meanMotion
        self.nu = nu
        self.E = E
        self.w = w
        self.M = M
        self.t0 = t0
        self.X = (r, v)
        self.energy = energy
        self.r_vec = r
        self.v_vec = v
        self.Omega = Omega

        self.MAX_NEXTON_ITERATIONS = 30
        self.NEWTON_TOLERANCE = 1E-14   #1E-15 fails sometimes, maybe test for stalling around this value.

#        if(t0 > 0):
        self.PrintElements(verbose)

    def PrintElements(self, verbose=True):
        if(verbose == True):
            print()
            print("r_vec:", self.r_vec)
            print("v_vec:", self.v_vec)
            print("mu:", self.mu)
            print("a:", self.a)
            print("e:", self.e)
            print("i:", self.i)
            print("h:", self.h)
            print("p:", self.p)
            print("M:", self.M)
            print("E:", self.E)
            print("nu:", self.nu)
            print("w:", self.w)
            print("Omega:", self.Omega)
            print("Period:", self.period)
            print("t:", self.t0)
            print("temp:", self.tempStore)
            print("mode:", self.mode)
