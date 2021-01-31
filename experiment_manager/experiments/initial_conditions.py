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
import calculate_orbits as Orbit
from matplotlib import pyplot as plt
import json

# This is where you want to set your value of G based upon whatever units you want.
def ConvertToRelativeMasses(mass_real):
    # Gravitational constant and solar mass from MERCURY integration package
    G_au_kg_dy = 1.4881806877180788e-34        
    
    mass_real *= G_au_kg_dy
    return mass_real    


def GetApophis1979():
    # Units of AU/DY
    #  These are the values obtained using LSODA in QP
    pos = np.array([
[ -3.66865553481315003192935220049211015E-0004,  8.68061903348902664675614983305406723E-0004, -6.57480629496675813387502045597680369E-0008],
[-0.681618689307206140812555007278134184, -0.743963925714684822346465722102937274, 5.54093722932563831906784519880195250E-0005],
[0.754400432555192513256059016501584560, -1.03275478519758113802166662683288968E-0003, 1.86677482779400751562592845902460873E-0002]])
    
    vel = np.array([
[-1.72400670703383871022118076202455437E-0008,-1.23856827171950553401116128699588587E-0008, 9.16910958998824043960198896327954771E-0013],
[1.24156263544141578553819498367599528E-002, -1.16838072791882536519626721298112578E-002,  8.91997989641411943317485660243687215E-007],
[1.72235066212605091797866733227687193E-003,  2.13892936850632403787793547765985795E-002, -1.08924811794808364682434351962720933E-003]])
    
    mass = np.array([
            1.9884158605722262e30, 
            5.97219E24, 
            2.699E10])
    
    pos, vel, mass, period, coords = CreateInitialConditions(pos, vel, mass)
    # Perform experiment in terms of Earth's orbital period.
    period = 365.24141698304345
    
    return pos, vel, mass, period, coords


def GetApophis2079():
    # Units of AU/DY
    pos = np.array([
[3.66392626256597328138763138508417811E-0004,-8.68270049698502360506398624832847831E-0004,6.57632992286878354588095748176690125E-0008],
[-0.995681074350545687452403150836864168      ,2.66292151344132726551329333265510313E-0003,-1.13965475735725477564900294140012876E-0006],
[-0.810532860073168858369775754232700806      ,-0.541815796359805393814609645989955287      ,6.89729176791859801316165549570176858E-0003]])
    
    vel = np.array([
[2.10563647378310193628812910329311122E-0008,4.40520832131313163436815379489439748E-0009,-3.12041832390222214547945535858613987E-0013],
[-3.35011587561404726138921941339289811E-0004,-1.72742647251350157091438509181831871E-0002,1.30117271821582419179884596154853254E-0006],
[1.23953631694646262658876701302565352E-0002,-1.36089724769463432106672919948336350E-0002,6.76995314038052418791457820632721928E-0004]])
    
    mass = np.array([
            1.9884158605722262e30, 
            5.97219E24, 
            2.699E10])
    
    pos, vel, mass, period, coords = CreateInitialConditions(pos, vel, mass)
    # Perform experiment in terms of Earth's orbital period.
    period = 365.24141698304345
    
    return pos, vel, mass, period, coords


def GetApophisAtFlyby():
    # Units of AU/DY
    #  These are the values taken from the Horizons system that are used to
    #  integrate backwards to get a reference trajectory.
    pos = np.array([
    [0.0, 0.0, 0.0],
    [-9.173794596324105E-01, -4.053012334792483E-01, 2.967128668305804E-05],
    [-9.175062761494387E-01, -4.050870189551113E-01, 6.995085009543227E-05]])
    
    vel = np.array([
    [0.0, 0.0, 0.0],
    [6.675617613903552E-03, -1.580756890748290E-02, 1.197279685894811E-06],
    [1.034166119416305E-02, -1.384354206066643E-02, 1.065794929006502E-03]])
    
    mass = np.array([
            1.9884158605722262e30, 
            5.97219E24, 
            2.699E10])
    
    pos, vel, mass, period, coords = CreateInitialConditions(pos, vel, mass)
    # Perform experiment in terms of Earth's orbital period.
    period = 365.24141698304345
    
    return pos, vel, mass, period, coords


# Call this on your initial conditions to convert to a fixed COM and get the orbital period etc.
def CreateInitialConditions(pos, vel, mass_real):
    period = 1E90

    mass = ConvertToRelativeMasses(mass_real)
    
    pos, vel, mass_real = ConvertToFixedFrames(pos, vel, mass_real)
    
    for i in range(1, len(mass)):
        orbit = Orbit.OrbitalElements(pos[i], vel[i], (mass[0]+mass[i]), verbose=False)

        if(orbit.period < period):
            period = orbit.period
                
    return pos, vel, mass, period, False


def ConvertToFixedFrames(pos, vel, mass):
    q = np.copy(pos)
    v = np.copy(vel)
    m = np.copy(mass)
    mFull = np.repeat(m, 3)
    mFull.shape = (len(m), 3)
    
    mom = mFull * v
    momTot = np.sum(mom, axis=0)
    vTot = momTot / np.sum(m)
    v -= vTot
    
    com = np.sum(q*mFull, axis=0) / np.sum(m)
    q -= com
    
    return q, v, mass





