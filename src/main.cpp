// TES is an open source integration package for modelling exoplanet evolution.
// Copyright (C) <2021>  <Peter Bartram, Alexander Wittig>

// TES is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>. 

#include <stdio.h>
#include "Simulation.h"
#include "dhem.h"
#include "radau.h"
#include "radau_step.h"
#include "UniversalVars.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "json.hpp"

static SIMULATION * sim;

using json = nlohmann::json;
using namespace std;

int main(int argc, char * argv[])
{
  // Default input file for if not provided via command line.
  char * file = (char*)"./input.json";

  // Extract input file name from the command line.
  if(argc > 1)
  {
    file = argv[1];
  }

  // Read the input JSON file
  std::ifstream i(file);
  json j;
  i >> j;

  // Get number of bodies.
  double n = j["n"].get<double>();

  // Allocate memory for input arrays.
  double * Q = new double[int(3*n)];
  double * V = new double[int(3*n)];
  double * m = new double[int(n)];

  // Extract all of the input configuration data.
  double t0 = j["t0"].get<double>();
  double period = j["period"].get<double>();
  double orbits = j["orbits"].get<double>();
  uint32_t output_spacing = uint32_t(j["output spacing"].get<double>());
  uint32_t output_samples = uint32_t(j["output samples"].get<double>());
  double rTol = j["rTol"].get<double>();
  double aTol = j["aTol"].get<double>();
  double rectisPerOrbit = j["rectisPerOrbit"].get<double>();
  double dQcutoff = j["dQcutoff"].get<double>();
  double dPcutoff = j["dPcutoff"].get<double>();
  double hInitial = j["hInitial"].get<double>();
  double timeOut = j["timeOut"].get<double>();
  string outputFileString = j["outputFile"].get<string>();
  char * outputFile = (char *)outputFileString.c_str();

  // Extract masses.
  for(uint32_t i = 0; i < n; i++)
  {
    m[i] = j["mass"][to_string(i)].get<double>();
  }

  // Extract positions and velocities.
  for(uint32_t i = 0; i < n; i++)
  {
    Q[3*i+0] = j["particles"][to_string(i)]["x"].get<double>();
    Q[3*i+1] = j["particles"][to_string(i)]["y"].get<double>();
    Q[3*i+2] = j["particles"][to_string(i)]["z"].get<double>();

    V[3*i+0] = j["particles"][to_string(i)]["vx"].get<double>();
    V[3*i+1] = j["particles"][to_string(i)]["vy"].get<double>();
    V[3*i+2] = j["particles"][to_string(i)]["vz"].get<double>();
  }


  sim = Simulation_Init(n);

  // Perform initialisation of our integrator objects.
  Sim_AddInitialConditions(Q, V, m);

  // Store configuration values in the simulation.
  sim->t0 = t0;
  sim->tEnd = period*orbits;
  sim->hInitial = hInitial;
  sim->orbits = orbits;
  sim->aTol = aTol;
  sim->rTol = rTol;
  sim->rectisPerOrbit = rectisPerOrbit;
  sim->orbits = orbits;
  sim->period = period;
  sim->dQcutoff = dQcutoff;
  sim->dPcutoff = dPcutoff;
  sim->outputFile = outputFile;
  sim->timeOut = timeOut;
  sim->output_spacing = output_spacing;
  sim->output_samples = output_samples;
  
  UniversalVars_Init(sim);
  dhem_Init(sim, period/rectisPerOrbit, 9);
  Radau_Init(sim);
  // Perform the integration.
  Radau_integrate();
  // Clean up after onesself.
  UniversalVars_Free();
  dhem_Free();
  Radau_Free();
  Simulation_Free();
}
