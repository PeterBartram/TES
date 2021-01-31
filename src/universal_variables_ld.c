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

#include "UniversalVars.h"
#include "radau.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "dhem.h"

#define PI (long double)3.141592653589793238462643383279L


long double invfactorial[35] = {1.0L,1.0L,0.5L,0.166666666666666666666666666666667L,0.0416666666666666666666666666666667L,
0.00833333333333333333333333333333333L,0.00138888888888888888888888888888889L,0.000198412698412698412698412698412698L,
0.0000248015873015873015873015873015873L,0.00000275573192239858906525573192239859L,0.000000275573192239858906525573192239859L,
0.0000000250521083854417187750521083854417L,0.00000000208767569878680989792100903212014L,0.000000000160590438368216145993923771701549L,
1.14707455977297247138516979786821e-11L,7.64716373181981647590113198578807e-13L,4.77947733238738529743820749111754e-14L,
2.81145725434552076319894558301032e-15L,1.56192069685862264622163643500573e-16L,8.22063524662432971695598123687228e-18L,
4.11031762331216485847799061843614e-19L,1.95729410633912612308475743735054e-20L,8.89679139245057328674889744250247e-22L,
3.86817017063068403771691193152281e-23L,1.61173757109611834904871330480117e-24L,6.44695028438447339619485321920469e-26L,
2.47959626322479746007494354584796e-27L,9.18368986379554614842571683647391e-29L,3.27988923706983791015204172731211e-30L,
1.13099628864477169315587645769383e-31L,3.76998762881590564385292152564611e-33L,1.21612504155351794962997468569229e-34L,
3.80039075485474359259367089278841e-36L,1.15163356207719502805868814932982e-37L,3.38715753552116184723143573332301e-39L};

#define NORM(x, i) sqrt2(x[3*i]*x[3*i]+x[3*i+1]*x[3*i+1]+x[3*i+2]*x[3*i+2])
#define DOT(x, y, i) x[3*i+0]*y[3*i+0]+x[3*i+1]*y[3*i+1]+x[3*i+2]*y[3*i+2]

#define CROSS(a,b,c,i) \
	(a)[3*i+0] = (b)[3*i+1] * (c)[3*i+2] - (c)[3*i+1] * (b)[3*i+2]; \
	(a)[3*i+1] = (b)[3*i+2] * (c)[3*i+0] - (c)[3*i+2] * (b)[3*i+0]; \
	(a)[3*i+2] = (b)[3*i+0] * (c)[3*i+1] - (c)[3*i+0] * (b)[3*i+1];

#define CROSS_LOCAL(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

#define MAX_NEWTON_ITERATIONS 50

static void RebasisOsculatingOrbitsMomenta_Internal(long double * z_Q, long double * z_V, double z_t, uint32_t i);
static void RebasisOsculatingOrbitsCore(double z_t, uint32_t i);
static void CalculateClassicalOrbitalElementsSingle(uint32_t i);
static void C_Stumpff(long double * cs, long double z);
static long double CalcUpdateValue(long double X, long double dt, uint32_t i, long double * C);
static long double SolveForUnivsersalAnomaly(long double dt, long double h, uint32_t i, long double * C);
static inline void add_cs(long double* out, long double* cs, long double inp);
static long double sqrt2(long double x);
static long double fabs2(long double x);
static long double fmod2(long double x, long double div);

static UNIVERSAL_VARS * p_uVars = NULL;
static SIMULATION * sim = NULL;


void CalculateOsculatingOrbitsForSingleStep(double **Xosc_map, 
                                            const double t0, const double h, double const * const h_array, 
                                            const uint32_t z_stagesPerStep, uint32_t z_rebasis)
{
  double * Qout;
  double * Pout;

  for(uint32_t stage = 0; stage < z_stagesPerStep; stage++)
  {
    volatile uint32_t temp = 1;
    double t = t0 + h*h_array[stage];

    for(uint32_t i = 1; i < sim->n; i++)
    {
      double * Qout = Xosc_map[stage];
      double * Pout = &Qout[3*sim->n];
   
      long double C[4] = {0.0, 0.0, 0.0, 0.0};

      // Calculate the dt value and wrap around the orbital period.
      long double dt = 0;

      // This will not allow us to maintain the osculating orbit for more than 
      // a step as we would need a dt = t_rebasis + h*(h_arry[] ...) version.
      if(stage == 0)
      {
        dt = 0;
      }
      else
      {
        dt = h*(h_array[stage]-h_array[stage-1]); 
      }

      dt = fmod2(dt, p_uVars->period[i]);
      
      // Calculate our step since last time we were called and update storage of tLast.
      long double h = (long double)t - (long double)p_uVars->tLast[i];
      p_uVars->tLast[i] = t;

      long double X = SolveForUnivsersalAnomaly(dt, h, i, C);

      long double G[4] = {0.0, 0.0, 0.0, 0.0};
      G[0] = C[0]*X;
      G[1] = C[1]*X;
      G[2] = C[2]*(X*X);
      G[3] = C[3]*(X*X*X);

      long double f = -p_uVars->mu*G[2] / p_uVars->Q0_norm[i];
      long double g = dt - p_uVars->mu*G[3];

      for(uint32_t j = 0; j < 3; j++)
      {
        p_uVars->Q1[3*i+j] = (f*p_uVars->Q0[i*3+j] +g*p_uVars->V0[i*3+j]) + p_uVars->Q0[i*3+j];
        Qout[i*3+j] = (double)p_uVars->Q1[3*i+j];
      }

      long double Qnorm = NORM(p_uVars->Q1, i);
      long double fp = -p_uVars->mu * G[1] / (p_uVars->Q0_norm[i]*Qnorm);
      long double gp =  -p_uVars->mu * G[2] / Qnorm;

      for(uint32_t j = 0; j < 3; j++)
      {
        p_uVars->P1[3*i+j] = ((fp*p_uVars->Q0[i*3+j] + gp*p_uVars->V0[i*3+j]) + p_uVars->V0[i*3+j])*p_uVars->mass[i];
        Pout[i*3+j] = (double)p_uVars->P1[3*i+j];  
      }

      if(stage+1 < z_stagesPerStep)
      {
        RebasisOsculatingOrbitsMomenta_Internal(p_uVars->Q1, p_uVars->P1, t, i);    
      }
    }
  } 
}

void ApplyCorrectorToOsculatingOrbitCalculation(double **Xosc_map, double t, uint32_t z_stagePerStep)
{ 
    double * Qout = Xosc_map[z_stagePerStep-1];
    double * Pout = &Qout[3*sim->n];      

    for(uint32_t i = 1; i < sim->n; i++)
    {
        for(uint32_t j = 0; j < 3; j++)
        {
          long double cs_q = sim->radau->cs_dq[3*i+j];
          long double cs_p = sim->radau->cs_dp[3*i+j];
          p_uVars->Q1[3*i+j] -= cs_q;
          p_uVars->P1[3*i+j] -= cs_p;
          p_uVars->Q1[3*i+j] -= p_uVars->uv_csq[3*i+j];
          p_uVars->P1[3*i+j] -= p_uVars->uv_csp[3*i+j];            

          Qout[i*3+j] = (double)p_uVars->Q1[3*i+j];
          Pout[i*3+j] = (double)p_uVars->P1[3*i+j];

          p_uVars->uv_csq[3*i+j] = (long double)Qout[3*i+j] - p_uVars->Q1[3*i+j];
          p_uVars->uv_csp[3*i+j] = (long double)Pout[3*i+j] - p_uVars->P1[3*i+j];

          // In this second operation we then take anything left over from the oscualting orbit corrector
          // and place this into the compensated sum for the delta.
          long double dq = (long double)sim->radau->dQ[3*i+j];
          long double dp = (long double)sim->radau->dP[3*i+j];
          dq -= p_uVars->uv_csq[3*i+j];
          dp -= p_uVars->uv_csp[3*i+j];

          sim->radau->dQ[3*i+j] = (double)dq;
          sim->radau->dP[3*i+j] = (double)dp;

          sim->radau->cs_dq[3*i+j] = (long double)sim->radau->dQ[3*i+j] - dq;
          sim->radau->cs_dp[3*i+j] = (long double)sim->radau->dP[3*i+j] - dp;

          sim->radau->cs_dq[3*i+j] = 0;
          sim->radau->cs_dp[3*i+j] = 0;

          p_uVars->uv_csq[3*i+j] = 0;
          p_uVars->uv_csp[3*i+j] = 0;
        }
        RebasisOsculatingOrbits_Momenta(Qout, Pout, t, i);  
    }
  CalculateClassicalOrbitalElements(); // Remove this (All classical elements can be removed actually)
}


void UniversalVars_Init(SIMULATION * z_sim)
{
  sim = z_sim;
  // Create the main control data structure for universal variables.
  p_uVars = (UNIVERSAL_VARS *)malloc(sizeof(UNIVERSAL_VARS));
  memset(p_uVars, 0, sizeof(UNIVERSAL_VARS));

  z_sim->uVars = p_uVars;

  p_uVars->stateVectorSize = 3 * z_sim->n * sizeof(long double);
  p_uVars->controlVectorSize = z_sim->n * sizeof(long double);

  // Allocate memory for all objects used in universal vars
  p_uVars->Q0 = (long double *)malloc(p_uVars->stateVectorSize);
  p_uVars->V0 = (long double *)malloc(p_uVars->stateVectorSize);
  p_uVars->Q1 = (long double *)malloc(p_uVars->stateVectorSize);
  p_uVars->P1 = (long double *)malloc(p_uVars->stateVectorSize);
  p_uVars->t0 = (double *)malloc(z_sim->n*sizeof(double));
  p_uVars->tLast = (double *)malloc(z_sim->n*sizeof(double));
  p_uVars->Q0_norm = (long double *)malloc(p_uVars->controlVectorSize);
  p_uVars->beta = (long double *)malloc(p_uVars->controlVectorSize);
  p_uVars->eta = (long double *)malloc(p_uVars->controlVectorSize);
  p_uVars->zeta = (long double *)malloc(p_uVars->controlVectorSize);
  p_uVars->period = (long double *)malloc(p_uVars->controlVectorSize);
  p_uVars->Xperiod = (long double *)malloc(p_uVars->controlVectorSize);
  p_uVars->X = (long double *)malloc(p_uVars->controlVectorSize);
  p_uVars->mu = (long double)sim->G*(long double)sim->mass[0];
  p_uVars->e = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->a = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->h = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->h_norm = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->peri = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->apo = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->mass = (long double *)malloc(p_uVars->controlVectorSize);
  p_uVars->dt = (long double *)malloc(p_uVars->controlVectorSize);

  // Take a local copy of the mass vector to save casting.
  for(uint32_t i = 0; i < z_sim->n; i++)
  {
    p_uVars->mass[i] = (long double)z_sim->mass[i];
  }

  memset(p_uVars->Q0, 0, p_uVars->stateVectorSize);
  memset(p_uVars->V0, 0, p_uVars->stateVectorSize);
  memset(p_uVars->Q1, 0, p_uVars->stateVectorSize);
  memset(p_uVars->P1, 0, p_uVars->stateVectorSize);
  memset(p_uVars->t0, 0, z_sim->n*sizeof(double));
  memset(p_uVars->tLast, 0, z_sim->n*sizeof(double));
  memset(p_uVars->Q0_norm, 0, p_uVars->controlVectorSize);
  memset(p_uVars->beta, 0, p_uVars->controlVectorSize);
  memset(p_uVars->eta, 0, p_uVars->controlVectorSize);
  memset(p_uVars->zeta, 0, p_uVars->controlVectorSize);
  memset(p_uVars->period, 0, p_uVars->controlVectorSize);
  memset(p_uVars->Xperiod, 0, p_uVars->controlVectorSize);
  memset(p_uVars->X, 0, p_uVars->controlVectorSize);
  memset(p_uVars->a, 0, p_uVars->controlVectorSize);
  memset(p_uVars->e, 0, p_uVars->controlVectorSize);
  memset(p_uVars->h, 0, p_uVars->stateVectorSize);
  memset(p_uVars->h_norm, 0, p_uVars->controlVectorSize);
  memset(p_uVars->peri, 0, p_uVars->controlVectorSize);
  memset(p_uVars->apo, 0, p_uVars->controlVectorSize);
  memset(p_uVars->dt, 0, p_uVars->controlVectorSize);

  // CS vars
  p_uVars->uv_csq = (double *)malloc(3 * z_sim->n * sizeof(double));
  p_uVars->uv_csv = (double *)malloc(3 * z_sim->n * sizeof(double));
  p_uVars->uv_csp = (double *)malloc(3 * z_sim->n * sizeof(double));
  memset(p_uVars->uv_csq, 0, 3 * z_sim->n * sizeof(double));
  memset(p_uVars->uv_csv, 0, 3 * z_sim->n * sizeof(double));
  memset(p_uVars->uv_csp, 0, 3 * z_sim->n * sizeof(double));  
}

void RebasisOsculatingOrbits_Momenta(double * z_Q, double * z_P, double z_t, uint32_t i)
{
  p_uVars->Q0[3*i+0] = z_Q[3*i+0];
  p_uVars->Q0[3*i+1] = z_Q[3*i+1];
  p_uVars->Q0[3*i+2] = z_Q[3*i+2];
  p_uVars->V0[3*i+0] = z_P[3*i+0] / p_uVars->mass[i];
  p_uVars->V0[3*i+1] = z_P[3*i+1] / p_uVars->mass[i];
  p_uVars->V0[3*i+2] = z_P[3*i+2] / p_uVars->mass[i];

  RebasisOsculatingOrbitsCore(z_t, i);
  CalculateClassicalOrbitalElementsSingle(i);
}


void RebasisOsculatingOrbits(double * z_Q, double * z_V, double z_t, uint32_t i)
{
  p_uVars->Q0[3*i+0] = z_Q[3*i+0];
  p_uVars->Q0[3*i+1] = z_Q[3*i+1];
  p_uVars->Q0[3*i+2] = z_Q[3*i+2];
  p_uVars->V0[3*i+0] = z_V[3*i+0];
  p_uVars->V0[3*i+1] = z_V[3*i+1];
  p_uVars->V0[3*i+2] = z_V[3*i+2];

  RebasisOsculatingOrbitsCore(z_t, i);
  CalculateClassicalOrbitalElementsSingle(i);
}

void RebasisOsculatingOrbitsMomenta_Internal(long double * z_Q, long double * z_P, double z_t, uint32_t i)
{
  p_uVars->Q0[3*i+0] = z_Q[3*i+0];
  p_uVars->Q0[3*i+1] = z_Q[3*i+1];
  p_uVars->Q0[3*i+2] = z_Q[3*i+2];
  p_uVars->V0[3*i+0] = z_P[3*i+0] / p_uVars->mass[i];
  p_uVars->V0[3*i+1] = z_P[3*i+1] / p_uVars->mass[i];
  p_uVars->V0[3*i+2] = z_P[3*i+2] / p_uVars->mass[i];

  RebasisOsculatingOrbitsCore(z_t, i);
  CalculateClassicalOrbitalElementsSingle(i);
}

static void RebasisOsculatingOrbitsCore(double z_t, uint32_t i)
{
  p_uVars->t0[i] = z_t;
  p_uVars->Q0_norm[i] = NORM(p_uVars->Q0, i);
  long double V0_norm2 = NORM(p_uVars->V0, i);
  V0_norm2 *= V0_norm2;

  p_uVars->beta[i] = ((long double)2.0L*p_uVars->mu / p_uVars->Q0_norm[i]) - V0_norm2;
  p_uVars->zeta[i] = p_uVars->mu - p_uVars->beta[i] * p_uVars->Q0_norm[i];
  p_uVars->eta[i] = DOT(p_uVars->Q0, p_uVars->V0, i);
  long double temp = p_uVars->mu / p_uVars->beta[i];
  p_uVars->period[i] = ((long double)4.0L*PI*PI / p_uVars->mu)*(temp*temp*temp);
  p_uVars->period[i] = sqrt2(p_uVars->period[i]);
  p_uVars->Xperiod[i] = (long double)2.0L*PI / sqrt2(p_uVars->beta[i]);
  p_uVars->X[i] = 0.0;
  p_uVars->dt[i] = 0;
}

void GetStateVectorAfterTime(double z_t, uint32_t i, double * Qout, double * Pout)
{
  long double C[4] = {0.0, 0.0, 0.0, 0.0};
  // Calculate the dt value and wrap around the orbital period.
   long double dt = (long double)z_t - (long double)p_uVars->t0[i];
   dt = fmod2(dt, p_uVars->period[i]); 

   // Calculate our step since last time we were called and update storage of tast.
   long double h = (long double)z_t - (long double)p_uVars->tLast[i];
   p_uVars->tLast[i] = z_t;

   long double X = SolveForUnivsersalAnomaly(dt, h, i, C);

   long double G[4] = {0.0, 0.0, 0.0, 0.0};
   G[0] = C[0]*X;
   G[1] = C[1]*X;
   G[2] = C[2]*(X*X);
   G[3] = C[3]*(X*X*X);

   long double f = -p_uVars->mu*G[2] / p_uVars->Q0_norm[i];
   long double g = dt - p_uVars->mu*G[3];

   for(uint32_t j = 0; j < 3; j++)
   {
      p_uVars->Q1[3*i+j] = (f*p_uVars->Q0[i*3+j] +g*p_uVars->V0[i*3+j]) + p_uVars->Q0[i*3+j];
      Qout[i*3+j] = (double)p_uVars->Q1[3*i+j];

      // Assuming our values of Q1 are good to either long or quad precision.
      // Qout_cs[3*i+j] = p_uVars->Q1[3*i+j] - (long double)Qout[3*i+j];
   }

   long double Qnorm = NORM(p_uVars->Q1, i);
   long double fp = -p_uVars->mu * G[1] / (p_uVars->Q0_norm[i]*Qnorm);
   long double gp =  -p_uVars->mu * G[2] / Qnorm;

   for(uint32_t j = 0; j < 3; j++)
   {
      p_uVars->P1[3*i+j] = ((fp*p_uVars->Q0[i*3+j] + gp*p_uVars->V0[i*3+j]) + p_uVars->V0[i*3+j])*p_uVars->mass[i];
      Pout[i*3+j] = (double)p_uVars->P1[3*i+j];

      // Assuming our values of P1 are good to either long or quad precision.
      // Pout_cs[3*i+j] = p_uVars->P1[3*i+j] - (long double)Pout[3*i+j];
   }


  RebasisOsculatingOrbitsMomenta_Internal(p_uVars->Q1, p_uVars->P1, z_t, i);
  // RebasisOsculatingOrbits_Momenta(Qout, Pout, z_t, i);
}


static long double SolveForUnivsersalAnomaly(long double dt, long double h, uint32_t i, long double * C)
{
  // First order guess as per Rein but should try the 2nd order one as well!
  long double X = p_uVars->X[i] + (h / p_uVars->Q0_norm[i]);
  X = fmod2(X, p_uVars->Xperiod[i]);

  long double prevVals[MAX_NEWTON_ITERATIONS];
  for(uint32_t k = 0; k < MAX_NEWTON_ITERATIONS; k++)
  {
    prevVals[i] = 1E37;
  }

  for(uint32_t k = 0; k < MAX_NEWTON_ITERATIONS; k++)
  {
    X = CalcUpdateValue(X, dt, i, C);

    for(uint32_t j = 0; j < k; j++)
    {
      if(X == prevVals[j])
      {
        p_uVars->X[i] = X;
        return X;
      }
    }
    prevVals[k] = X;
  }
  return -(long double)1.0L;
}

static long double CalcUpdateValue(long double X, long double dt, uint32_t i, long double * C)
{
  long double X2 = X*X;
  long double X3 = X2*X;
  long double Z = p_uVars->beta[i]*X2;
  C_Stumpff(C, Z);

  long double num = (X*(p_uVars->eta[i]*X*C[1] + p_uVars->zeta[i]*X2*C[2]) -
              p_uVars->eta[i]*X2*C[2] -
              p_uVars->zeta[i]*X3*C[3] + dt);

  long double denom = p_uVars->Q0_norm[i] + p_uVars->eta[i]*X*C[1] + p_uVars->zeta[i]*X2*C[2];

  return (num / denom);
}

/**
Note: this function is taken from the rebound code base.
*/
static void C_Stumpff(long double * cs, long double z)
{
    unsigned int n = 0;
    while(fabs2(z)>0.1){
        z = z/(long double)4.0L;
        n++;
    }
    
    uint32_t nmax = 0; 
    nmax = 13;    

     // This must be an odd number. 13 is correct for double precision and 21 is correct for quadruple.
    long double c_odd  = invfactorial[nmax];
    long double c_even = invfactorial[nmax-1];
    for(int np=nmax-2;np>=3;np-=2){
        c_odd  = invfactorial[np]    - z *c_odd;
        c_even = invfactorial[np-1]  - z *c_even;
    }

    cs[3] = c_odd;
    cs[2] = c_even;
    cs[1] = invfactorial[1]  - z *c_odd;
    cs[0] = invfactorial[0]  - z *c_even;

    for (;n>0;n--){
        cs[3] = (cs[2]+cs[0]*cs[3])*(long double)0.25L;
        cs[2] = cs[1]*cs[1]*(long double)0.5L;
        cs[1] = cs[0]*cs[1];
        cs[0] = (long double)2.0L*cs[0]*cs[0]-(long double)1.0L;
    }
}

static void CalculateClassicalOrbitalElementsSingle(uint32_t i)
{
  if(sim->termination_check_enable)
  {
    // Semimajor axis
    p_uVars->a[i] = p_uVars->mu / p_uVars->beta[i];

    // Angular momentum
    CROSS(p_uVars->h, p_uVars->Q0, p_uVars->V0, i);
    p_uVars->h_norm[i] = NORM(p_uVars->h, i);

    // Eccentricity
    long double eccVec[3] = {0,0,0};
    CROSS_LOCAL(eccVec, &p_uVars->V0[3*i], &p_uVars->h[3*i]);
    eccVec[0] /= p_uVars->mu;
    eccVec[1] /= p_uVars->mu;
    eccVec[2] /= p_uVars->mu;
    long double Q_norm = NORM(p_uVars->Q0, i);
    eccVec[0] -= p_uVars->Q0[3*i+0] / Q_norm;
    eccVec[1] -= p_uVars->Q0[3*i+1] / Q_norm;
    eccVec[2] -= p_uVars->Q0[3*i+2] / Q_norm;
    p_uVars->e[i] = NORM(eccVec, 0);

    // Apogee and perigee.
    p_uVars->peri[i] = p_uVars->a[i]*(1-p_uVars->e[i]);
    p_uVars->apo[i] = p_uVars->a[i]*(1+p_uVars->e[i]);
  }
}

void CalculateClassicalOrbitalElements(void)
{
  if(sim->termination_check_enable)
  {
    for(uint32_t i = 1; i < sim->n; i++)
    {
      CalculateClassicalOrbitalElementsSingle(i);
    }
  }


  // printf("\n\ne: ");
  // for(uint32_t i = 1; i < sim->n; i++)
  // {
  //   printf(" %.5E", (double)p_uVars->e[i]);
  // }

  // printf("\na: ");
  // for(uint32_t i = 1; i < sim->n; i++)
  // {
  //   printf(" %.5E", (double)p_uVars->a[i]);
  // }


}

void UniversalVars_Free(void)
{
  free(p_uVars->Q0);
  free(p_uVars->V0);
  free(p_uVars->Q1);
  free(p_uVars->P1);
  free(p_uVars->t0);
  free(p_uVars->tLast);
  free(p_uVars->Q0_norm);
  free(p_uVars->beta);
  free(p_uVars->eta);
  free(p_uVars->zeta);
  free(p_uVars->period);
  free(p_uVars->Xperiod);
  free(p_uVars->X);
  free(p_uVars);
  sim->uVars = NULL;

}


static long double sqrt2(long double x)
{
  #if TYPE_INTERNAL == TYPE_DOUBLE
    return sqrt(x);
  #elif TYPE_INTERNAL == TYPE_LONG_DOUBLE
    return sqrtl(x);
  #elif TYPE_INTERNAL == TYPE_FLOAT128
    return sqrtq(x);
  #endif
}


static long double fmod2(long double x, long double div)
{
  #if TYPE_INTERNAL == TYPE_DOUBLE
    return fmod(x, div);
  #elif TYPE_INTERNAL == TYPE_LONG_DOUBLE
    return fmodl(x, div);
  #elif TYPE_INTERNAL == TYPE_FLOAT128
    return fmodq(x, div);
  #endif
}

static long double fabs2(long double x)
{
  #if TYPE_INTERNAL == TYPE_DOUBLE
    return fabs(x);
  #elif TYPE_INTERNAL == TYPE_LONG_DOUBLE
    return fabsl(x);
  #elif TYPE_INTERNAL == TYPE_FLOAT128
    return fabsq(x);
  #endif
}


static inline void add_cs(long double * out, long double * cs, long double inp)
{
    // if(sim->enable_cs == 0)
    if(1)
    {
      out[0] += inp;
    }
    else
    {
      const long double y = inp - cs[0];
      const long double t = out[0] + y;
      cs[0] = (t - out[0]) - y;
      out[0] = t;      
    }
}
