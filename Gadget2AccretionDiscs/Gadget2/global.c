#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file global.c
 *  \brief Computes global physical properties of the system
 */


/*! This routine computes various global properties of the particle
 *  distribution and stores the result in the struct `SysState'.
 *  Currently, not all the information that's computed here is actually
 *  used (e.g. momentum is not really used anywhere), just the energies are
 *  written to a log-file every once in a while.
 */
void compute_global_quantities_of_system(void)
{
  int i, j, n;
  struct state_of_system sys;
  double a1, a2, a3;
  double entr = 0, egyspec, vel[3];
  double dt_entr, dt_gravkick, dt_hydrokick;


  //We want to add the extra little bits needed to output the radiated energy, total angular/linear momentum...

  if(All.ComovingIntegrationOn)
    {
      a1 = All.Time;
      a2 = All.Time * All.Time;
      a3 = All.Time * All.Time * All.Time;
    }
  else
    {
      a1 = a2 = a3 = 1;
    }


  for(n = 0; n < 6; n++)
    {
      sys.MassComp[n] = sys.EnergyKinComp[n] = sys.EnergyPotComp[n] = sys.EnergyIntComp[n] = 0;
#ifdef EXTRA_STATS
      sys.EnergyRadComp[n]=0;
#endif

      for(j = 0; j < 4; j++)
	sys.CenterOfMassComp[n][j] = sys.MomentumComp[n][j] = sys.AngMomentumComp[n][j] = 0;
    }

  for(i = 0; i < NumPart; i++)
    {
      sys.MassComp[P[i].Type] += P[i].Mass;

      //Potential is only stored in the orbiting particles if we have self gravity turned off, so we don't need the factor of a half
#ifdef ADD_CENTRAL_GRAVITY
      sys.EnergyPotComp[P[i].Type] += P[i].Mass * P[i].Potential / a1;
#else
      sys.EnergyPotComp[P[i].Type] += 0.5 * P[i].Mass * P[i].Potential / a1;
#endif

      if(All.ComovingIntegrationOn)
	{
	  dt_entr = (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;
	  dt_gravkick = get_gravkick_factor(P[i].Ti_begstep, All.Ti_Current) -
	    get_gravkick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
	  dt_hydrokick = get_hydrokick_factor(P[i].Ti_begstep, All.Ti_Current) -
	    get_hydrokick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
	}
      else
	dt_entr = dt_gravkick = dt_hydrokick =
	  (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

      for(j = 0; j < 3; j++)
	{
	  vel[j] = P[i].Vel[j] + P[i].GravAccel[j] * dt_gravkick;
	  if(P[i].Type == 0)
	    vel[j] += SphP[i].HydroAccel[j] * dt_hydrokick;
	}
      if(P[i].Type == 0)
	entr = SphP[i].Entropy + SphP[i].DtEntropy * dt_entr;

#ifdef PMGRID
      if(All.ComovingIntegrationOn)
	dt_gravkick = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
	  get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
      else
	dt_gravkick = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;

      for(j = 0; j < 3; j++)
	vel[j] += P[i].GravPM[j] * dt_gravkick;
#endif

      sys.EnergyKinComp[P[i].Type] +=
	0.5 * P[i].Mass * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) / a2;

      if(P[i].Type == 0)
	{
#ifdef ISOTHERM_EQS
          egyspec = entr;
#else
	  egyspec = entr / (GAMMA_MINUS1) * pow(SphP[i].Density / a3, GAMMA_MINUS1);
#endif
	  sys.EnergyIntComp[0] += P[i].Mass * egyspec;
#ifdef EXTRA_STATS
#if defined BETA_COOLING && defined OUTPUTRADIATEDENERGY
     sys.EnergyRadComp[0] += P[i].Mass * SphP[i].RadiatedEnergy;
#endif
#endif
	}



      for(j = 0; j < 3; j++)
	{
	  sys.MomentumComp[P[i].Type][j] += P[i].Mass * vel[j];
	  sys.CenterOfMassComp[P[i].Type][j] += P[i].Mass * P[i].Pos[j];
	}

      sys.AngMomentumComp[P[i].Type][0] += P[i].Mass * (P[i].Pos[1] * vel[2] - P[i].Pos[2] * vel[1]);
      sys.AngMomentumComp[P[i].Type][1] += P[i].Mass * (P[i].Pos[2] * vel[0] - P[i].Pos[0] * vel[2]);
      sys.AngMomentumComp[P[i].Type][2] += P[i].Mass * (P[i].Pos[0] * vel[1] - P[i].Pos[1] * vel[0]);
    }


  /* some the stuff over all processors */
  MPI_Reduce(&sys.MassComp[0], &SysState.MassComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyPotComp[0], &SysState.EnergyPotComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyIntComp[0], &SysState.EnergyIntComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sys.EnergyKinComp[0], &SysState.EnergyKinComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef EXTRA_STATS
  MPI_Reduce(&sys.EnergyRadComp[0], &SysState.EnergyRadComp[0], 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  MPI_Reduce(&sys.MomentumComp[0][0], &SysState.MomentumComp[0][0], 6 * 4, MPI_DOUBLE, MPI_SUM, 0,
	     MPI_COMM_WORLD);
  MPI_Reduce(&sys.AngMomentumComp[0][0], &SysState.AngMomentumComp[0][0], 6 * 4, MPI_DOUBLE, MPI_SUM, 0,
	     MPI_COMM_WORLD);
  MPI_Reduce(&sys.CenterOfMassComp[0][0], &SysState.CenterOfMassComp[0][0], 6 * 4, MPI_DOUBLE, MPI_SUM, 0,
	     MPI_COMM_WORLD);


  if(ThisTask == 0)
    {
      for(i = 0; i < 6; i++)
      {
	SysState.EnergyTotComp[i] = SysState.EnergyKinComp[i] +
	  SysState.EnergyPotComp[i] + SysState.EnergyIntComp[i];
#ifdef EXTRA_STATS
   SysState.EnergyTotComp[i] += SysState.EnergyRadComp[i];
#endif
      }

      SysState.Mass = SysState.EnergyKin = SysState.EnergyPot = SysState.EnergyInt = SysState.EnergyTot = 0;
#ifdef EXTRA_STATS
      SysState.EnergyRad = 0;
#endif


      for(j = 0; j < 3; j++)
	SysState.Momentum[j] = SysState.AngMomentum[j] = SysState.CenterOfMass[j] = 0;

      for(i = 0; i < 6; i++)
	{
	  SysState.Mass += SysState.MassComp[i];
	  SysState.EnergyKin += SysState.EnergyKinComp[i];
	  SysState.EnergyPot += SysState.EnergyPotComp[i];
	  SysState.EnergyInt += SysState.EnergyIntComp[i];
	  SysState.EnergyTot += SysState.EnergyTotComp[i];
#ifdef EXTRA_STATS
     SysState.EnergyRad += SysState.EnergyRadComp[i];
#endif

	  for(j = 0; j < 3; j++)
	    {
	      SysState.Momentum[j] += SysState.MomentumComp[i][j];
	      SysState.AngMomentum[j] += SysState.AngMomentumComp[i][j];
	      SysState.CenterOfMass[j] += SysState.CenterOfMassComp[i][j];
	    }
	}

      for(i = 0; i < 6; i++)
	for(j = 0; j < 3; j++)
	  if(SysState.MassComp[i] > 0)
	    SysState.CenterOfMassComp[i][j] /= SysState.MassComp[i];

      for(j = 0; j < 3; j++)
	if(SysState.Mass > 0)
	  SysState.CenterOfMass[j] /= SysState.Mass;

      for(i = 0; i < 6; i++)
	{
	  SysState.CenterOfMassComp[i][3] = SysState.MomentumComp[i][3] = SysState.AngMomentumComp[i][3] = 0;
	  for(j = 0; j < 3; j++)
	    {
	      SysState.CenterOfMassComp[i][3] +=
		SysState.CenterOfMassComp[i][j] * SysState.CenterOfMassComp[i][j];
	      SysState.MomentumComp[i][3] += SysState.MomentumComp[i][j] * SysState.MomentumComp[i][j];
	      SysState.AngMomentumComp[i][3] +=
		SysState.AngMomentumComp[i][j] * SysState.AngMomentumComp[i][j];
	    }
	  SysState.CenterOfMassComp[i][3] = sqrt(SysState.CenterOfMassComp[i][3]);
	  SysState.MomentumComp[i][3] = sqrt(SysState.MomentumComp[i][3]);
	  SysState.AngMomentumComp[i][3] = sqrt(SysState.AngMomentumComp[i][3]);
	}

      SysState.CenterOfMass[3] = SysState.Momentum[3] = SysState.AngMomentum[3] = 0;

      for(j = 0; j < 3; j++)
	{
	  SysState.CenterOfMass[3] += SysState.CenterOfMass[j] * SysState.CenterOfMass[j];
	  SysState.Momentum[3] += SysState.Momentum[j] * SysState.Momentum[j];
	  SysState.AngMomentum[3] += SysState.AngMomentum[j] * SysState.AngMomentum[j];
	}

      SysState.CenterOfMass[3] = sqrt(SysState.CenterOfMass[3]);
      SysState.Momentum[3] = sqrt(SysState.Momentum[3]);
      SysState.AngMomentum[3] = sqrt(SysState.AngMomentum[3]);
    }

  /* give everyone the result, maybe the want to do something with it */
  MPI_Bcast(&SysState, sizeof(struct state_of_system), MPI_BYTE, 0, MPI_COMM_WORLD);
}

#ifdef DEAD_GAS
//Disables the particles that aren't ready to have forces calculated on them yet
void kill_gas(void)
{
  int i,j;
  double r2,d2;

  d2=All.Doom_radius*All.Doom_radius;
  for(i=0;i<N_gas;i++)
  {
    r2=(P[i].Pos[0]-SysState.CenterOfMass[0])*(P[i].Pos[0]-SysState.CenterOfMass[0])+(P[i].Pos[1]-SysState.CenterOfMass[1])*(P[i].Pos[1]-SysState.CenterOfMass[1])+(P[i].Pos[2]-SysState.CenterOfMass[2])*(P[i].Pos[2]-SysState.CenterOfMass[2]);
    if(r2>d2 && P[i].Ti_endstep>=0)
    {
      P[i].Ti_endstep = -P[i].Ti_endstep -1;
    }
    if(r2>d2)
    {
      for(j=0;j<3;j++)
      {
        P[i].GravAccel[j]=0;
        SphP[i].HydroAccel[j]=0;
      }
    }
    else
    {
      //printf("Particle %d at %g,%g,%g is alive with drift %g!\n",i,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],sqrt(P[i].Vel[0]*P[i].Vel[0]+P[i].Vel[1]*P[i].Vel[1]+P[i].Vel[2]*P[i].Vel[2]));
    }
  }
}

//Sometimes there's no ressurection built into the routines
void ressurect(void)
{
  int i;
  for(i=0;i<N_gas;i++)
    if(P[i].Ti_endstep<0)
      P[i].Ti_endstep = -P[i].Ti_endstep -1;
}
#endif
