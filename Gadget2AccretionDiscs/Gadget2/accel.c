#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*! \file accel.c
 *  \brief driver routine to carry out force computation
 */


/*! This routine computes the accelerations for all active particles.
 *  First, the long-range PM force is computed if the TreePM algorithm is
 *  used and a "big" PM step is done.  Next, the gravitational tree forces
 *  are computed. This also constructs the tree, if needed.
 *
 *  If gas particles are present, the density-loop for active SPH particles
 *  is carried out. This includes an iteration on the correct number of
 *  neighbours.  Finally, the hydrodynamical forces are added.
 */
void compute_accelerations(int mode)
{
  double tstart, tend;
#ifdef SURFACE
  int i;
  double dt_entr;
#endif

  if(ThisTask == 0)
    {
      printf("Start force computation...\n");
      fflush(stdout);
    }

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      tstart = second();
      long_range_force();
      tend = second();
      All.CPU_PM += timediff(tstart, tend);
    }
#endif

  tstart = second();		/* measure the time for the full force computation */
#if defined TWODIMS || defined SURFACE
  //Make sure we have a currentish estimate of the centre of mass and the energy of the particles
  //Needs to happen before gravity if we're going to smooth using H which depends on the entropy
  if(All.Ti_Current ==0)
  {
    density();
    compute_global_quantities_of_system();
  }
#endif

  gravity_tree();		/* computes gravity accel. */

  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
    gravity_tree();		/* For the first timestep, we redo it
				 * to allow usage of relative opening
				 * criterion for consistent accuracy.
				 */
  tend = second();
  All.CPU_Gravity += timediff(tstart, tend);

#ifdef FORCETEST
  gravity_forcetest();
#endif
#ifdef DIRTY_SHAMEFUL_SECRET
  if(NumPart-N_gas){
    P[N_gas].GravAccel[0]=0;
    P[N_gas].GravAccel[1]=0;
    P[N_gas].GravAccel[2]=0;
  }
#endif

  if(All.TotN_gas > 0)
    {
      if(ThisTask == 0)
	{
	  printf("Start density computation...\n");
	  fflush(stdout);
	}

      tstart = second();
      density();		/* computes density, and pressure */
      tend = second();
#ifdef SURFACE
      double tot,bigtot;
      //Want to do it after density, so we have accurate density values for calculating energy
      sur_density();
      //Pre update value
      tot=0;
      for(i=0;i<N_gas;i++)
      {
        dt_entr = (All.Ti_Current - .5*(P[i].Ti_endstep+P[i].Ti_begstep)) * All.Timebase_interval;
        tot += P[i].Mass * (SphP[i].Entropy+dt_entr*SphP[i].DtEntropy)*pow(SphP[i].Density,GAMMA_MINUS1) / GAMMA_MINUS1;
      }
      MPI_Allreduce(&tot,&bigtot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      if(ThisTask==0)
        printf("Pre update total energy = %g.\n",bigtot);

      for(i=0;i<N_gas;i++)
      {
        if(P[i].Ti_endstep == All.Ti_Current)
        {
          //Set the pressure to the new value
          SphP[i].Pressure = SphP[i].SurEntropy * pow(SphP[i].Density,GAMMA);
          //Set the entropy such that K+dK*dt = K required for current internal energy
          dt_entr = (All.Ti_Current - .5*(P[i].Ti_endstep+P[i].Ti_begstep)) * All.Timebase_interval;
          SphP[i].Entropy = (SphP[i].SurEntropy * GAMMA_MINUS1 / pow(SphP[i].Density,GAMMA_MINUS1)) - SphP[i].DtEntropy*dt_entr;
        }
      }

      //Post update value
      tot=0;
      for(i=0;i<N_gas;i++)
      {
        dt_entr = (All.Ti_Current - .5*(P[i].Ti_endstep+P[i].Ti_begstep)) * All.Timebase_interval;
        tot += P[i].Mass * (SphP[i].Entropy+dt_entr*SphP[i].DtEntropy)*pow(SphP[i].Density,GAMMA_MINUS1) / GAMMA_MINUS1;
      }
      MPI_Allreduce(&tot,&bigtot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      if(ThisTask==0)
        printf("Post update total energy = %g.\n",bigtot);
#endif

      All.CPU_Hydro += timediff(tstart, tend);

      tstart = second();
      force_update_hmax();      /* tell the tree nodes the new SPH smoothing length such that they are guaranteed to hold the correct max(Hsml) */
      tend = second();
      All.CPU_Predict += timediff(tstart, tend);


      if(ThisTask == 0)
	{
	  printf("Start hydro-force computation...\n");
	  fflush(stdout);
	}

      tstart = second();
      hydro_force();		/* adds hydrodynamical accelerations and computes viscous entropy injection  */
      tend = second();
      All.CPU_Hydro += timediff(tstart, tend);

    }

  if(ThisTask == 0)
    {
      printf("force computation done.\n");
      fflush(stdout);
    }

}
