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

#ifdef DEAD_GAS
  if(All.Ti_Current ==0)
    compute_global_quantities_of_system();
  kill_gas();
#endif
  gravity_tree();		/* computes gravity accel. */

  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
  {
#ifdef DEAD_GAS
    kill_gas();
#endif
    gravity_tree();		/* For the first timestep, we redo it
				 * to allow usage of relative opening
				 * criterion for consistent accuracy.
				 */
  }
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
#ifdef SURFACE
      //Calculate the surface quantities and smooth out the energy
      sur_density();
      //What is the total energy, pre spread
      for(i=0;i<N_gas;i++)
        SphP[i].Entropy=SphP[i].SurEntropy;
#endif
#ifdef DEAD_GAS
      kill_gas();
#endif
      density();		/* computes density, and pressure */
      tend = second();
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
#ifdef DEAD_GAS
      kill_gas();
#endif
      hydro_force();		/* adds hydrodynamical accelerations and computes viscous entropy injection  */
#ifdef DEAD_GAS
      ressurect();
#endif
      tend = second();
      All.CPU_Hydro += timediff(tstart, tend);
    }

  if(ThisTask == 0)
    {
      printf("force computation done.\n");
      fflush(stdout);
    }

}
