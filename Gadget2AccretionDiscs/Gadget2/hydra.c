#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include "allvars.h"
#include "proto.h"

/*! \file hydra.c
 *  \brief Computation of SPH forces and rate of entropy generation
 *
 *  This file contains the "second SPH loop", where the SPH forces are
 *  computed, and where the rate of change of entropy due to the shock heating
 *  (via artificial viscosity) is computed.
 */


static double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;

#ifdef PERIODIC
static double boxSize, boxHalf;

#ifdef LONG_X
static double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
static double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
static double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif

#ifdef BETA_COOLING
static double starData[7];
#endif



/*! This function is the driver routine for the calculation of hydrodynamical
 *  force and rate of change of entropy due to shock heating for all active
 *  particles .
 */
void hydro_force(void)
{
  long long ntot, ntotleft;
  int i, j, k, n, ngrp, maxfill, source, ndone;
  int *nbuffer, *noffset, *nsend_local, *nsend, *numlist, *ndonelist;
  int level, sendTask, recvTask, nexport, place;
  double soundspeed_i;
  double tstart, tend, sumt, sumcomm;
  double timecomp = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  MPI_Status status;
#ifdef BETA_COOLING
  double tdyn,E,R,v2;
  int numsinks,root,globalroot;
#endif

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) + All.OmegaLambda;

      hubble_a = All.Hubble * sqrt(hubble_a);
      hubble_a2 = All.Time * All.Time * hubble_a;

      fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;

      fac_egy = pow(All.Time, 3 * (GAMMA - 1));

      fac_vsic_fix = hubble_a * pow(All.Time, 3 * GAMMA_MINUS1);

      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    hubble_a = hubble_a2 = atime = fac_mu = fac_vsic_fix = a3inv = fac_egy = 1.0;


  /* `NumSphUpdate' gives the number of particles on this processor that want a force update */
  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      if(P[n].Ti_endstep == All.Ti_Current)
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);

#ifdef BETA_COOLING
  /* Get the position and mass of the central object and send it to everyone */
  numsinks=NumPart - N_gas;
  starData[0]=starData[1]=starData[2]=starData[3]= -1.0;
  root=-1;
  for(i=0; i<numsinks;i++)
  {
    if(P[i+N_gas].ID==All.StarID)
    {
      starData[0] = P[i+N_gas].Pos[0];
      starData[1] = P[i+N_gas].Pos[1];
      starData[2] = P[i+N_gas].Pos[2];
      starData[3] = P[i+N_gas].Mass;
      //This needs to be changed to the PREDICTED velocity
      starData[4] = P[i+N_gas].Vel[0];
      starData[5] = P[i+N_gas].Vel[1];
      starData[6] = P[i+N_gas].Vel[2];
      root = ThisTask;
    }
  }
  /* Get the node that has the data */
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&root,&globalroot,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  /* Broadcast it. */
  MPI_Bcast(&starData,7,MPI_DOUBLE,globalroot,MPI_COMM_WORLD);
  //printf("The star ID is %d. The position is (%g,%g,%g) and the mass is %g.\n",All.StarID,starData[0],starData[1],starData[2],starData[3]);
#endif

  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);


  i = 0;			/* first particle for this task */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeHydro - NTask; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    ndone++;

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;

	    hydro_evaluate(i, 0);

	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    for(k = 0; k < 3; k++)
		      {
			HydroDataIn[nexport].Pos[k] = P[i].Pos[k];
			HydroDataIn[nexport].Vel[k] = SphP[i].VelPred[k];
		      }
		    HydroDataIn[nexport].Hsml = SphP[i].Hsml;
		    HydroDataIn[nexport].Mass = P[i].Mass;
		    HydroDataIn[nexport].DhsmlDensityFactor = SphP[i].DhsmlDensityFactor;
		    HydroDataIn[nexport].Density = SphP[i].Density;
		    HydroDataIn[nexport].Pressure = SphP[i].Pressure;
		    HydroDataIn[nexport].Timestep = P[i].Ti_endstep - P[i].Ti_begstep;

		    /* calculation of F1 */
		    soundspeed_i = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
#ifdef NOBALSARA
          HydroDataIn[nexport].F1 = 1.0;
#else
		    HydroDataIn[nexport].F1 = fabs(SphP[i].DivVel) /
		      (fabs(SphP[i].DivVel) + SphP[i].CurlVel +
		       0.0001 * soundspeed_i / SphP[i].Hsml / fac_mu);
#endif

		    HydroDataIn[nexport].Index = i;
		    HydroDataIn[nexport].Task = j;
        
#if defined MMAV || defined CDAV
        HydroDataIn[nexport].Alpha = SphP[i].Alpha;
#endif
#ifdef SINK_PARTICLES
        HydroDataIn[nexport].AccretionTarget = SphP[i].AccretionTarget;
#endif	
#ifdef PRICE_GRAV_SOFT
        HydroDataIn[nexport].Zeta = SphP[i].Zeta;
#endif
        
		    nexport++;
		    nsend_local[j]++;
		  }
	      }
	  }
      tend = second();
      timecomp += timediff(tstart, tend);

      qsort(HydroDataIn, nexport, sizeof(struct hydrodata_in), hydro_compare_key);

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	     noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

      tend = second();
      timeimbalance += timediff(tstart, tend);



      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	   {
	     tstart = second();
	     for(j = 0; j < NTask; j++)
	       nbuffer[j] = 0;
	     for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	     {
	       maxfill = 0;
	       for(j = 0; j < NTask; j++)
		    {
		      if((j ^ ngrp) < NTask)
		        if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		          maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		    }
	       if(maxfill >= All.BunchSizeHydro)
		      break;

	       sendTask = ThisTask;
	       recvTask = ThisTask ^ ngrp;

	       if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		      {
		        /* get the particles */
		        MPI_Sendrecv(&HydroDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A,
				   &HydroDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, &status);
		      }
		    }

	       for(j = 0; j < NTask; j++)
		      if((j ^ ngrp) < NTask)
		        nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	     }
	     tend = second();
	     timecommsumm += timediff(tstart, tend);

	     /* now do the imported particles */
	     tstart = second();
	     for(j = 0; j < nbuffer[ThisTask]; j++)
          hydro_evaluate(j, 1);
	     tend = second();
	     timecomp += timediff(tstart, tend);

	     /* do a block to measure imbalance */
	     tstart = second();
	     MPI_Barrier(MPI_COMM_WORLD);
	     tend = second();
	     timeimbalance += timediff(tstart, tend);

	     /* get the result */
	     tstart = second();
	     for(j = 0; j < NTask; j++)
	       nbuffer[j] = 0;
	     for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	     {
	       maxfill = 0;
	       for(j = 0; j < NTask; j++)
		    {
		      if((j ^ ngrp) < NTask)
		        if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		          maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		    }
	       if(maxfill >= All.BunchSizeHydro)
		      break;

	       sendTask = ThisTask;
	       recvTask = ThisTask ^ ngrp;

	       if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		      {
		        /* send the results */
		        MPI_Sendrecv(&HydroDataResult[nbuffer[ThisTask]],
		    	   nsend[recvTask * NTask + ThisTask] * sizeof(struct hydrodata_out),
		    	   MPI_BYTE, recvTask, TAG_HYDRO_B,
		    	   &HydroDataPartialResult[noffset[recvTask]],
		    	   nsend_local[recvTask] * sizeof(struct hydrodata_out),
		    	   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, &status);

		        /* add the result to the particles */
		        for(j = 0; j < nsend_local[recvTask]; j++)
		        {
		          source = j + noffset[recvTask];
		          place = HydroDataIn[source].Index;

		          for(k = 0; k < 3; k++)
		            SphP[place].HydroAccel[k] += HydroDataPartialResult[source].Acc[k];
  
     	    	   SphP[place].DtEntropy += HydroDataPartialResult[source].DtEntropy;
    
		          if(SphP[place].MaxSignalVel < HydroDataPartialResult[source].MaxSignalVel)
		            SphP[place].MaxSignalVel = HydroDataPartialResult[source].MaxSignalVel;
		        }
		      }
		    }

	       for(j = 0; j < NTask; j++)
		      if((j ^ ngrp) < NTask)
		        nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	     }
	     tend = second();
	     timecommsumm += timediff(tstart, tend);

	     level = ngrp - 1;
	  }

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);



  /* do final operations on results */
  tstart = second();

  for(i = 0; i < N_gas; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      {
	SphP[i].DtEntropy *= GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].Density, GAMMA_MINUS1));
#ifdef BETA_COOLING
   //The conversion factor from dK/dt to du/dt cancel with the conversion factor from u to K in du/dt = -u / beta/Omega.  This is why we do the above line before the cooling...
   for(j=0,R=0,v2=0;j<3;j++)
   {
     R+=(P[i].Pos[j]-starData[j])*(P[i].Pos[j]-starData[j]);
     v2+=SphP[i].VelPred[j]*SphP[i].VelPred[j];
   }
   R=sqrt(R);
//E = (v2/2.0) - ((All.G*starData[3]) / (R * P[i].Mass));
   //tdyn=sqrt(-8.0*E*E*E)/(All.G*starData[3]*All.G*starData[3]);
   tdyn=sqrt(R*R*R/All.G/starData[3]);
   SphP[i].DtEntropy -= SphP[i].Entropy / All.CoolingRate / tdyn;
#ifdef OUTPUTRADIATEDENERGY
   //We want to track the ENERGY radiated, not the ENTROPY
   SphP[i].DtRadiatedEnergy = (SphP[i].Entropy*pow(SphP[i].Density,GAMMA_MINUS1)) / (All.CoolingRate * tdyn * (GAMMA_MINUS1));
#endif
#endif
   //These values will be needed by the density loop the next time around...
   //Not done in the density loop since the density loop could in principal be run multiple times per particle and these values should only be set once to the final results...
#ifdef CDAV
   SphP[i].DivVelOld = SphP[i].DivVel;
   P[i].GravAccelOld[0] = P[i].GravAccel[0];
   P[i].GravAccelOld[1] = P[i].GravAccel[1];
   P[i].GravAccelOld[2] = P[i].GravAccel[2];
#endif
#ifdef SPH_BND_PARTICLES
	if(P[i].ID == 0)
	  {
	    SphP[i].DtEntropy = 0;
	    for(k = 0; k < 3; k++)
	      SphP[i].HydroAccel[k] = 0;
	  }
#endif
      }

  tend = second();
  timecomp += timediff(tstart, tend);

  /* collect some timing information */

  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += sumimbalance / NTask;
    }
}


/*! This function is the 'core' of the SPH force computation. A target
 *  particle is specified which may either be local, or reside in the
 *  communication buffer.
 */
void hydro_evaluate(int target, int mode)
{
  int j, k, n, timestep, startnode, numngb;
  FLOAT *pos, *vel;
  FLOAT mass, h_i, dhsmlDensityFactor, rho, pressure, f1, f2;
#if defined MMAV || defined CDAV
  FLOAT alpha_visc,alpha_visc_j;
  double tmp;
#endif
  double acc[3], dtEntropy, maxSignalVel;
  double dx, dy, dz, dvx, dvy, dvz;
  double h_i2, hinv, hinv4;
  double p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;
  double hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;
  double h_j, dwk_j, r, r2, u, hfc_visc;

#ifndef NOVISCOSITYLIMITER
  double dt;
#endif
#ifdef PRICE_GRAV_SOFT
  double zeta=0;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = SphP[target].Hsml;
      mass = P[target].Mass;
      dhsmlDensityFactor = SphP[target].DhsmlDensityFactor;
      rho = SphP[target].Density;
      pressure = SphP[target].Pressure;
      timestep = P[target].Ti_endstep - P[target].Ti_begstep;
      soundspeed_i = sqrt(GAMMA * pressure / rho);
#if defined MMAV || defined CDAV
      alpha_visc = SphP[target].Alpha;
#endif
#ifdef PRICE_GRAV_SOFT
      zeta = SphP[target].Zeta;
#endif
#ifdef NOBALSARA
      f1 = 1.0;
#else
      f1 = fabs(SphP[target].DivVel) /
	(fabs(SphP[target].DivVel) + SphP[target].CurlVel +
	 0.0001 * soundspeed_i / SphP[target].Hsml / fac_mu);
#endif
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;
      mass = HydroDataGet[target].Mass;
      dhsmlDensityFactor = HydroDataGet[target].DhsmlDensityFactor;
      rho = HydroDataGet[target].Density;
      pressure = HydroDataGet[target].Pressure;
      timestep = HydroDataGet[target].Timestep;
      soundspeed_i = sqrt(GAMMA * pressure / rho);
      f1 = HydroDataGet[target].F1;
#if defined MMAV || defined CDAV
      alpha_visc = HydroDataGet[target].Alpha;
#endif
#ifdef PRICE_GRAV_SOFT
      zeta = HydroDataGet[target].Zeta;
#endif
    }


  /* initialize variables before SPH loop is started */
  acc[0] = acc[1] = acc[2] = dtEntropy = 0;
  maxSignalVel = 0;
//  //I think this should really always be here, even for the default GADGET viscosity, but it's unlikely to ever really make much difference.  Because it's not here in the default gadget version, we'll leave it out altogether for consistency...
//#if defined CDAV || defined MMAV
//  //Because the r=0 case will not be evaluated below...
//  maxSignalVel = soundspeed_i;
//#endif

  p_over_rho2_i = pressure / (rho * rho) * dhsmlDensityFactor;
#ifdef PRICE_GRAV_SOFT
  p_over_rho2_i = dhsmlDensityFactor * ((pressure/(rho*rho))+(zeta*All.G*.5));
#endif
  h_i2 = h_i * h_i;
  

  /* Now start the actual SPH computation for this particle */
  startnode = All.MaxPart;
  do
    {
      numngb = ngb_treefind_pairs(&pos[0], h_i, &startnode);

      for(n = 0; n < numngb; n++)
	{
	  j = Ngblist[n];

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#if defined MMAV || defined CDAV
     alpha_visc_j = SphP[j].Alpha;
#endif

#ifdef PERIODIC			/*  find the closest image in the given box size  */
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;
	  h_j = SphP[j].Hsml;
	  if(r2 < h_i2 || r2 < h_j * h_j)
	    {
	      r = sqrt(r2);
	      if(r > 0)
		{
		  p_over_rho2_j = SphP[j].Pressure / (SphP[j].Density * SphP[j].Density);
		  soundspeed_j = sqrt(GAMMA * p_over_rho2_j * SphP[j].Density);
#ifdef PRICE_GRAV_SOFT
        p_over_rho2_j += 0.5 * All.G * SphP[j].Zeta;
#endif
		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];
		  vdotr = dx * dvx + dy * dvy + dz * dvz;

		  if(All.ComovingIntegrationOn)
		    vdotr2 = vdotr + hubble_a2 * r2;
		  else
		    vdotr2 = vdotr;

		  if(r2 < h_i2)
		    {
		      hinv = 1.0 / h_i;
#ifndef  TWODIMS
		      hinv4 = hinv * hinv * hinv * hinv;
#else
		      //hinv4 = hinv * hinv * hinv / boxSize_Z;
		      hinv4 = hinv * hinv * hinv;
#endif
		      u = r * hinv;
		      if(u < 0.5)
            {
			     dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
            }
		      else
            {
			     dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
            }
		    }
		  else
		    {
		      dwk_i = 0;
		    }

		  if(r2 < h_j * h_j)
		    {
		      hinv = 1.0 / h_j;
#ifndef  TWODIMS
		      hinv4 = hinv * hinv * hinv * hinv;
#else
		      //hinv4 = hinv * hinv * hinv / boxSize_Z;
		      hinv4 = hinv * hinv * hinv;
#endif
		      u = r * hinv;
		      if(u < 0.5)
			dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		      else
			dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		    }
		  else
		    {
		      dwk_j = 0;
		    }

		  if(soundspeed_i + soundspeed_j > maxSignalVel)
		    maxSignalVel = soundspeed_i + soundspeed_j;

        //This is the standard AV block, unchanged from original GADGET...

//#ifndef CDAV
//#ifndef MMAV
		  if(vdotr2 < 0)	/* ... artificial viscosity */
		    {
		      mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */

//#if defined CDAV || defined MMAV
//            vsig = soundspeed_i + soundspeed_j - 2.0 * All.ArtViscPropConst*mu_ij;
//#else
                      vsig = soundspeed_i + soundspeed_j - 3 * mu_ij;
//#endif
//            vsig = soundspeed_i + soundspeed_j - 2.0 * All.ArtViscPropConst*mu_ij;


		      if(vsig > maxSignalVel)
			maxSignalVel = vsig;

		      rho_ij = 0.5 * (rho + SphP[j].Density);
#ifdef NOBALSARA
            f2 = 1.0;
#else
		      f2 =
			fabs(SphP[j].DivVel) / (fabs(SphP[j].DivVel) + SphP[j].CurlVel +
						0.0001 * soundspeed_j / fac_mu / SphP[j].Hsml);
#endif
#if defined CDAV || defined MMAV
            //printf("[%d] Alpha i/j = %g/%g.\n",ThisTask,alpha_visc,alpha_visc_j);
            visc = 0.25 * (alpha_visc + alpha_visc_j) * (soundspeed_i + soundspeed_j-2.0*All.ArtViscPropConst*mu_ij) * (-mu_ij) / rho_ij;
#else
		      visc = 0.25 * All.ArtBulkViscConst * (soundspeed_i + soundspeed_j-3*mu_ij) * (-mu_ij) / rho_ij * (f1 + f2);
#endif

		      /* .... end artificial viscosity evaluation */
#ifndef NOVISCOSITYLIMITER
		      /* make sure that viscous acceleration is not too large */
		      dt = imax(timestep, (P[j].Ti_endstep - P[j].Ti_begstep)) * All.Timebase_interval;
		      if(dt > 0 && (dwk_i + dwk_j) < 0)
			{
			  visc = dmin(visc, 0.5 * fac_vsic_fix * vdotr2 /
				      (0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
			}
#endif
		    }
		  else
		    visc = 0;
#ifdef ART_COND
        visc -= All.ArtCondConst * sqrt(fabs(pressure - SphP[j].Pressure)/rho_ij) * GAMMA
          * ((pressure/rho) - (SphP[j].Pressure/SphP[j].Density));
#endif

         // if(ThisTask==1)
         // {
         //   printf("Default Visc is %g\n",visc);
         // }

		  p_over_rho2_j *= SphP[j].DhsmlDensityFactor;

		  hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;

		  hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;

		  acc[0] -= hfc * dx;
		  acc[1] -= hfc * dy;
		  acc[2] -= hfc * dz;
		  dtEntropy += 0.5 * hfc_visc * vdotr2;
//#endif
//#endif
//The block for CD Method
//This block should only be used if you want to use the explicit CD form for the AV and for the signal velocity/dissipation
//#if defined CDAV || defined CDAV_DRIFTUPDATE
//        //Set the Signal velocity, only if within h_i
//        if(dwk_i!=0)
//        {
//          if(maxSignalVel < .5*(soundspeed_i + soundspeed_j) - dmin(0,vdotr2)/r)
//          {
//            maxSignalVel = .5*(soundspeed_i + soundspeed_j) -dmin(0,vdotr2)/r;
//          }
//        }
//        if(vdotr2 < 0)
//        {
//            //We are explicitly putting in the form of the C&D paper...
//            mu_ij = (4*vdotr2*h_j*h_j*h_i*h_i)/(r*r*(h_i*h_i+h_j*h_j)*(h_i+h_j));
//            //mu_ij is now Pi_ij
//            mu_ij = -mu_ij*(0.5*(soundspeed_i+soundspeed_j)-All.ArtViscPropConst*mu_ij);
//            //Calculate the acceleration...
//            hfc_visc = .5*mu_ij*P[j].Mass*((alpha_visc*dhsmlDensityFactor*dwk_i/rho)+(alpha_visc_j*SphP[j].DhsmlDensityFactor*dwk_j/SphP[j].Density))/r;
//            //Calculate the dispersion...
//            hfc = .5*vdotr2*P[j].Mass * mu_ij*alpha_visc*dhsmlDensityFactor*dwk_i/(rho*r);
//#ifndef NOVISCOSITYLIMITER
//		      dt = imax(timestep, (P[j].Ti_endstep - P[j].Ti_begstep)) * All.Timebase_interval;
//		      if(dt > 0 && (dwk_i + dwk_j) < 0)
//            {
//              tmp = .25 * fac_vsic_fix * vdotr2 / (r*r*dt);
//              if(hfc_visc<tmp)
//              {
//                printf("Viscosity limiter invoked for particle %d!\n",j);
//                hfc_visc = tmp;
//                hfc = vdotr2 * hfc_visc*.5;
//              }
//            }
//#endif
//        }
//        else
//        {
//          hfc_visc=hfc=0;
//        }
//        dtEntropy+=hfc;
//#endif
//Block for the MM method
//#if defined MMAV || defined CDAV
//
//        if(vdotr2 <0)
//        {
//          mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */
//          //ArtViscPropConst is 3/2 in original implementation...
//          //vsig = soundspeed_i + soundspeed_j - All.ArtViscPropConst*2.0 * mu_ij;
//          vsig = soundspeed_i + soundspeed_j - 3 * mu_ij;
//
//	       if(vsig > maxSignalVel)
//	         maxSignalVel = vsig;
//
//		    rho_ij = 0.5 * (rho + SphP[j].Density);
//
//
//          //if(ThisTask==0)
//          //{
//          //  printf("Average Alpha is %g\n",.5*(alpha_visc+alpha_visc_j));
//          //}
//
//          //We don't just use the signal velocity here as we want to allow beta to be something other than 3/2 alpha.
//          visc = 0.25 * (alpha_visc + alpha_visc_j) * (soundspeed_i+soundspeed_j-2.0*All.ArtViscPropConst*mu_ij) * (-mu_ij) / rho_ij;
//#ifndef NOVISCOSITYLIMITER
//          /* make sure that viscous acceleration is not too large */
//		    dt = imax(timestep, (P[j].Ti_endstep - P[j].Ti_begstep)) * All.Timebase_interval;
//		    if(dt > 0 && (dwk_i + dwk_j) < 0)
//          {
//			  visc = dmin(visc, 0.5 * fac_vsic_fix * vdotr2 /
//				      (0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
//
//            //Limiter is designed to set the magnitude of the viscous acceleration to be no more than a quarter of the line of sight velocity connecting i and j, divided by dt
//            //tmp=0.5*fac_vsic_fix * vdotr2 / (0.5 * (mass+P[j].Mass)*(dwk_i+dwk_j)*r*dt);
//            //if(visc>tmp)
//            //{
//            //  visc = tmp;
//            //}
//          }
//#endif
//        }
//        else
//        {
//          visc=0;
//        }
//        //if(ThisTask==0)
//        //{
//        //  printf("MM Visc alpha_i,alpha_j =  %g  %g  %g\n",visc,alpha_visc,alpha_visc_j);
//        //}
//
//		  p_over_rho2_j *= SphP[j].DhsmlDensityFactor;
//
//		  hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;
//
//		  hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;
//
//		  acc[0] -= hfc * dx;
//		  acc[1] -= hfc * dy;
//		  acc[2] -= hfc * dz;
//		  dtEntropy += 0.5 * hfc_visc * vdotr2;
//#endif
		}
	    }
	}
    }
  while(startnode >= 0);

  /* Now collect the result at the right place */
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	     SphP[target].HydroAccel[k] = acc[k];
      SphP[target].DtEntropy = dtEntropy;
      SphP[target].MaxSignalVel = maxSignalVel;
    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].Acc[k] = acc[k];
      HydroDataResult[target].DtEntropy = dtEntropy;
      HydroDataResult[target].MaxSignalVel = maxSignalVel;
    }
}




/*! This is a comparison kernel for a sort routine, which is used to group
 *  particles that are going to be exported to the same CPU.
 */
int hydro_compare_key(const void *a, const void *b)
{
  if(((struct hydrodata_in *) a)->Task < (((struct hydrodata_in *) b)->Task))
    return -1;
  if(((struct hydrodata_in *) a)->Task > (((struct hydrodata_in *) b)->Task))
    return +1;
  return 0;
}
