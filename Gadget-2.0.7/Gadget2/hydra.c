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



/*! This function is the driver routine for the calculation of hydrodynamical
 *  force and rate of change of entropy due to shock heating for all active
 *  particles .
 *  If individual particle viscosity is turned on, then this function performs 
 *  two loops, the first to calculate the updated viscosities, the final one to 
 *  calculate the viscous acceleration.
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
#ifdef INDIVIDUALAV
  double fac,Tinv[9],V[9],divv,tmp,trSSt,xi_i,dt,divdotv;
  double alpha_loc,tau_i,A_i,R_i;
#endif
  MPI_Status status;

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

  //Reasonably sure you could malloc this to be just NTask * sizeof(int)
  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


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
		    HydroDataIn[nexport].F1 = fabs(SphP[i].DivVel) /
		      (fabs(SphP[i].DivVel) + SphP[i].CurlVel +
		       0.0001 * soundspeed_i / SphP[i].Hsml / fac_mu);

		    HydroDataIn[nexport].Index = i;
		    HydroDataIn[nexport].Task = j;
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
           {
			    SphP[place].HydroAccel[k] += HydroDataPartialResult[source].Acc[k];
#ifdef ARTVISCTEST
             SphP[place].OldArtViscAccel[k] += HydroDataPartialResult[source].OldArtViscAccel[k];
#endif
           }
           for(k=0;k<9;k++)
           {
             SphP[place].MatrixD[k] += HydroDataPartialResult[source].MatrixD[k];
             SphP[place].MatrixT[k] += HydroDataPartialResult[source].MatrixT[k];
           }

           SphP[place].R_i += HydroDataPartialResult[source].R_i;
			  SphP[place].DtEntropy += HydroDataPartialResult[source].DtEntropy;

			  if(SphP[place].MaxSignalVel < HydroDataPartialResult[source].MaxSignalVel)
			    SphP[place].MaxSignalVel = HydroDataPartialResult[source].MaxSignalVel;
           if(SphP[place].SignalVel < HydroDataPartialResult[source].SignalVel)
             SphP[place].SignalVel = HydroDataPartialResult[source].SignalVel;
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

#ifdef INDIVIDUALAV
   //Calculate the per particle artificial viscosity and the change in entropy.  Needs to be done here as calculating the artificial viscosity parameter requires a loop over neighbours AFTER calculating density.
   //Invert the matrix T
   fac = SphP[i].MatrixT[0]*(SphP[i].MatrixT[4]*SphP[i].MatrixT[8]-SphP[i].MatrixT[5]*SphP[i].MatrixT[7])+SphP[i].MatrixT[1]*(SphP[i].MatrixT[5]*SphP[i].MatrixT[6]-SphP[i].MatrixT[8]*SphP[i].MatrixT[3])+SphP[i].MatrixT[2]*(SphP[i].MatrixT[3]*SphP[i].MatrixT[7]-SphP[i].MatrixT[4]*SphP[i].MatrixT[6]);
   Tinv[0]=(SphP[i].MatrixT[4]*SphP[i].MatrixT[8]-SphP[i].MatrixT[5]*SphP[i].MatrixT[7])/fac;
   Tinv[1]=(SphP[i].MatrixT[2]*SphP[i].MatrixT[7]-SphP[i].MatrixT[1]*SphP[i].MatrixT[8])/fac;
   Tinv[2]=(SphP[i].MatrixT[1]*SphP[i].MatrixT[5]-SphP[i].MatrixT[2]*SphP[i].MatrixT[4])/fac;
   Tinv[3]=(SphP[i].MatrixT[5]*SphP[i].MatrixT[6]-SphP[i].MatrixT[3]*SphP[i].MatrixT[8])/fac;
   Tinv[4]=(SphP[i].MatrixT[0]*SphP[i].MatrixT[8]-SphP[i].MatrixT[2]*SphP[i].MatrixT[6])/fac;
   Tinv[5]=(SphP[i].MatrixT[2]*SphP[i].MatrixT[3]-SphP[i].MatrixT[0]*SphP[i].MatrixT[5])/fac;
   Tinv[6]=(SphP[i].MatrixT[3]*SphP[i].MatrixT[7]-SphP[i].MatrixT[4]*SphP[i].MatrixT[6])/fac;
   Tinv[7]=(SphP[i].MatrixT[6]*SphP[i].MatrixT[1]-SphP[i].MatrixT[0]*SphP[i].MatrixT[7])/fac;
   Tinv[8]=(SphP[i].MatrixT[0]*SphP[i].MatrixT[4]-SphP[i].MatrixT[1]*SphP[i].MatrixT[3])/fac;
   //Estimate the velocity matrix as D.T^-1
   V[0]=SphP[i].MatrixD[0]*Tinv[0]+SphP[i].MatrixD[1]*Tinv[3]+SphP[i].MatrixD[2]*Tinv[6];
   V[1]=SphP[i].MatrixD[0]*Tinv[1]+SphP[i].MatrixD[1]*Tinv[4]+SphP[i].MatrixD[2]*Tinv[7];
   V[2]=SphP[i].MatrixD[0]*Tinv[2]+SphP[i].MatrixD[1]*Tinv[5]+SphP[i].MatrixD[2]*Tinv[8];
   V[3]=SphP[i].MatrixD[3]*Tinv[0]+SphP[i].MatrixD[4]*Tinv[3]+SphP[i].MatrixD[5]*Tinv[6];
   V[4]=SphP[i].MatrixD[3]*Tinv[1]+SphP[i].MatrixD[4]*Tinv[4]+SphP[i].MatrixD[5]*Tinv[7];
   V[5]=SphP[i].MatrixD[3]*Tinv[2]+SphP[i].MatrixD[4]*Tinv[5]+SphP[i].MatrixD[5]*Tinv[8];
   V[6]=SphP[i].MatrixD[6]*Tinv[0]+SphP[i].MatrixD[7]*Tinv[3]+SphP[i].MatrixD[8]*Tinv[6];
   V[7]=SphP[i].MatrixD[6]*Tinv[1]+SphP[i].MatrixD[7]*Tinv[4]+SphP[i].MatrixD[8]*Tinv[7];
   V[8]=SphP[i].MatrixD[6]*Tinv[2]+SphP[i].MatrixD[7]*Tinv[5]+SphP[i].MatrixD[8]*Tinv[8];
   //Estimate divv as trace(V)
   divv=SphP[i].MatrixD[0]*Tinv[0]+SphP[i].MatrixD[1]*Tinv[3]+SphP[i].MatrixD[2]*Tinv[6]+SphP[i].MatrixD[3]*Tinv[1]+SphP[i].MatrixD[4]*Tinv[4]+SphP[i].MatrixD[5]*Tinv[7]+SphP[i].MatrixD[6]*Tinv[2]+SphP[i].MatrixD[7]*Tinv[5]+SphP[i].MatrixD[8]*Tinv[8];
   //Calculate the switch
   tmp=1.0/NUMDIMS;
   R_i = SphP[i].R_i / SphP[i].Density;
   trSSt = (V[0]-tmp*divv)*(V[0]-tmp*divv)+(V[4]-tmp*divv)*(V[4]-tmp*divv)+(V[8]-tmp*divv)*(V[8]-tmp*divv)+0.5*(V[1]+V[3])*(V[1]+V[3])+.5*(V[2]+V[6])*(V[2]+V[6])+.5*(V[5]+V[7])*(V[5]+V[7]);
   xi_i = pow(2.0*pow(1.0-R_i,4.0)*divv,2.0);
   xi_i = xi_i / (xi_i+trSSt);
   //Calculate the new indicator, using the old value for DivVel
   dt = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;
   divdotv=(divv-SphP[i].DivVelNew)/dt;
   if(divdotv>=0.0)
   {
     A_i = -1.0 *xi_i *divdotv;
     alpha_loc = All.AlphaMax * SphP[i].Hsml*SphP[i].Hsml*A_i / (SphP[i].SignalVel*SphP[i].SignalVel + SphP[i].Hsml*SphP[i].Hsml*A_i);
   }
   else
   {
     alpha_loc = 0.0;
   }
   //Not sure if it's OK to do the time integration here...
   if(alpha_loc > SphP[i].ArtVisc)
   {
     SphP[i].ArtVisc = alpha_loc;
   }
   else
   {
     tau_i = SphP[i].Hsml / (2.0 * SphP[i].SignalVel * All.ArtViscDecayLength);
     SphP[i].ArtVisc = alpha_loc - (alpha_loc - SphP[i].ArtVisc) * exp( -1.0 *dt / tau_i);
   }
   //Now you can calculate the true dEntropy
   SphP[i].DtEntropy *= SphP[i].ArtVisc *SphP[i].DhsmlDensityFactor / SphP[i].Density / pow(SphP[i].Hsml,NUMDIMS+2);
   SphP[i].DivVelNew=divv;
   SphP[i].DivDotVel=divdotv;
#endif

	SphP[i].DtEntropy *= GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].Density, GAMMA_MINUS1));
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
  double acc[3], dtEntropy, maxSignalVel,Di[9],Ti[9];
  double dx, dy, dz, dvx, dvy, dvz;
  double h_i2, h_i5,hinv, hinv3,hinv4;
  double p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;
  double hfc, wk_i, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;
  double h_j, dwk_j, r, r2, u, hfc_visc;
#ifdef ARTVISCTEST
  double oldArtVisc[3];
#endif
#ifdef INDIVIDUALAV
  double vsig_i,R_i,tmp,fac,mass_j;
#endif

#ifndef NOVISCOSITYLIMITER
  double dt;
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
      f1 = fabs(SphP[target].DivVel) /
	(fabs(SphP[target].DivVel) + SphP[target].CurlVel +
	 0.0001 * soundspeed_i / SphP[target].Hsml / fac_mu);
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
    }


  /* initialize variables before SPH loop is started */
  acc[0] = acc[1] = acc[2] = dtEntropy = 0;
  Di[0]=Di[1]=Di[2]=Di[3]=Di[4]=Di[5]=Di[6]=Di[7]=Di[8]=Di[9]=0;
  Ti[0]=Ti[1]=Ti[2]=Ti[3]=Ti[4]=Ti[5]=Ti[6]=Ti[7]=Ti[8]=Ti[9]=0;
  R_i=vsig_i=0.0;
  maxSignalVel = 0;

  p_over_rho2_i = pressure / (rho * rho) * dhsmlDensityFactor;
  h_i2 = h_i * h_i;
  h_i5 = pow(h_i,NUMDIMS+2);

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
            hinv3 = hinv * hinv * hinv;
#else
		      hinv4 = hinv * hinv * hinv / boxSize_Z;
            hinv3 = hinv * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      if(u < 0.5)
            {
#ifdef INDIVIDUALAV
              wk_i = hinv3 * (KERNEL_COEFF_1 * KERNEL_COEFF_2 * (u-1)*u*u);
#endif
	       	  dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
            }
		      else
            {
#ifdef INDIVIDUALAV
              wk_i = hinv3*(KERNEL_COEFF_5 * (1.0-u)*(1.0-u)*(1.0-u));
#endif
			     dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
            }
		    }
		  else
		    {
#ifdef INDIVIDUALAV
            wk_i = 0;
#endif
		      dwk_i = 0;
		    }

		  if(r2 < h_j * h_j)
		    {
		      hinv = 1.0 / h_j;
#ifndef  TWODIMS
		      hinv4 = hinv * hinv * hinv * hinv;
#else
		      hinv4 = hinv * hinv * hinv / boxSize_Z;
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

#ifdef INDIVIDUALAV
        mass_j=P[j].Mass;
        //Artificial viscosity stuff..
        fac = h_i5*mass_j * dwk_i / r / SphP[j].Density;
        Di[0] += fac * (dvx*dx);
        Di[1] += fac * (dvx*dy);
        Di[2] += fac * (dvx*dz);
        Di[3] += fac * (dvy*dx);
        Di[4] += fac * (dvy*dy);
        Di[5] += fac * (dvy*dz);
        Di[6] += fac * (dvz*dx);
        Di[7] += fac * (dvz*dy);
        Di[8] += fac * (dvz*dz);

        Ti[0] += fac * (dx*dx);
        Ti[1] += fac * (dx*dy);
        Ti[2] += fac * (dx*dz);
        Ti[3] += fac * (dy*dx);
        Ti[4] += fac * (dy*dy);
        Ti[5] += fac * (dy*dz);
        Ti[6] += fac * (dz*dx);
        Ti[7] += fac * (dz*dy);
        Ti[8] += fac * (dz*dz);
        //Switch
        if(SphP[j].DivVel>=0)
          R_i += mass_j*wk_i;
        else
          R_i -= mass_j*wk_i;
        //The signal velocity....
        tmp=0.5 * (soundspeed_i + soundspeed_j);
        if(vdotr2<=0)
          tmp -=vdotr2;
        if(tmp>vsig_i)
          vsig_i=tmp;
        //Calculate all the j dependent factors
        if(vdotr2 <0)
        {
          mu_ij = vdotr2 / ((h_i+h_j)*r2*((1.0/h_i)*(1.0/h_i)+(1.0/h_j)*(1.0/h_j)));
          mu_ij = -1.0 * mu_ij * (0.5*(soundspeed_i+soundspeed_j)-All.ArtViscbparam*mu_ij);
          dtEntropy += mass_j*vdotr2*mu_ij*0.5*dwk_i*h_i5/r;
        }

#endif

#ifdef ARTVISCTEST
		  if(vdotr2 < 0)	/* ... artificial viscosity */
		    {
		      mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */

                      vsig = soundspeed_i + soundspeed_j - 3 * mu_ij;

		      if(vsig > maxSignalVel)
			maxSignalVel = vsig;

		      rho_ij = 0.5 * (rho + SphP[j].Density);
		      f2 =
			fabs(SphP[j].DivVel) / (fabs(SphP[j].DivVel) + SphP[j].CurlVel +
						0.0001 * soundspeed_j / fac_mu / SphP[j].Hsml);


		      visc = 0.25 * All.ArtBulkViscConst * vsig * (-mu_ij) / rho_ij * (f1 + f2);

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
#else
        visc=0;
#endif

		  p_over_rho2_j *= SphP[j].DhsmlDensityFactor;

		  hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;
#ifdef ARTVISCTEST
        oldArtVisc[0] -= hfc_visc *dx;
        oldArtVisc[1] -= hfc_visc *dy;
        oldArtVisc[2] -= hfc_visc *dz;
#endif

        //Only include this artificial viscosity in the hydro acceleration if we're not calculating the more sophisticated one
        //If we're doing the per particle art visc, that will be added to the acceleration in a separate loop
#ifndef INDIVIDUALAV
		  hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;
		  dtEntropy += 0.5 * hfc_visc * vdotr2;
#else
        hfc = P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;
#endif

		  acc[0] -= hfc * dx;
		  acc[1] -= hfc * dy;
		  acc[2] -= hfc * dz;
		}
	    }
	}
    }
  while(startnode >= 0);



  /* Now collect the result at the right place */
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
      {
	     SphP[target].HydroAccel[k] = acc[k];
#ifdef ARTVISCTEST
        SphP[target].OldArtViscAccel[k] = oldArtVisc[k];
#endif
      }
#ifdef INDIVIDUALAV
      for(k=0; k<9; k++)
      {
        SphP[target].MatrixD[k] = Di[k];
        SphP[target].MatrixT[k] = Ti[k];
      }
      SphP[target].R_i = R_i;
      SphP[target].SignalVel = vsig_i;
#endif
      SphP[target].DtEntropy = dtEntropy;
      SphP[target].MaxSignalVel = maxSignalVel;
    }
  else
    {
      for(k = 0; k < 3; k++)
      {
	     HydroDataResult[target].Acc[k] = acc[k];
#ifdef ARTVISCTEST
        HydroDataResult[target].OldArtViscAccel[k] = oldArtVisc[k];
#endif
      }
#ifdef INDIVIDUALAV
      for(k=0; k<9; k++)
      {
        HydroDataResult[target].MatrixD[k] = Di[k];
        HydroDataResult[target].MatrixT[k] = Ti[k];
      }
      HydroDataResult[target].R_i = R_i;
      HydroDataResult[target].SignalVel = vsig_i;
#endif
      HydroDataResult[target].DtEntropy = dtEntropy;
      HydroDataResult[target].MaxSignalVel = maxSignalVel;
    }
}



/*! This function calculates the acceleration due to artificial viscosity*/
void hydro_visc_force(void)
{
  long long ntot, ntotleft;
  int i, j, k, n, ngrp, maxfill, source, ndone;
  int *nbuffer, *noffset, *nsend_local, *nsend, *numlist, *ndonelist;
  int level, sendTask, recvTask, nexport, place;
  double tstart, tend, sumt, sumcomm;
  double timecomp = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  MPI_Status status;

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

  //Reasonably sure you could malloc this to be just NTask * sizeof(int)
  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


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

	    visc_evaluate(i, 0);

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
		    HydroDataIn[nexport].DhsmlDensityFactor = SphP[i].DhsmlDensityFactor;
		    HydroDataIn[nexport].Density = SphP[i].Density;
          HydroDataIn[nexport].ArtVisc = SphP[i].ArtVisc;
		    HydroDataIn[nexport].Pressure = SphP[i].Pressure;

		    HydroDataIn[nexport].Index = i;
		    HydroDataIn[nexport].Task = j;
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
	    visc_evaluate(j, 1);
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
           {
			    SphP[place].HydroAccel[k] += HydroDataPartialResult[source].Acc[k];
			    SphP[place].ArtViscAccel[k] += HydroDataPartialResult[source].Acc[k];
           }
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

#ifdef SPH_BND_PARTICLES
  for(i = 0; i < N_gas; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
    {
      if(P[i].ID == 0)
	   {
	     for(k = 0; k < 3; k++)
	       SphP[i].HydroAccel[k] = 0;
	   }
    }
#endif

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

/* The heart of the loop that calculates the acceleration due to artificial viscosity */
void visc_evaluate(int target, int mode)
{
  FLOAT *pos, *vel;
  FLOAT h_i,f_i, rho, alpha_i, pressure;
  int startnode, numngb,n,j;
  double soundspeed_i,p_over_rho2_i,h_i2,dx,dy,dz,r2,r,h_j,dvx,dvy,dvz,vdotr;
  double p_over_rho2_j,soundspeed_j,hinv,hinv4,u,dwk_i,dwk_j,hinv_j,mu_ij,visc;
  double acc[3];

  if(mode == 0)
  {
    pos = P[target].Pos;
    vel = SphP[target].VelPred;
    h_i = SphP[target].Hsml;
    f_i = SphP[target].DhsmlDensityFactor;
    rho = SphP[target].Density;
    alpha_i = SphP[target].ArtVisc;
    pressure = SphP[target].Pressure;
  }
  else
  {
    pos = HydroDataGet[target].Pos;
    vel = HydroDataGet[target].Vel;
    h_i = HydroDataGet[target].Hsml;
    f_i = HydroDataGet[target].DhsmlDensityFactor;
    rho = HydroDataGet[target].Density;
    alpha_i = HydroDataGet[target].ArtVisc;
    pressure = HydroDataGet[target].Pressure;
  }
  soundspeed_i = sqrt(GAMMA * pressure / rho);


  /* initialize variables before SPH loop is started */
  acc[0] = acc[1] = acc[2] ;
  p_over_rho2_i = pressure / (rho * rho) * f_i;
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
	       dvx = vel[0] - SphP[j].VelPred[0];
	       dvy = vel[1] - SphP[j].VelPred[1];
	       dvz = vel[2] - SphP[j].VelPred[2];
		    vdotr = dx * dvx + dy * dvy + dz * dvz;
          if(All.ComovingIntegrationOn)
            vdotr += hubble_a2 *r2;
          if(vdotr<0)
          {
		      p_over_rho2_j = SphP[j].Pressure / (SphP[j].Density * SphP[j].Density);
		      soundspeed_j = sqrt(GAMMA * p_over_rho2_j * SphP[j].Density);
            hinv=1.0/h_i;
            hinv_j=1.0/h_j;

		      if(r2 < h_i2)
		      {
#ifndef  TWODIMS
		        hinv4 = hinv * hinv * hinv * hinv;
#else
		        hinv4 = hinv * hinv * hinv / boxSize_Z;
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
#ifndef  TWODIMS
		        hinv4 = hinv_j * hinv_j * hinv_j * hinv_j;
#else
		        hinv4 = hinv_j * hinv_j * hinv_j / boxSize_Z;
#endif
		        u = r * hinv_j;
		        if(u < 0.5)
			       dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		        else
			       dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		      }
		      else
		      {
		        dwk_j = 0;
		      }
            mu_ij = vdotr / ((h_i+h_j)*r2*(hinv*hinv+hinv_j*hinv_j));
            mu_ij = -1.0 * mu_ij * (0.5*(soundspeed_i+soundspeed_j)-All.ArtViscbparam*mu_ij);
            visc = mu_ij *P[j].Mass*(((alpha_i*f_i*dwk_i)/(rho))+((SphP[j].ArtVisc*SphP[j].DhsmlDensityFactor*dwk_j)/(SphP[j].Density)))/r;
            acc[0] -= visc*dx;
            acc[1] -= visc*dy;
            acc[2] -= visc*dz;
          }
 		  }
      }
    }
  }
  while(startnode >= 0);



  /* Now collect the result at the right place */
  if(mode == 0)
  {
    for(j = 0; j < 3; j++)
    {
      SphP[target].ArtViscAccel[j] = acc[j];
      SphP[target].HydroAccel[j] += acc[j];
    }
  }
  else
  {
    for(j = 0; j < 3; j++)
    {
	   HydroDataResult[target].Acc[j] = acc[j];
    }
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
