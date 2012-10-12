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

#if defined BETA_COOLING || defined NK_AV
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
#endif
#if defined BETA_COOLING || defined NK_AV
  int numsinks,root,globalroot;
#endif
#ifdef CDAV_DRIFTUPDATE
  double divv,diva,xi,A,fac;
  double Tinv[9],V[9];
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

#if defined BETA_COOLING || defined NK_AV
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
#ifdef CDAV_DRIFTUPDATE
          HydroDataIn[nexport].oldAccel[0] = SphP[i].oldAccel[0];
          HydroDataIn[nexport].oldAccel[1] = SphP[i].oldAccel[1];
          HydroDataIn[nexport].oldAccel[2] = SphP[i].oldAccel[2];
          HydroDataIn[nexport].divvsign = (SphP[i].DivVel >0) - (SphP[i].DivVel<0);
          HydroDataIn[nexport].gradRho[0] = SphP[i].gradRho[0];
          HydroDataIn[nexport].gradRho[1] = SphP[i].gradRho[1];
          HydroDataIn[nexport].gradRho[2] = SphP[i].gradRho[2];
#endif
#ifdef VAR_H_TEST
          HydroDataIn[nexport].htest_f[0] = SphP[i].htest_f[0];
          HydroDataIn[nexport].htest_f[1] = SphP[i].htest_f[1];
          HydroDataIn[nexport].htest_f[2] = SphP[i].htest_f[2];
          HydroDataIn[nexport].htest_g = SphP[i].htest_g;
#endif

		    HydroDataIn[nexport].Index = i;
		    HydroDataIn[nexport].Task = j;
        
#if defined MMAV_DRIFTUPDATE || defined MMAV || defined CDAV || defined CDAV_DRIFTUPDATE
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
#ifdef NK_AV
           SphP[place].NumNK += HydroDataPartialResult[source].NumNK;
           SphP[place].NumN += HydroDataPartialResult[source].NumN;
#endif
#ifdef CDAV_DRIFTUPDATE
           for(k=0;k<9;k++)
           {
             SphP[place].T[k]+=HydroDataPartialResult[source].T[k];
             SphP[place].D[k]+=HydroDataPartialResult[source].D[k];
             SphP[place].E[k]+=HydroDataPartialResult[source].E[k];
           }
           SphP[place].R+=HydroDataPartialResult[source].R;
           if(HydroDataPartialResult[source].CDAVSignalVel > SphP[place].CDAVSignalVel)
           {
             SphP[place].CDAVSignalVel = HydroDataPartialResult[source].CDAVSignalVel;
           }
#endif
#ifdef VAR_H_TEST
           SphP[place].htest_t += HydroDataPartialResult[source].htest_t;
#endif
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
#endif
#ifdef SPH_BND_PARTICLES
	if(P[i].ID == 0)
	  {
	    SphP[i].DtEntropy = 0;
	    for(k = 0; k < 3; k++)
	      SphP[i].HydroAccel[k] = 0;
	  }
#endif
#ifdef CDAV_DRIFTUPDATE
   SphP[i].R /= SphP[i].Density;
   //Invert T
   fac = 1/(SphP[i].T[0]*(SphP[i].T[4]*SphP[i].T[8]-SphP[i].T[5]*SphP[i].T[7])-SphP[i].T[1]*(SphP[i].T[8]*SphP[i].T[3]-SphP[i].T[5]*SphP[i].T[6])+SphP[i].T[2]*(SphP[i].T[3]*SphP[i].T[7]-SphP[i].T[4]*SphP[i].T[6]));
   Tinv[0]=fac*(SphP[i].T[4]*SphP[i].T[8]-SphP[i].T[5]*SphP[i].T[7]);
   Tinv[1]=fac*(SphP[i].T[5]*SphP[i].T[6]-SphP[i].T[3]*SphP[i].T[8]);
   Tinv[2]=fac*(SphP[i].T[3]*SphP[i].T[7]-SphP[i].T[4]*SphP[i].T[6]);
   Tinv[3]=fac*(SphP[i].T[2]*SphP[i].T[7]-SphP[i].T[1]*SphP[i].T[8]);
   Tinv[4]=fac*(SphP[i].T[0]*SphP[i].T[8]-SphP[i].T[2]*SphP[i].T[6]);
   Tinv[5]=fac*(SphP[i].T[6]*SphP[i].T[1]-SphP[i].T[0]*SphP[i].T[7]);
   Tinv[6]=fac*(SphP[i].T[1]*SphP[i].T[5]-SphP[i].T[2]*SphP[i].T[4]);
   Tinv[7]=fac*(SphP[i].T[2]*SphP[i].T[3]-SphP[i].T[0]*SphP[i].T[5]);
   Tinv[8]=fac*(SphP[i].T[0]*SphP[i].T[4]-SphP[i].T[1]*SphP[i].T[3]);
   //Multiply D by Tinv to calculate V
   V[0]=SphP[i].D[0]*Tinv[0]+SphP[i].D[1]*Tinv[3]+SphP[i].D[2]*Tinv[6];
   V[3]=SphP[i].D[0]*Tinv[1]+SphP[i].D[1]*Tinv[4]+SphP[i].D[2]*Tinv[7];
   V[6]=SphP[i].D[0]*Tinv[2]+SphP[i].D[1]*Tinv[5]+SphP[i].D[2]*Tinv[8];
   V[1]=SphP[i].D[3]*Tinv[0]+SphP[i].D[4]*Tinv[3]+SphP[i].D[5]*Tinv[6];
   V[4]=SphP[i].D[3]*Tinv[1]+SphP[i].D[4]*Tinv[4]+SphP[i].D[5]*Tinv[7];
   V[7]=SphP[i].D[3]*Tinv[2]+SphP[i].D[4]*Tinv[5]+SphP[i].D[5]*Tinv[8];
   V[2]=SphP[i].D[6]*Tinv[0]+SphP[i].D[7]*Tinv[3]+SphP[i].D[8]*Tinv[6];
   V[5]=SphP[i].D[6]*Tinv[1]+SphP[i].D[7]*Tinv[4]+SphP[i].D[8]*Tinv[7];
   V[8]=SphP[i].D[6]*Tinv[2]+SphP[i].D[7]*Tinv[5]+SphP[i].D[8]*Tinv[8];
   //Now we can calculate DivVel
   divv = V[0] + V[4] + V[8];
   if(ThisTask==1)
   {
     //printf("The ratio of the old to the new divv is %g for pcl %d\n",SphP[i].DivVel/divv,i);
   }
   //Next, calculate Div acceleration = tr(E.T^-1)-tr(V^2)
   diva = SphP[i].E[0]*Tinv[0]+SphP[i].E[1]*Tinv[3]+SphP[i].E[2]*Tinv[6] + 
     SphP[i].E[3]*Tinv[1]+SphP[i].E[4]*Tinv[4]+SphP[i].E[5]*Tinv[7] + 
     SphP[i].E[6]*Tinv[2]+SphP[i].E[7]*Tinv[5]+SphP[i].E[8]*Tinv[8];
   //Subtract the V^2 part
   diva = diva-(V[0]*V[0]+V[4]*V[4]+V[8]*V[8]+2*(V[1]*V[3]+V[2]*V[6]+V[5]*V[7]));
   //Alternative way of calculating it...
   //dt_alpha = (P[i].Ti_endstep-P[i].Ti_begstep)*All.Timebase_interval;
   //if(SphP[i].DtDrift==0 || SphP[i].oldDivVel==0)
   //{
   //  diva=0;
   //}
   //else
   //{
   //  diva=((divv-SphP[i].oldDivVel)/SphP[i].DtDrift);
   //}
   //SphP[i].DivVel=divv;
   SphP[i].oldDivVel=divv;
   //Calculate xi
   xi=2*(1-SphP[i].R)*(1-SphP[i].R)*(1-SphP[i].R)*(1-SphP[i].R)*divv;
   xi *= xi;
   //The added bit is tr(S.S^T)
   //A here is a temporary variable
   A=V[0]*V[0]+V[4]*V[4]+V[8]*V[8] + 
        0.5*((V[1]+V[3])*(V[1]+V[3])+(V[2]+V[6])*(V[2]+V[6])+(V[5]+V[7])*(V[5]+V[7])) -
        divv*divv/3;
   if(xi+A!=0)
   {
     xi /=xi+A;
   }
#ifdef NOBALSARA
   xi=1;
#endif
   //Now calculate A_i
   A=xi*dmax(-diva,0);
   //Calculate alpha_loc, store it in DtAlpha
   SphP[i].DtAlpha = SphP[i].CDAVSignalVel*SphP[i].CDAVSignalVel + SphP[i].Hsml*SphP[i].Hsml*A;
   if(SphP[i].DtAlpha!=0)
   {
     SphP[i].DtAlpha = All.ArtBulkViscConst*SphP[i].Hsml*SphP[i].Hsml*A/SphP[i].DtAlpha;
   }
#endif
#ifdef VAR_H_TEST
   //Multiply by the final factors
   SphP[i].htest_t *= SphP[i].Hsml*SphP[i].DhsmlDensityFactor/(NUMDIMS*SphP[i].Density*SphP[i].Density);
   if(ThisTask==1)
   {
     printf("The relative contribution of the extra term (term/div.v) to div.v is %g\n",SphP[i].htest_t/SphP[i].DivVel);
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
  FLOAT mass, h_i, dhsmlDensityFactor, rho, pressure, f1, alpha_visc;
#if defined MMAV_DRIFTUPDATE || defined MMAV || defined CDAV || defined CDAV_DRIFTUPDATE
  FLOAT alpha_visc_j;
#else
  FLOAT f2;
  alpha_visc = All.ArtBulkViscConst;
#endif
#ifdef CDAV_DRIFTUPDATE
  FLOAT *oldacc,*gradRho;
  FLOAT divvsign;
  double fac_1,fac_2,R,T[9],D[9],E[9],CDAVvsig,wk_i,mass_j,dax,day,daz;
#endif
#ifdef VAR_H_TEST
  double wk_i,rinv,htest_term[3],htest_t;
  FLOAT *htest_f;
  FLOAT htest_g;
#endif

  double acc[3], dtEntropy, maxSignalVel;
  double dx, dy, dz, dvx, dvy, dvz;
  double h_i2, hinv, hinv4;
  double p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;
  double hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;
  double h_j, dwk_j, r, r2, u, hfc_visc;
  double tmp;
#ifndef NOVISCOSITYLIMITER
  double dt;
#endif
#ifdef NK_AV
  double v2r_i,v2r_j,r2_j,v2_j,NK_test;
  int numN,numNK;
#endif
#ifdef PRICE_GRAV_SOFT
  double zeta=0;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
#ifdef CDAV_DRIFTUPDATE
      oldacc = SphP[target].oldAccel;
      divvsign = (SphP[target].DivVel>0) - (SphP[target].DivVel<0);
      gradRho = SphP[target].gradRho;
#endif
      h_i = SphP[target].Hsml;
      mass = P[target].Mass;
      dhsmlDensityFactor = SphP[target].DhsmlDensityFactor;
      rho = SphP[target].Density;
      pressure = SphP[target].Pressure;
      timestep = P[target].Ti_endstep - P[target].Ti_begstep;
      soundspeed_i = sqrt(GAMMA * pressure / rho);
#if defined MMAV_DRIFTUPDATE || defined MMAV || defined CDAV || defined CDAV_DRIFTUPDATE
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
#ifdef VAR_H_TEST
      htest_f = SphP[target].htest_f;
      htest_g = SphP[target].htest_g;
#endif
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
#ifdef CDAV_DRIFTUPDATE
      oldacc = HydroDataGet[target].oldAccel;
      divvsign = HydroDataGet[target].divvsign;
      gradRho = HydroDataGet[target].gradRho;
#endif
      h_i = HydroDataGet[target].Hsml;
      mass = HydroDataGet[target].Mass;
      dhsmlDensityFactor = HydroDataGet[target].DhsmlDensityFactor;
      rho = HydroDataGet[target].Density;
      pressure = HydroDataGet[target].Pressure;
      timestep = HydroDataGet[target].Timestep;
      soundspeed_i = sqrt(GAMMA * pressure / rho);
      f1 = HydroDataGet[target].F1;
#if defined MMAV_DRIFTUPDATE || defined MMAV || defined CDAV || defined CDAV_DRIFTUPDATE
      alpha_visc = HydroDataGet[target].Alpha;
#endif
#ifdef PRICE_GRAV_SOFT
      zeta = HydroDataGet[target].Zeta;
#endif
#ifdef VAR_H_TEST
      htest_f = HydroDataGet[target].htest_f;
      htest_g = HydroDataGet[target].htest_g;
#endif
    }


  /* initialize variables before SPH loop is started */
  acc[0] = acc[1] = acc[2] = dtEntropy = 0;
  maxSignalVel = 0;
#ifdef CDAV
  //Because the r=0 case will not be evaluated below...
  maxSignalVel = soundspeed_i;
#endif

  p_over_rho2_i = pressure / (rho * rho) * dhsmlDensityFactor;
#ifdef PRICE_GRAV_SOFT
  p_over_rho2_i = dhsmlDensityFactor * ((pressure/(rho*rho))+(zeta*All.G*.5));
#endif
  h_i2 = h_i * h_i;
  
#ifdef NK_AV
  //Calculate the local values
  v2r_i = ((starData[4]-vel[0])*(starData[4]-vel[0])+(starData[5]-vel[1])*(starData[5]-vel[1])+(starData[6]-vel[2])*(starData[6]-vel[2]))*sqrt((starData[0]-pos[0])*(starData[0]-pos[0])+(starData[1]-pos[1])*(starData[1]-pos[1])+(starData[2]-pos[2])*(starData[2]-pos[2]));
  numN=numNK=0;
#endif

#ifdef CDAV_DRIFTUPDATE
  CDAVvsig=0;
  //Need to include the self-reference (i.e. b=a) in this calculation
  R=divvsign*mass*KERNEL_COEFF_1/(h_i*h_i*h_i);
  for(k=0;k<9;k++)
  {
    E[k]=D[k]=T[k]=0;
  }
#endif

#ifdef VAR_H_TEST
  htest_t=0;
#endif

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

#if defined MMAV_DRIFTUPDATE || defined MMAV || defined CDAV || defined CDAV_DRIFTUPDATE
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
#ifdef NK_AV
            numN++;
            v2_j = ((starData[4]-SphP[j].VelPred[0])*(starData[4]-SphP[j].VelPred[0])+(starData[5]-SphP[j].VelPred[1])*(starData[5]-SphP[j].VelPred[1])+(starData[6]-SphP[j].VelPred[2])*(starData[6]-SphP[j].VelPred[2]));
            r2_j = ((starData[0]-P[j].Pos[0])*(starData[0]-P[j].Pos[0])+(starData[1]-P[j].Pos[1])*(starData[1]-P[j].Pos[1])+(starData[2]-P[j].Pos[2])*(starData[2]-P[j].Pos[2]));
            NK_test =fabs(((SphP[j].VelPred[0]-starData[4])*(P[j].Pos[0]-starData[0])+(SphP[j].VelPred[1]-starData[5])*(P[j].Pos[1]-starData[1])+(SphP[j].VelPred[2]-starData[6])*(P[j].Pos[2]-starData[2]))/sqrt(v2_j*r2_j));
            v2r_j = v2_j * sqrt(r2_j);
            //if(NK_test > All.NKtollerence || fabs(1-(v2r_i/v2r_j)) < All.NKtollerence)
            if(fabs(1-(v2r_i/v2r_j)) < All.NKtollerence)
            {
              numNK++;
              //printf("The particle had ratio %g and dot prod %g\n",NK_test,1-(v2r_i/v2r_j));
            }
#endif
		      hinv = 1.0 / h_i;
#ifndef  TWODIMS
		      hinv4 = hinv * hinv * hinv * hinv;
#else
		      hinv4 = hinv * hinv * hinv / boxSize_Z;

#endif
		      u = r * hinv;
		      if(u < 0.5)
            {
			     dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
#if defined CDAV_DRIFTUPDATE || defined VAR_H_TEST
              wk_i = hinv*hinv*hinv*(KERNEL_COEFF_1+KERNEL_COEFF_2 * (u-1)*u*u);
#endif
            }
		      else
            {
			     dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
#if defined CDAV_DRIFTUPDATE || defined VAR_H_TEST
              wk_i = hinv*hinv*hinv*KERNEL_COEFF_5*(1.0 -u)*(1.0-u)*(1.0-u);
#endif
            }
		    }
		  else
		    {
		      dwk_i = 0;
#if defined VAR_H_TEST || defined CDAV_DRIFTUPDATE 
            wk_i = 0;
#endif
		    }
#ifdef CDAV_DRIFTUPDATE
        mass_j = P[j].Mass;
        hinv=1/h_i;
        u=r*hinv;
        fac_1 = mass_j * dwk_i / r;
        fac_2 = mass_j * (1/rho)* (((dwk_i *r*h_i)/NUMDIMS)+wk_i);
        dax = acc[0] - SphP[j].oldAccel[0];
        day = acc[1] - SphP[j].oldAccel[1];
        daz = acc[2] - SphP[j].oldAccel[2];

        D[0] += dvx * (fac_1 * dx + fac_2 * gradRho[0]);
        D[1] += dvx * (fac_1 * dy + fac_2 * gradRho[1]);
        D[2] += dvx * (fac_1 * dz + fac_2 * gradRho[2]);
        D[3] += dvy * (fac_1 * dx + fac_2 * gradRho[0]);
        D[4] += dvy * (fac_1 * dy + fac_2 * gradRho[1]);
        D[5] += dvy * (fac_1 * dz + fac_2 * gradRho[2]);
        D[6] += dvz * (fac_1 * dx + fac_2 * gradRho[0]);
        D[7] += dvz * (fac_1 * dy + fac_2 * gradRho[1]);
        D[8] += dvz * (fac_1 * dz + fac_2 * gradRho[2]);

        E[0] += dax * (fac_1 * dx  +fac_2 * gradRho[0]);
        E[1] += dax * (fac_1 * dy  +fac_2 * gradRho[1]);
        E[2] += dax * (fac_1 * dz  +fac_2 * gradRho[2]);
        E[3] += day * (fac_1 * dx  +fac_2 * gradRho[0]);
        E[4] += day * (fac_1 * dy + fac_2 * gradRho[1]);
        E[5] += day * (fac_1 * dz + fac_2 * gradRho[2]);
        E[6] += daz * (fac_1 * dx + fac_2 * gradRho[0]);
        E[7] += daz * (fac_1 * dy + fac_2 * gradRho[1]);
        E[8] += daz * (fac_1 * dz + fac_2 * gradRho[2]);

        T[0] += dx * (fac_1 * dx + fac_2 * gradRho[0]);
        T[1] += dx * (fac_1 * dy + fac_2 * gradRho[1]);
        T[2] += dx * (fac_1 * dz + fac_2 * gradRho[2]);
        T[3] += dy * (fac_1 * dx + fac_2 * gradRho[0]);
        T[4] += dy * (fac_1 * dy + fac_2 * gradRho[1]);
        T[5] += dy * (fac_1 * dz + fac_2 * gradRho[2]);
        T[6] += dz * (fac_1 * dx + fac_2 * gradRho[0]);
        T[7] += dz * (fac_1 * dy + fac_2 * gradRho[1]);
        T[8] += dz * (fac_1 * dz + fac_2 * gradRho[2]);

        //This "divv" is the only part of the whole business which is not exact...
        //This misses out on the r=0 contribution...
        if(P[j].Type==0)
        {
          divvsign=0;
          if(SphP[j].DivVel>0)
          {
            divvsign=1;
          }
          else
          {
            if(SphP[j].DivVel<0)
            {
              divvsign=-1;
            }
          }
        }
        R += divvsign * mass_j * wk_i;
        //Calculate the signal velocity, only want to do this for things local to i...
        if(dwk_i!=0)
        {
          tmp=0.5*(soundspeed_i+soundspeed_j)-(1/r)*dmin(0,dx*dvx+dy*dvy+dz*dvz);
          if(tmp > CDAVvsig)
          {
            CDAVvsig=tmp;
          }
        }
#endif
#ifdef VAR_H_TEST
        rinv=1/r;
        hinv=1/h_i;
        u=r*hinv;
        htest_term[0]=-htest_f[0]*((NUMDIMS*wk_i*hinv)+htest_g*dwk_i*dx*rinv);
        htest_term[1]=-htest_f[1]*((NUMDIMS*wk_i*hinv)+htest_g*dwk_i*dy*rinv);
        htest_term[2]=-htest_f[2]*((NUMDIMS*wk_i*hinv)+htest_g*dwk_i*dz*rinv);
        htest_t -= P[j].Mass*htest_g*dwk_i*rinv*(dvx*dx+dvy*dy+dvz*dz);
        htest_t -= P[j].Mass*(dvx*htest_f[0]+dvy*htest_f[1]+dvz*htest_f[2])*(NUMDIMS*wk_i*hinv+dwk_i*u);
        //htest_t += P[j].Mass*(htest_term[0]*dvx + htest_term[1]*dvy+htest_term[2]*dvz);
        if(ThisTask==1)
        {
          //printf("Term_ab ~ (%g,%g,%g)\n",htest_term[0],htest_term[1],htest_term[2]);
        }
#endif

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

        //Each AV method is given a separate block.  At the end of the block hfc_visc should contain the magnitude of the viscous acceleration (divided by r) and dtEntropy should be updated.
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
#if defined MMAV || defined MMAV_DRIFTUPDATE || defined CDAV || defined CDAV_DRIFTUPDATE
        if(soundspeed_i+soundspeed_j >maxSignalVel)
          maxSignalVel = soundspeed_i+soundspeed_j;


        if(vdotr2 <0)
        {
		    rho_ij = 0.5 * (rho + SphP[j].Density);
          mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */
          //ArtViscPropConst is 3/2 in original implementation...
          vsig = soundspeed_i + soundspeed_j - All.ArtViscPropConst*2.0 * mu_ij;

	       if(vsig > maxSignalVel)
	         maxSignalVel = vsig;


          visc = 0.25 * (alpha_visc + alpha_visc_j) * vsig * (-mu_ij) / rho_ij;
#ifndef NOVISCOSITYLIMITER
          /* make sure that viscous acceleration is not too large */
		    dt = imax(timestep, (P[j].Ti_endstep - P[j].Ti_begstep)) * All.Timebase_interval;
		    if(dt > 0 && (dwk_i + dwk_j) < 0)
          {
            //Limiter is designed to set the magnitude of the viscous acceleration to be no more than a quarter of the line of sight velocity connecting i and j, divided by dt
            tmp=0.5*fac_vsic_fix * vdotr2 / (0.5 * (mass+P[j].Mass)*(dwk_i+dwk_j)*r*dt);
            if(visc>tmp)
            {
              visc = tmp;
            }
          }
#endif
        }
        else
        {
          visc=0;
        }
		  hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;
        dtEntropy += 0.5 * hfc_visc*vdotr2;
#endif
//Block for the standard method
#if !defined CDAV && !defined CDAV_DRIFTUPDATE && !defined MMAV && !defined MMAV_DRIFTUPDATE
        if(soundspeed_i+soundspeed_j >maxSignalVel)
          maxSignalVel = soundspeed_i+soundspeed_j;

		  if(vdotr2 < 0)	/* ... artificial viscosity */
		    {
		      mu_ij = fac_mu * vdotr2 / r;	/* note: this is negative! */
            //ArtViscPropConst is 3/2 in original implementation...
            vsig = soundspeed_i + soundspeed_j - All.ArtViscPropConst*2.0 * mu_ij;

		      if(vsig > maxSignalVel)
			maxSignalVel = vsig;

		      rho_ij = 0.5 * (rho + SphP[j].Density);
            //The standard thing...
#ifdef NOBALSARA
            f2=1.0;
#else
            f2 =
          fabs(SphP[j].DivVel) / (fabs(SphP[j].DivVel) + SphP[j].CurlVel +
                                  0.0001 * soundspeed_j / fac_mu / SphP[j].Hsml);
#endif
            visc = 0.25 * alpha_visc * vsig * (-mu_ij) / rho_ij * (f1 + f2);
		      /* .... end artificial viscosity evaluation */
#ifndef NOVISCOSITYLIMITER
		      /* make sure that viscous acceleration is not too large */
		      dt = imax(timestep, (P[j].Ti_endstep - P[j].Ti_begstep)) * All.Timebase_interval;
		      if(dt > 0 && (dwk_i + dwk_j) < 0)
			{
           //Limiter is designed to set the magnitude of the viscous acceleration to be no more than a quarter of the line of sight velocity connecting i and j, divided by dt
			  visc = dmin(visc, 0.5 * fac_vsic_fix * vdotr2 /
				      (0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
			}
#endif
//End no viscosity limiter
		    }
		  else
        {
		    visc = 0;
        }
		  hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;
        dtEntropy += 0.5 * hfc_visc * vdotr2;
#endif
//End AV blocks.  Entropy updated, need to update acceleration after adding pressure terms
		  p_over_rho2_j *= SphP[j].DhsmlDensityFactor;
		  hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;
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
	     SphP[target].HydroAccel[k] = acc[k];
      SphP[target].DtEntropy = dtEntropy;
      SphP[target].MaxSignalVel = maxSignalVel;
#ifdef NK_AV
      SphP[target].NumNK = numNK;
      SphP[target].NumN = numN;
#endif
#ifdef CDAV_DRIFTUPDATE
      for(k=0;k<9;k++)
      {
        SphP[target].T[k]=T[k];
        SphP[target].D[k]=D[k];
        SphP[target].E[k]=E[k];
      }
      SphP[target].R=R;
      SphP[target].CDAVSignalVel=CDAVvsig;
#endif
#ifdef VAR_H_TEST
      SphP[target].htest_t =htest_t;
#endif
    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].Acc[k] = acc[k];
      HydroDataResult[target].DtEntropy = dtEntropy;
      HydroDataResult[target].MaxSignalVel = maxSignalVel;
#ifdef NK_AV
      HydroDataResult[target].NumNK = numNK;
      HydroDataResult[target].NumN = numN;
#endif
#ifdef CDAV_DRIFTUPDATE
      for(k=0;k<9;k++)
      {
        HydroDataResult[target].T[k]=T[k];
        HydroDataResult[target].D[k]=D[k];
        HydroDataResult[target].E[k]=E[k];
      }
      HydroDataResult[target].R=R;
      HydroDataResult[target].CDAVSignalVel=CDAVvsig;
#endif
#ifdef VAR_H_TEST
      HydroDataResult[target].htest_t = htest_t;
#endif
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
