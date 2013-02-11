#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file density.c 
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and
 *  some auxiliary quantities are computed.  If the number of neighbours
 *  obtained falls outside the target range, the correct smoothing
 *  length is determined iteratively, if needed.
 */


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


/*! This function computes the local density for each active SPH particle,
 *  the number of neighbours in the current smoothing radius, and the
 *  divergence and curl of the velocity field.  The pressure is updated as
 *  well.  If a particle with its smoothing region is fully inside the
 *  local domain, it is not exported to the other processors. The function
 *  also detects particles that have a number of neighbours outside the
 *  allowed tolerance range. For these particles, the smoothing length is
 *  adjusted accordingly, and the density computation is executed again.
 *  Note that the smoothing length is not allowed to fall below the lower
 *  bound set by MinGasHsml.
 */
void density(void)
{
  long long ntot, ntotleft;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, n, ndone, npleft, maxfill, source, iter = 0;
  int level, ngrp, sendTask, recvTask, place, nexport;
  double dt_entr, tstart, tend, tstart_ngb = 0, tend_ngb = 0;
  double sumt, sumcomm, timengb, sumtimengb;
  double timecomp = 0, timeimbalance = 0, timecommsumm = 0, sumimbalance;
  MPI_Status status;
#ifdef CDAV
  int k;
  double Tinv[6],V[9],fac,dt_alpha,alphaloc;
  double divv,diva,xi,A;
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


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      SphP[n].Left = SphP[n].Right = 0;
#if defined CDAV
      SphP[n].AlphaOld=-1;
#endif

      if(P[n].Ti_endstep == All.Ti_Current)
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);



  /* we will repeat the whole thing for those particles where we didn't
   * find enough neighbours
   */
  do
    {
      i = 0;			/* begin with this index */
      ntotleft = ntot;		/* particles left for all tasks together */

      while(ntotleft > 0)
	{
	  for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;

	  /* do local particles and prepare export list */
	  tstart = second();
	  for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++)
	    if(P[i].Ti_endstep == All.Ti_Current)
	      {
		ndone++;


		for(j = 0; j < NTask; j++)
		  Exportflag[j] = 0;

		density_evaluate(i, 0);

		for(j = 0; j < NTask; j++)
		  {
		    if(Exportflag[j])
		      {
			DensDataIn[nexport].Pos[0] = P[i].Pos[0];
			DensDataIn[nexport].Pos[1] = P[i].Pos[1];
			DensDataIn[nexport].Pos[2] = P[i].Pos[2];
			DensDataIn[nexport].Vel[0] = SphP[i].VelPred[0];
			DensDataIn[nexport].Vel[1] = SphP[i].VelPred[1];
			DensDataIn[nexport].Vel[2] = SphP[i].VelPred[2];
#ifdef CDAV
         DensDataIn[nexport].Accel[0] = SphP[i].HydroAccel[0]+P[i].GravAccelOld[0];
         DensDataIn[nexport].Accel[1] = SphP[i].HydroAccel[1]+P[i].GravAccelOld[1];
         DensDataIn[nexport].Accel[2] = SphP[i].HydroAccel[2]+P[i].GravAccelOld[2];
#endif
			DensDataIn[nexport].Hsml = SphP[i].Hsml;
			DensDataIn[nexport].Index = i;
			DensDataIn[nexport].Task = j;
			nexport++;
			nsend_local[j]++;
		      }
		  }
	      }
	  tend = second();
	  timecomp += timediff(tstart, tend);

	  qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

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
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* get the particles */
			  MPI_Sendrecv(&DensDataIn[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				       recvTask, TAG_DENS_A,
				       &DensDataGet[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
				       MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);


	      tstart = second();
	      for(j = 0; j < nbuffer[ThisTask]; j++)
		density_evaluate(j, 1);
	      tend = second();
	      timecomp += timediff(tstart, tend);

	      /* do a block to explicitly measure imbalance */
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
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* send the results */
			  MPI_Sendrecv(&DensDataResult[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B,
				       &DensDataPartialResult[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
			  for(j = 0; j < nsend_local[recvTask]; j++)
			    {
			      source = j + noffset[recvTask];
			      place = DensDataIn[source].Index;

			      SphP[place].NumNgb += DensDataPartialResult[source].Ngb;
			      SphP[place].Density += DensDataPartialResult[source].Rho;
			      SphP[place].DivVel += DensDataPartialResult[source].Div;

			      SphP[place].DhsmlDensityFactor += DensDataPartialResult[source].DhsmlDensity;
#ifdef PRICE_GRAV_SOFT
               SphP[place].Zeta += DensDataPartialResult[source].Zeta;
#endif
#ifdef CDAV
               for(k=0;k<9;k++)
               {
                 if(k<6)
                 {
                   SphP[place].T[k] += DensDataPartialResult[source].T[k];
                 }
                 SphP[place].D[k] += DensDataPartialResult[source].D[k];
                 SphP[place].E[k] += DensDataPartialResult[source].E[k];
               }
               SphP[place].R += DensDataPartialResult[source].R;
#endif
			      SphP[place].Rot[0] += DensDataPartialResult[source].Rot[0];
			      SphP[place].Rot[1] += DensDataPartialResult[source].Rot[1];
			      SphP[place].Rot[2] += DensDataPartialResult[source].Rot[2];
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



      /* do final operations on results */
      tstart = second();
      for(i = 0, npleft = 0; i < N_gas; i++)
	{
	  if(P[i].Ti_endstep == All.Ti_Current)
	    {
         {
		SphP[i].DhsmlDensityFactor =
		  1 / (1 + SphP[i].Hsml * SphP[i].DhsmlDensityFactor / (NUMDIMS * SphP[i].Density));
#ifdef PRICE_GRAV_SOFT
      SphP[i].Zeta = -SphP[i].Hsml * SphP[i].Zeta / ( NUMDIMS * SphP[i].Density);
#endif

		SphP[i].CurlVel = sqrt(SphP[i].Rot[0] * SphP[i].Rot[0] +
				       SphP[i].Rot[1] * SphP[i].Rot[1] +
				       SphP[i].Rot[2] * SphP[i].Rot[2]) / SphP[i].Density;

		SphP[i].DivVel /= SphP[i].Density;
		dt_entr = (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

		SphP[i].Pressure =
		  (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);

#ifdef CDAV
      SphP[i].R /= SphP[i].Density;
      //Inverse of determinate of T
      fac = 1/(SphP[i].T[0]*(SphP[i].T[3]*SphP[i].T[5]-SphP[i].T[4]*SphP[i].T[4]) +
        SphP[i].T[1]*(SphP[i].T[2]*SphP[i].T[4]-SphP[i].T[1]*SphP[i].T[5]) +
        SphP[i].T[2]*(SphP[i].T[1]*SphP[i].T[4]-SphP[i].T[2]*SphP[i].T[3]));
      //The inverse of the matrix T
      Tinv[0]=fac*(SphP[i].T[3]*SphP[i].T[5]-SphP[i].T[4]*SphP[i].T[4]);
      Tinv[1]=fac*(SphP[i].T[2]*SphP[i].T[4]-SphP[i].T[1]*SphP[i].T[5]);
      Tinv[2]=fac*(SphP[i].T[1]*SphP[i].T[4]-SphP[i].T[2]*SphP[i].T[3]);
      Tinv[3]=fac*(SphP[i].T[0]*SphP[i].T[5]-SphP[i].T[2]*SphP[i].T[2]);
      Tinv[4]=fac*(SphP[i].T[1]*SphP[i].T[2]-SphP[i].T[0]*SphP[i].T[4]);
      Tinv[5]=fac*(SphP[i].T[0]*SphP[i].T[3]-SphP[i].T[1]*SphP[i].T[1]);
      //The velocity matrix D.  There's actually a mistake in the paper.  D.Tinv = V^t not V
      V[0]=SphP[i].D[0]*Tinv[0]+SphP[i].D[1]*Tinv[1]+SphP[i].D[2]*Tinv[2];
      V[3]=SphP[i].D[0]*Tinv[1]+SphP[i].D[1]*Tinv[3]+SphP[i].D[2]*Tinv[4];
      V[6]=SphP[i].D[0]*Tinv[2]+SphP[i].D[1]*Tinv[4]+SphP[i].D[2]*Tinv[5];
      V[1]=SphP[i].D[3]*Tinv[0]+SphP[i].D[4]*Tinv[1]+SphP[i].D[5]*Tinv[2];
      V[4]=SphP[i].D[3]*Tinv[1]+SphP[i].D[4]*Tinv[3]+SphP[i].D[5]*Tinv[4];
      V[7]=SphP[i].D[3]*Tinv[2]+SphP[i].D[4]*Tinv[4]+SphP[i].D[5]*Tinv[5];
      V[2]=SphP[i].D[6]*Tinv[0]+SphP[i].D[7]*Tinv[1]+SphP[i].D[8]*Tinv[2];
      V[5]=SphP[i].D[6]*Tinv[1]+SphP[i].D[7]*Tinv[3]+SphP[i].D[8]*Tinv[4];
      V[8]=SphP[i].D[6]*Tinv[2]+SphP[i].D[7]*Tinv[4]+SphP[i].D[8]*Tinv[5];
      //DivVel is now trivially estimated
      divv = V[0]+V[4]+V[8];
      //DivAccel = tr(E.T^-1)-tr(V^2)
      //This is the first part
      diva = (SphP[i].E[1]+SphP[i].E[3])*Tinv[1]+(SphP[i].E[2]+SphP[i].E[6])*Tinv[2]+
        (SphP[i].E[5]+SphP[i].E[7])*Tinv[4] + SphP[i].E[0]*Tinv[0] + 
        SphP[i].E[4]*Tinv[3]+SphP[i].E[8]*Tinv[5];
      //now subtract the second part...
      diva =diva-(V[0]*V[0]+V[4]*V[4]+V[8]*V[8]+2*(V[1]*V[3]+V[2]*V[6]+V[5]*V[7]));
      //div.v is only calculated at the force updates, hence the appropriate dt is the size of the current timestep
      dt_alpha=(P[i].Ti_endstep-P[i].Ti_begstep)*All.Timebase_interval;
      //Alternative method for calculating diva
      //if(dt_alpha!=0)
      //{
      //  A=(divv-SphP[i].DivVelOld)/dt_alpha;
      //}
      //else
      //{
      //  A=0;
      //}
      //Test if the newly calculated quantities return more or less what we expect from the old estimates.  In the case of diva the "old" estimate is from the difference between timesteps.
      //if(ThisTask==0)
      //{
      //  printf("(old,new,new/old) divv = (%g,%g,%g), diva = (%g,%g,%g)\n",SphP[i].DivVel,divv,divv/SphP[i].DivVel,A,diva,diva/A);
      //  //Outputs old,new,new/old for divv,diva,curl_x,curl_y,curl_z
      //  printf("sheemba%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",SphP[i].DivVel,divv,divv/SphP[i].DivVel,A,diva,diva/A,SphP[i].Rot[0]/SphP[i].Density,V[5]-V[7],(SphP[i].Density*(V[5]-V[7]))/SphP[i].Rot[0],SphP[i].Rot[1]/SphP[i].Density,V[6]-V[2],(SphP[i].Density*(V[6]-V[2]))/SphP[i].Rot[1],SphP[i].Rot[2]/SphP[i].Density,V[1]-V[3],(SphP[i].Density*(V[1]-V[3]))/SphP[i].Rot[2]);
      //}
      //diva=A;
      //SphP[i].DivVel = divv;

#ifdef NOBALSARA
      xi=1;
#else
      //Now calculate xi
      xi=2*(1-SphP[i].R)*(1-SphP[i].R)*(1-SphP[i].R)*(1-SphP[i].R)*divv;
      xi *= xi;
      //The added bit is tr(S.S^T)
      //This is just a temp variable, not really A
      A=(V[0]*V[0]+V[4]*V[4]+V[8]*V[8]) + 
        0.5*((V[1]+V[3])*(V[1]+V[3])+(V[2]+V[6])*(V[2]+V[6])+(V[5]+V[7])*(V[5]+V[7])) -
        (divv*divv/3);
      //If this is true, xi=0 already...
      if(xi+A!=0)
      {
        xi = xi/(xi+A);
      }
#endif
      //Now calculate A_i
      A=dmax(0,-1*diva)*xi;
      //Now want to adapt alpha.  Need to be careful here as there is a possibility
      //that this entire loop will have to repeat itself, in which case we want this
      //adaption step to see the same current alpha every time.  This is the purpose
      //of "alphaold" which is initialised to 0 at the start of each particle.
      alphaloc = SphP[i].MaxSignalVel*SphP[i].MaxSignalVel + SphP[i].Hsml*SphP[i].Hsml*A;
      //If everything is zero, leave alpha at 0
      if(alphaloc!=0)
      {
        alphaloc = (All.ArtBulkViscConst*SphP[i].Hsml*SphP[i].Hsml*A)/alphaloc;
      }
      //This is so alpha isn't "set" twice if the density routine has to be rerun...
      if(SphP[i].AlphaOld==-1)
      {
        SphP[i].AlphaOld=SphP[i].Alpha;
      }
      else
      {
        SphP[i].Alpha=SphP[i].AlphaOld;
      }
      //Finally, advance the artificial viscosity value to the new value, if it's the first time step, don't.
      if(All.Ti_Current)
      {
        if(SphP[i].Alpha < alphaloc)
        {
          SphP[i].Alpha = alphaloc;
          SphP[i].DtAlpha = 0;
        }
        else
        {
          //SphP[i].Alpha = alphaloc +(SphP[i].Alpha-alphaloc)*exp(-2*All.VariableViscDecayLength*SphP[i].MaxSignalVel*dt_alpha/SphP[i].Hsml);
          SphP[i].DtAlpha = (alphaloc-SphP[i].Alpha)*2*All.VariableViscDecayLength * SphP[i].MaxSignalVel / SphP[i].Hsml;
        }
      }
#endif
         }


	      /* now check whether we had enough neighbours */

	      if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation) ||
		 (SphP[i].NumNgb > (All.DesNumNgb + All.MaxNumNgbDeviation)
		  && SphP[i].Hsml > (1.01 * All.MinGasHsml)))
		{
		  /* need to redo this particle */
		  npleft++;

		  if(SphP[i].Left > 0 && SphP[i].Right > 0)
		    if((SphP[i].Right - SphP[i].Left) < 1.0e-3 * SphP[i].Left)
		      {
			/* this one should be ok */
			npleft--;
			P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
			continue;
		      }

		  if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation))
		    SphP[i].Left = dmax(SphP[i].Hsml, SphP[i].Left);
		  else
		    {
		      if(SphP[i].Right != 0)
			{
			  if(SphP[i].Hsml < SphP[i].Right)
			    SphP[i].Right = SphP[i].Hsml;
			}
		      else
			SphP[i].Right = SphP[i].Hsml;
		    }

		  if(iter >= MAXITER - 10)
		    {
		      printf
			("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (int) P[i].ID, SphP[i].Hsml, SphP[i].Left, SphP[i].Right,
			 (float) SphP[i].NumNgb, SphP[i].Right - SphP[i].Left, P[i].Pos[0], P[i].Pos[1],
			 P[i].Pos[2]);
		      fflush(stdout);
		    }

		  if(SphP[i].Right > 0 && SphP[i].Left > 0)
		    SphP[i].Hsml = pow(0.5 * (pow(SphP[i].Left, 3) + pow(SphP[i].Right, 3)), 1.0 / 3);
		  else
		    {
		      if(SphP[i].Right == 0 && SphP[i].Left == 0)
			endrun(8188);	/* can't occur */

		      if(SphP[i].Right == 0 && SphP[i].Left > 0)
			{
			  if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
			    {
			      SphP[i].Hsml *=
				1 - (SphP[i].NumNgb -
				     All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;
			    }
			  else
			    SphP[i].Hsml *= 1.26;
			}

		      if(SphP[i].Right > 0 && SphP[i].Left == 0)
			{
			  if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
			    {
			      SphP[i].Hsml *=
				1 - (SphP[i].NumNgb -
				     All.DesNumNgb) / (NUMDIMS * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;
			    }
			  else
			    SphP[i].Hsml /= 1.26;
			}
		    }

		  if(SphP[i].Hsml < All.MinGasHsml)
		    SphP[i].Hsml = All.MinGasHsml;
		}
	      else
		P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
	    }
	}
      tend = second();
      timecomp += timediff(tstart, tend);


      numlist = malloc(NTask * sizeof(int) * NTask);
      MPI_Allgather(&npleft, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
      for(i = 0, ntot = 0; i < NTask; i++)
	ntot += numlist[i];
      free(numlist);

      if(ntot > 0)
	{
	  if(iter == 0)
	    tstart_ngb = second();

	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
      else
	tend_ngb = second();
    }
  while(ntot > 0);


  /* mark as active again */
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);


  /* collect some timing information */
  if(iter > 0)
    timengb = timediff(tstart_ngb, tend_ngb);
  else
    timengb = 0;

  MPI_Reduce(&timengb, &sumtimengb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += sumimbalance / NTask;
      All.CPU_EnsureNgb += sumtimengb / NTask;
    }
}



/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
void density_evaluate(int target, int mode)
{
  int j, n, startnode, numngb, numngb_inbox;
  double h, h2, fac, hinv, hinv3, hinv4;
  double rho, divv, wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  double dvx, dvy, dvz, rotv[3];
  double weighted_numngb, dhsmlrho;
#ifdef PRICE_GRAV_SOFT
  double zeta,dphi,hinv2;
#endif
#ifdef CDAV
  int i;
  double D[9],E[9],T[6];
  double dax,day,daz;
  double R,divvsign;
  double acc[3];
#endif
  FLOAT *pos, *vel;

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h = SphP[target].Hsml;
#ifdef CDAV
      for(i=0;i<3;i++)
      {
        acc[i] = SphP[target].HydroAccel+P[target].GravAccelOld[i];
      }
#endif
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;
#ifdef CDAV
      for(i=0;i<3;i++)
      {
        acc[i] = DensDataGet[target].Accel[i];
      }
#endif
    }

  h2 = h * h;
  hinv = 1.0 / h;
#ifdef PRICE_GRAV_SOFT
  hinv2 = hinv * hinv;
#endif
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;

  rho = divv = rotv[0] = rotv[1] = rotv[2] = 0;
  weighted_numngb = 0;
  dhsmlrho = 0;
#ifdef PRICE_GRAV_SOFT
  zeta = 0;
#endif
#ifdef CDAV
  for(i=0;i<9;i++)
  {
    if(i<6)
    {
      T[i]=0;
    }
    E[i]=D[i]=0;
  }
  R=0;
#endif

  startnode = All.MaxPart;
  numngb = 0;
  do
    {
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
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

	  if(r2 < h2)
	    {
	      numngb++;

	      r = sqrt(r2);

	      u = r * hinv;

	      if(u < 0.5)
		{
		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		  dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
#ifdef PRICE_GRAV_SOFT
        dphi = hinv2 * (-16.0 * u*u +48.0*u*u*u*u - 38.4*u*u*u*u*u +2.8);
#endif
		}
	      else
		{
		  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		  dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
#ifdef PRICE_GRAV_SOFT
        dphi = hinv2 * (-32.0 * u*u +64.0*u*u*u-48.0*u*u*u*u+12.8*u*u*u*u*u+3.2);
#endif
		}

	      mass_j = P[j].Mass;

	      rho += mass_j * wk;

	      weighted_numngb += NORM_COEFF * wk / hinv3;

	      dhsmlrho += -mass_j * (NUMDIMS * hinv * wk + u * dwk);
#ifdef CDAV
         if(P[j].Type==0)
         {
           divvsign=0;
           if(SphP[j].DivVel>0)
           {
             divvsign=1;
           }
           if(SphP[j].DivVel<0)
           {
             divvsign=-1;
           }
           R += divvsign * mass_j * wk;
         }
#endif

#ifdef PRICE_GRAV_SOFT
         zeta += mass_j * dphi;
#endif


	      if(r > 0)
		{
		  fac = mass_j * dwk / r;

		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];
		  divv -= fac * (dx * dvx + dy * dvy + dz * dvz);

		  rotv[0] += fac * (dz * dvy - dy * dvz);
		  rotv[1] += fac * (dx * dvz - dz * dvx);
		  rotv[2] += fac * (dy * dvx - dx * dvy);

#ifdef CDAV
        dax = acc[0] - SphP[j].HydroAccel[0]-P[j].GravAccelOld[0];
        day = acc[1] - SphP[j].HydroAccel[1]-P[j].GravAccelOld[1];
        daz = acc[2] - SphP[j].HydroAccel[2]-P[j].GravAccelOld[2];
        //The factors of h in fac are irrelevant as they appear in both D and T and so will cancel each other out
        D[0] += fac * (dvx*dx);
        D[1] += fac * (dvx*dy);
        D[2] += fac * (dvx*dz);
        D[3] += fac * (dvy*dx);
        D[4] += fac * (dvy*dy);
        D[5] += fac * (dvy*dz);
        D[6] += fac * (dvz*dx);
        D[7] += fac * (dvz*dy);
        D[8] += fac * (dvz*dz);

        E[0] += fac * (dax*dx);
        E[1] += fac * (dax*dy);
        E[2] += fac * (dax*dz);
        E[3] += fac * (day*dx);
        E[4] += fac * (day*dy);
        E[5] += fac * (day*dz);
        E[6] += fac * (daz*dx);
        E[7] += fac * (daz*dy);
        E[8] += fac * (daz*dz);

        T[0] += fac * (dx*dx);
        T[1] += fac * (dx*dy);
        T[2] += fac * (dx*dz);
        T[3] += fac * (dy*dy);
        T[4] += fac * (dy*dz);
        T[5] += fac * (dz*dz);
#endif

		}
	    }
	}
    }
  while(startnode >= 0);

  if(mode == 0)
    {
      SphP[target].NumNgb = weighted_numngb;
      SphP[target].Density = rho;
      SphP[target].DivVel = divv;
      SphP[target].DhsmlDensityFactor = dhsmlrho;
      SphP[target].Rot[0] = rotv[0];
      SphP[target].Rot[1] = rotv[1];
      SphP[target].Rot[2] = rotv[2];
#ifdef PRICE_GRAV_SOFT
      SphP[target].Zeta = zeta;
#endif
#ifdef CDAV
      for(i=0;i<9;i++)
      {
        if(i<6)
        {
          SphP[target].T[i]=T[i];

        }
        SphP[target].D[i]=D[i];
        SphP[target].E[i]=E[i];
      }
      SphP[target].R=R;
#endif
    }
  else
    {
      DensDataResult[target].Rho = rho;
      DensDataResult[target].Ngb = weighted_numngb;
      DensDataResult[target].DhsmlDensity = dhsmlrho;
      DensDataResult[target].Div = divv;
      DensDataResult[target].Rot[0] = rotv[0];
      DensDataResult[target].Rot[1] = rotv[1];
      DensDataResult[target].Rot[2] = rotv[2];
#ifdef PRICE_GRAV_SOFT
      DensDataResult[target].Zeta = zeta;
#endif
#ifdef CDAV
      for(i=0; i<9; i++)
      {
        if(i<6)
        {
          DensDataResult[target].T[i]=T[i];
        }
        DensDataResult[target].D[i]=D[i];
        DensDataResult[target].E[i]=E[i];
      }
      DensDataResult[target].R = R;
#endif
    }
}




/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int dens_compare_key(const void *a, const void *b)
{
  if(((struct densdata_in *) a)->Task < (((struct densdata_in *) b)->Task))
    return -1;

  if(((struct densdata_in *) a)->Task > (((struct densdata_in *) b)->Task))
    return +1;

  return 0;
}
