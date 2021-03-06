#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*! \file gravtree.c 
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for
 *  all active local particles, and particles are exported to other
 *  processors if needed, where they can receive additional force
 *  contributions. If the TreePM algorithm is enabled, the force computed
 *  will only be the short-range part.
 */

/*! This function computes the gravitational forces for all active
 *  particles.  If needed, a new tree is constructed, otherwise the
 *  dynamically updated tree is used.  Particles are only exported to other
 *  processors when really needed, thereby allowing a good use of the
 *  communication buffer.
 */
void gravity_tree(void)
{
  long long ntot;
  int numnodes, nexportsum = 0;
  int i, j, iter = 0;
  int *numnodeslist, maxnumnodes, nexport, *numlist, *nrecv, *ndonelist;
  double tstart, tend, timetree = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  double ewaldcount;
  double costtotal, ewaldtot, *costtreelist, *ewaldlist;
  double maxt, sumt, *timetreelist, *timecommlist;
  double fac, plb, plb_max, sumcomm;

#ifndef NOGRAVITY
  int *noffset, *nbuffer, *nsend, *nsend_local;
  long long ntotleft;
  int ndone, maxfill, ngrp;
  int k, place;
  int level, sendTask, recvTask;
  double ax, ay, az;
  MPI_Status status;
#endif
#ifdef ADD_CENTRAL_GRAVITY
  int numsinks,root,globalroot,liveStar,liveStarGlobal;
  double starData[4],r,h,h_inv,h3_inv,u,starGrav[3],starGravGlobal[3];
#endif
  /* set new softening lengths */
  if(All.ComovingIntegrationOn)
    set_softenings();


  /* contruct tree if needed */
  tstart = second();
  if(TreeReconstructFlag)
    {
      if(ThisTask == 0)
	printf("Tree construction.\n");

      force_treebuild(NumPart);

      TreeReconstructFlag = 0;

      if(ThisTask == 0)
	printf("Tree construction done.\n");
    }
  tend = second();
  All.CPU_TreeConstruction += timediff(tstart, tend);

  costtotal = ewaldcount = 0;

  /* Note: 'NumForceUpdate' has already been determined in find_next_sync_point_and_drift() */
  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


#ifndef NOGRAVITY
  if(ThisTask == 0)
    printf("Begin tree force.\n");


#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(((1 << P[i].Type) & (SELECTIVE_NO_GRAVITY)))
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      iter++;

      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeForce - NTask; i++)
      {
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    ndone++;
	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;
#ifndef PMGRID
	    costtotal += force_treeevaluate(i, 0, &ewaldcount);
#else
	    costtotal += force_treeevaluate_shortrange(i, 0);
#endif
	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    for(k = 0; k < 3; k++)
		      GravDataGet[nexport].u.Pos[k] = P[i].Pos[k];
#ifdef UNEQUALSOFTENINGS
		    GravDataGet[nexport].Type = P[i].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		    if(P[i].Type == 0)
		      GravDataGet[nexport].Soft = SphP[i].Hsml;
#endif
#endif
		    GravDataGet[nexport].w.OldAcc = P[i].OldAcc;
		    GravDataIndexTable[nexport].Task = j;
		    GravDataIndexTable[nexport].Index = i;
		    GravDataIndexTable[nexport].SortIndex = nexport;
		    nexport++;
		    nexportsum++;
		    nsend_local[j]++;
		  }
	      }
	  }
      }
      tend = second();
      timetree += timediff(tstart, tend);

      qsort(GravDataIndexTable, nexport, sizeof(struct gravdata_index), grav_tree_compare_key);

      for(j = 0; j < nexport; j++)
	GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

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
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&GravDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A,
				   &GravDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
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
	    {
#ifndef PMGRID
	      costtotal += force_treeevaluate(j, 1, &ewaldcount);
#else
	      costtotal += force_treeevaluate_shortrange(j, 1);
#endif
	    }
	  tend = second();
	  timetree += timediff(tstart, tend);

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
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&GravDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B,
				   &GravDataOut[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  place = GravDataIndexTable[noffset[recvTask] + j].Index;

			  for(k = 0; k < 3; k++)
			    P[place].GravAccel[k] += GravDataOut[j + noffset[recvTask]].u.Acc[k];

			  P[place].GravCost += GravDataOut[j + noffset[recvTask]].w.Ninteractions;
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
  /* now add things for comoving integration */

#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      {
#ifdef PMGRID
	ax = P[i].GravAccel[0] + P[i].GravPM[0] / All.G;
	ay = P[i].GravAccel[1] + P[i].GravPM[1] / All.G;
	az = P[i].GravAccel[2] + P[i].GravPM[2] / All.G;
#else
	ax = P[i].GravAccel[0];
	ay = P[i].GravAccel[1];
	az = P[i].GravAccel[2];
#endif
	P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
      }


  if(All.TypeOfOpeningCriterion == 1)
    All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */

  /*  muliply by G */
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      for(j = 0; j < 3; j++)
      {
	     P[i].GravAccel[j] *= All.G;
        //printf("Gravity! %g\n",P[i].GravAccel[j]);
      }


  /* Finally, the following factor allows a computation of a cosmological simulation 
     with vacuum energy in physical coordinates */
#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn == 0)
    {
      fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif

  if(ThisTask == 0)
    printf("tree is done.\n");
#else /* gravity is switched off */

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;

#endif

#ifdef SINK_GRAV_ONLY
  sink_grav();
#endif


#ifdef ADD_CENTRAL_GRAVITY
  /* Get the position and mass of the central object and send it to everyone */
  numsinks=NumPart - N_gas;
  starData[0]=starData[1]=starData[2]=starData[3]= -1.0;
  root=-1;
  liveStar=0;
  for(i=0; i<3; i++)
  {
    starGrav[i]=0.0;
    starGravGlobal[i]=0.0;
  }
  for(i=0; i<numsinks;i++)
  {
    if(P[i+N_gas].ID==All.StarID)
    {
      starData[0] = P[i+N_gas].Pos[0];
      starData[1] = P[i+N_gas].Pos[1];
      starData[2] = P[i+N_gas].Pos[2];
      starData[3] = P[i+N_gas].Mass;
      root = ThisTask;
      //Do we need to update the star's gravity?
      if(P[i+N_gas].Ti_endstep == All.Ti_Current)
      {
        liveStar=1;
      }
    }
  }
  /* Get the node that has the data */
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&root,&globalroot,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&liveStar,&liveStarGlobal,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  /* Broadcast it. */
  MPI_Bcast(&starData,4,MPI_DOUBLE,globalroot,MPI_COMM_WORLD);
  //We have the central object mass and position, add its gravity, it is softened by the type 1 softening...
  h = All.ForceSoftening[1];
  h_inv = 1.0 / h;
  h3_inv = h_inv * h_inv * h_inv;
  for(i = 0; i < NumPart; i++)
  {
    if(P[i].ID != All.StarID)
    {
      //If we need to update the star's gravity we need to calculate this for all particles...
      if(liveStarGlobal)
      {
        r=sqrt((starData[0]-P[i].Pos[0])*(starData[0]-P[i].Pos[0])+(starData[1]-P[i].Pos[1])*(starData[1]-P[i].Pos[1])+(starData[2]-P[i].Pos[2])*(starData[2]-P[i].Pos[2]));
        if(r >= h)
        {
	         fac = 1 / (r*r*r);
        }
        else
        {
  	       u = r * h_inv;
  	       if(u < 0.5)
          {
  	         fac = h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
          }
  	       else
          {
  	         fac = h3_inv * (21.333333333333 - 48.0 * u +
  		  	   38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
          }
        }
        for(j=0;j<3;j++)
        {
          starGrav[j]+=(P[i].Pos[j]-starData[j])*All.G*P[i].Mass*fac;
        }
      }
      //Otherwise, just give the star's gravity to those that need it
      if(P[i].Ti_endstep == All.Ti_Current)
      {
        r=sqrt((starData[0]-P[i].Pos[0])*(starData[0]-P[i].Pos[0])+(starData[1]-P[i].Pos[1])*(starData[1]-P[i].Pos[1])+(starData[2]-P[i].Pos[2])*(starData[2]-P[i].Pos[2]));
        if(r >= h)
        {
	         fac = 1 / (r*r*r);
        }
        else
        {
  	       h_inv = 1.0 / h;
  	       h3_inv = h_inv * h_inv * h_inv;
  	       u = r * h_inv;
  	       if(u < 0.5)
          {
  	         fac = h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
          }
  	       else
          {
  	         fac = h3_inv * (21.333333333333 - 48.0 * u +
  		  	   38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
          }
        }
        for(j=0;j<3;j++)
        {
          P[i].GravAccel[j]+=(starData[j]-P[i].Pos[j])*All.G*starData[3]*fac;
        }
      }
    }
  }
  //Gather the forces of the pcles on the star together and add them to the star
  if(liveStarGlobal)
  {
    //Finally we need to combine all the starGrav values for the star...
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&starGrav[0],&starGravGlobal[0],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&starGrav[1],&starGravGlobal[1],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&starGrav[2],&starGravGlobal[2],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //Finally, find the actual star and add it
    if(globalroot==ThisTask)
    {
      for(i=0; i<numsinks;i++)
      {
        if(P[i+N_gas].ID==All.StarID)
        {
          for(j=0;j<3;j++)
          {
            P[i+N_gas].GravAccel[j]+=starGravGlobal[j];
          }
        }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif



  /* Now the force computation is finished */

  /*  gather some diagnostic information */

  timetreelist = malloc(sizeof(double) * NTask);
  timecommlist = malloc(sizeof(double) * NTask);
  costtreelist = malloc(sizeof(double) * NTask);
  numnodeslist = malloc(sizeof(int) * NTask);
  ewaldlist = malloc(sizeof(double) * NTask);
  nrecv = malloc(sizeof(int) * NTask);

  numnodes = Numnodestree;

  MPI_Gather(&costtotal, 1, MPI_DOUBLE, costtreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&numnodes, 1, MPI_INT, numnodeslist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&timecommsumm, 1, MPI_DOUBLE, timecommlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&NumPart, 1, MPI_INT, nrecv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&ewaldcount, 1, MPI_DOUBLE, ewaldlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Reduce(&nexportsum, &nexport, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.TotNumOfForces += ntot;

      fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
      fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g  iter= %d\n",
	      (int) (ntot / 1000000000), (int) (ntot % 1000000000),
	      (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
	      nexport / ((double) ntot), iter);
      /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */

      fac = NTask / ((double) All.TotNumPart);

      for(i = 0, maxt = timetreelist[0], sumt = 0, plb_max = 0,
	  maxnumnodes = 0, costtotal = 0, sumcomm = 0, ewaldtot = 0; i < NTask; i++)
	{
	  costtotal += costtreelist[i];

	  sumcomm += timecommlist[i];

	  if(maxt < timetreelist[i])
	    maxt = timetreelist[i];
	  sumt += timetreelist[i];

	  plb = nrecv[i] * fac;

	  if(plb > plb_max)
	    plb_max = plb;

	  if(numnodeslist[i] > maxnumnodes)
	    maxnumnodes = numnodeslist[i];

	  ewaldtot += ewaldlist[i];
	}
      fprintf(FdTimings, "work-load balance: %g  max=%g avg=%g PE0=%g\n",
	      maxt / (sumt / NTask), maxt, sumt / NTask, timetreelist[0]);
      fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
      fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
	      maxnumnodes / (All.TreeAllocFactor * All.MaxPart));
      fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", ntot / (sumt + 1.0e-20),
	      ntot / (maxt * NTask), ((double) (costtotal)) / ntot, ((double) ewaldtot) / ntot);
      fprintf(FdTimings, "\n");

      fflush(FdTimings);

      All.CPU_TreeWalk += sumt / NTask;
      All.CPU_Imbalance += sumimbalance / NTask;
      All.CPU_CommSum += sumcomm / NTask;
    }

  free(nrecv);
  free(ewaldlist);
  free(numnodeslist);
  free(costtreelist);
  free(timecommlist);
  free(timetreelist);
}



/*! This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].  We check that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
  int i;

  if(All.ComovingIntegrationOn)
    {
      if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
        All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
      else
        All.SofteningTable[0] = All.SofteningGas;
      
      if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
        All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
      else
        All.SofteningTable[1] = All.SofteningHalo;
      
      if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
        All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
      else
        All.SofteningTable[2] = All.SofteningDisk;
      
      if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
        All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
      else
        All.SofteningTable[3] = All.SofteningBulge;
      
      if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
        All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
      else
        All.SofteningTable[4] = All.SofteningStars;
      
      if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
        All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
      else
        All.SofteningTable[5] = All.SofteningBndry;
    }
  else
    {
      All.SofteningTable[0] = All.SofteningGas;
      All.SofteningTable[1] = All.SofteningHalo;
      All.SofteningTable[2] = All.SofteningDisk;
      All.SofteningTable[3] = All.SofteningBulge;
      All.SofteningTable[4] = All.SofteningStars;
      All.SofteningTable[5] = All.SofteningBndry;
    }

  for(i = 0; i < 6; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  All.MinGasHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];
}


/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU.
 */
int grav_tree_compare_key(const void *a, const void *b)
{
  if(((struct gravdata_index *) a)->Task < (((struct gravdata_index *) b)->Task))
    return -1;

  if(((struct gravdata_index *) a)->Task > (((struct gravdata_index *) b)->Task))
    return +1;

  return 0;
}

#ifdef SINK_GRAV_ONLY
void sink_grav(void)
{
  int *noffset;
  int *nbuffer;
  int nsinks;
  struct particle_data * sinks;
  int i,j;
  double dx,dy,dz,r,fac,u,h;
  noffset = malloc(NTask * sizeof(int));
  nbuffer = malloc(NTask * sizeof(int));

  nsinks = NumPart-N_gas;
  MPI_Allgather(&nsinks,1,MPI_INT,nbuffer,1,MPI_INT,MPI_COMM_WORLD);
  nsinks=0;
  for(i=0;i<NTask;i++)
    nsinks+=nbuffer[i];
  //Get all the type 1 particles on all processors
  sinks = malloc(nsinks*sizeof(struct particle_data));
  for(j=1,noffset[0]=0,nbuffer[0]=nbuffer[0]*sizeof(struct particle_data);j<NTask;j++)
  {
    nbuffer[j] = nbuffer[j] * sizeof(struct particle_data);
    noffset[j] = noffset[j-1] + nbuffer[j-1];
  }

  MPI_Allgatherv(&P[N_gas],sizeof(struct particle_data) * (NumPart-N_gas),MPI_BYTE,
      sinks,nbuffer,noffset,MPI_BYTE,MPI_COMM_WORLD);

  //The force softening for type 1 particles...
  //The usual GADGET factor of 2.8 is included (even though I don't
  //know why it's used)
  h=2.8*All.SofteningHalo;
  for(i=0;i<NumPart;i++)
  {
    if(P[i].Ti_endstep == All.Ti_Current)
    {
      //This particle needs doing
      P[i].GravAccel[0] = P[i].GravAccel[1] = P[i].GravAccel[2] = 0;
      //Add the gravitational force due to each type 1 particle
      for(j=0;j<nsinks;j++)
      {
        if(P[i].ID==sinks[j].ID)
        {
          continue;
        }
        dx = P[i].Pos[0]-sinks[j].Pos[0];
        dy = P[i].Pos[1]-sinks[j].Pos[1];
#ifdef TWODIMS
        dz = 0;
#else
        dz = P[i].Pos[2]-sinks[j].Pos[2];
#endif
        r = sqrt(dx*dx+dy*dy+dz*dz);
        if(r >= h)
        {
	         fac = 1.0 / (r*r*r);
        }
        else
        {
          //Plumber sphere!
          fac = 1.0 / (h*h*h*pow(1.0+(r/h)*(r/h),1.5));
          //Usual crap
  	       //u = r / h;
  	       //if(u < 0.5)
          //{
  	       //  fac = (1.0/(h*h*h)) * (10.666666666667 + u * u * (32.0 * u - 38.4));
          //}
  	       //else
          //{
  	       //  fac = (1.0/(h*h*h)) * (21.333333333333 - 48.0 * u +
  		  	 //  38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
          //}
        }
        if(P[i].Type==1 && P[i].Mass>1)
        {
          if(sinks[j].Mass>1)
          {
            printf("Self gravity factor %g! close approach = %d \n",fac,r<h);
          }
          else
          {
            printf("factor due to companion %g. close approach = %d \n",fac,r<h);
          }
        }

        fac *= All.G*sinks[j].Mass;
        P[i].GravAccel[0] -= fac * dx;
        P[i].GravAccel[1] -= fac * dy;
        P[i].GravAccel[2] -= fac * dz;
      }
      if(P[i].Type==1)
        printf("[%d] P[%d] gravity = (%g,%g,%g).\n",ThisTask,i,P[i].GravAccel[0],P[i].GravAccel[1],P[i].GravAccel[2]);
    }
  }
  //Resurrect gas
  for(i = 0; i < NumPart; i++)
  {
    if(P[i].Ti_endstep < 0)
    {
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  //Freedom!
  free(sinks);
  free(nbuffer);
  free(noffset);
}
#endif
