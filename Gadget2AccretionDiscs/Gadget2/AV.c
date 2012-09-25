/*
 * =====================================================================================
 *
 *       Filename:  AV.c
 *
 *    Description:  This contains the functions needed to perform an extra loop over neighbours as needed for the correct calculation of div.v and curl(v) required by the improved artificial viscosity prescription.
 *
 *        Version:  1.0
 *        Created:  21/09/12 14:05:19
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Matthew Young, 
 *        Company:  
 *
 * =====================================================================================
 */

/*! This function computes the matricies needed for a second order accurate estimate 
 * of div.v and div.a, both of which are used by the improved AV implementation.
 * Anything to do with periodic boundary conditions or cosmological integration 
 * is ignored.
 */

void artificialViscosity(void)
{
  /* The number of particles needing a new force.  If they don't need a new force, no need to calculate a new AV for them... */
  for(n=0, NumSphUpdate = 0; n < N_gas; n++)
  {
    if(P[n].Ti_endstep == All.Ti_Current)
    {
      NumSpdUpdate++;
    }
  }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist,1,MPI_INT,MPI_COMM_WORLD);
  for(i=0, ntot=0; i<NTask; i++)
  {
    ntot +=numlist[i];
  }
  free(numlist);

  noffset = malloc(sizeof(int) * NTask);
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  i=0; 
  ntotleft = ntot;

  while(ntotleft > 0)
  {
    for(j=0; j< NTask; j++)
    {
      nsend_local[j] = 0;
    }


    /*  do local particles and prepare the export list */
    tstart = second();
    for(nexport = 0, ndone = 0;i < N_gas && nexport < All.BunchSizeHydro - NTask; i++)
    {
      if(P[i].Ti_endstep == All.Ti_Current)
      {
        ndone++;

        for(j=0;j<NTask;j++)
        {
          Exportflag[j] = 0;
        }

        av_evaluate(i,0);

        for(j=0; j<NTask; j++)
        {
          if(Exportflag[j])
          {
            for(k = 0; k <3; k++)
            {
              AvDataIn[nexport].Pos[k] = P[i].Pos[k];
              AvDataIn[nexport].Vel[k] = SphP[i].VelPred[k];
              AvDataIn[nexport].Acc[k] = SphP[i].HydroAccel[k];
              AvDataIn[nexport].GradW[k] = SphP[i].GradW[k];
            }
            nexport++;
            nsend_local[j]++;
          }
        }
        tend = second();
        timecomp += timediff(tstart,tend);

        qsort(AvDataIn, nexport, sizeof(struct avdata_in), av_compare_key);

        for(j = 1, noffset[0] = 0; j<NTask; j++)
        {
          noffset[j] = noffset[j-1] + nsend_local[j-1];
        }

        tstart = second();

        MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

        tend = second();

        timeimbalance += timediff(tstart,tend);

        /*Now do the exported particles*/

        for(level = 1 ; level < (1 << PTask); level++)
        {
          tstart = second();
          for(j=0; j< NTask; j++)
          {
            nbuffer[j] = 0;
          }
          for(ngrp = level; ngrp < (1 << PTask); ngrp++)
          {
            maxfill = 0;
            for(j=0; j <NTask; j++)
            {
              if((j ^ ngrp) < NTask)
              {
                if(maxfill < nbuffer[j] + nsend[(j^ngrp) * NTask +j])
                {
                  maxfill = nbuffer[j] + nsend[(j^ngrp) * NTask +j];
                }
              }
            }
            if(maxfill >= All.BunchSizeHydro)
            {
              break;
            }

            sendTask = ThisTask;
            recvTask = ThisTask ^ ngrp;

            if(recvTask < NTask)
            {
              if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] >0)
              {
                /* Get the particle */
                MPI_Sendrecv(&AvDataIn[noffset[recvTask]], 
                nsend_local[recvTask] * sizeof(struct avdata_in), MPI_BYTE, 
                recvTask, TAG_AV_A, &AvDataGet[nbuffer[ThisTask]], 
                nsend[recvTask * NTask + ThisTask] * sizeof(struct avdata_in), 
                MPI_BYTE, recvTask, TAG_AV_A, MPI_COMM_WORLD, &status);
              }
            }

            for(j=0; j < NTask; j++)
            {
              if((j^ngrp) < NTask)
              {
                nbuffer[j] += nsend[(j ^ ngrp) *NTask + j];
              }
            }
          }
          tend = second();
          timecommsumm += timediff(tstart,tend);

          /* Finally calculate imported particles */
          tstart = second();
          for(j=0; j < nbuffer[ThisTask]; j++)
          {
            av_evaluate(j, 1);
          }
          tend=second();
          timecomp += timediff(tstart,tend);
          
          tstart = second();
          MPI_Barrier(MPI_COMM_WORLD);
          tend=second();
          timeimbalance += timediff(tstart, tend);

          /* Bring the results back locally */
          tstart = second();
          for(j=0; j< NTask; j++)
          {
            nbuffer[j] = 0;
          }
          for(ngrp = level; ngrp < (1 <<PTask); ngrp++)
          {
            maxfill = 0;
            for(j=0; j<NTask; j++)
            {
              if((j^ngrp) < NTask)
              {
                if(maxfill < nbuffer[j] + nsend[(j^ngrp)*NTask+j])
                {
                  maxfill = nbuffer[j]+nsend[(j^ngrp)*NTask+j];
                }
              }
            }
            if(maxfill >= All.BunchSizeHydro)
            {
              break;
            }

            sendTask = ThisTask;
            recvTask = ThisTask ^ ngrp;

            if(recvTask < NTask)
            {
              if(nsend[ThisTask * NTask +recvTask] >0 || nsend[recvTask *NTask + ThisTask] >0)
              {
                /*  send the results */
                MPI_Sendrecv(&AvDataResult[nbuffer[ThisTask]],
                    nsend[recvTask * NTask + ThisTask] * sizeof(struct avdata_out),
                    MPI_BYTE, recvTask, TAG_AV_B,
                    &AvDataPartialResult[noffset[recvTask]],
                    nsend_local[recvTask] * sizeof(struct avdata_out),
                    MPI_BYTE, recvTask, TAG_AV_B, MPI_COMM_WORLD, &status);

                /* add the reslt to particles */
                for(j= 0; j<nsend_local[recvTask]; j++)
                {
                  source = j+ noffset[recvTask];
                  place = AvDataIn[source].Index;
                  /* stuff goes here... */
                }
              }
            }

            for(j=0; j<NTask; j++)
            {
              if((j^ngrp) <NTask)
              {
                nbuffer[j] += nsend[(j^ngrp) * NTask +j];
              }
            }
          }

          tend = second();
          timecommsumm += timediff(tstart, tend);

          level = ngrp -1;
        }

        MPI_Allgather(&ndone, 1, MPI_INT, ndonelist,1,MPI_INT,MPI_COMM_WORLD);
        for(j=0;j<NTask;j++)
        {
          ntotleft -= ndonelist[j];
        }

        free(ndonelist);
        free(nsend);
        free(nbuffer);
        free(noffset);

        /*do final ops */
        tstart = second();

        /* stuff */

        tend = second();
        timecomp += timediff(tstart, tend);



      }
    }
  }
}




