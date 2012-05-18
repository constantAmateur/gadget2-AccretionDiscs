#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


#ifdef SINK_PARTICLES
static double *list_sink_posx;
static double *list_sink_posy;
static double *list_sink_posz;
static double *list_sink_velx;
static double *list_sink_vely;
static double *list_sink_velz;
static int *list_sink_ID;
static double *list_sink_mass;
#endif


/*! \file predict.c
 *  \brief drift particles by a small time interval
 *
 *  This function contains code to implement a drift operation on all the
 *  particles, which represents one part of the leapfrog integration scheme.
 */


/*! This function drifts all particles from the current time to the future:
 *  time0 - > time1
 *
 *  If there is no explicit tree construction in the following timestep, the
 *  tree nodes are also drifted and updated accordingly. Note: For periodic
 *  boundary conditions, the mapping of coordinates onto the interval
 *  [0,All.BoxSize] is only done before the domain decomposition, or for
 *  outputs to snapshot files.  This simplifies dynamic tree updates, and
 *  allows the domain decomposition to be carried out only every once in a
 *  while.
 */
void move_particles(int time0, int time1)
{
  int i, j;
  double dt_drift, dt_gravkick, dt_hydrokick, dt_entr;
  double t0, t1;


  t0 = second();

  if(All.ComovingIntegrationOn)
    {
      dt_drift = get_drift_factor(time0, time1);
      dt_gravkick = get_gravkick_factor(time0, time1);
      dt_hydrokick = get_hydrokick_factor(time0, time1);
    }
  else
    {
      dt_drift = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	P[i].Pos[j] += P[i].Vel[j] * dt_drift;

      if(P[i].Type == 0)
	{
#ifdef PMGRID
	  for(j = 0; j < 3; j++)
	    SphP[i].VelPred[j] +=
	      (P[i].GravAccel[j] + P[i].GravPM[j]) * dt_gravkick + SphP[i].HydroAccel[j] * dt_hydrokick;
#else
	  for(j = 0; j < 3; j++)
	    SphP[i].VelPred[j] += P[i].GravAccel[j] * dt_gravkick + SphP[i].HydroAccel[j] * dt_hydrokick;
#endif
	  SphP[i].Density *= exp(-SphP[i].DivVel * dt_drift);
	  SphP[i].Hsml *= exp(0.333333333333 * SphP[i].DivVel * dt_drift);

	  if(SphP[i].Hsml < All.MinGasHsml)
	    SphP[i].Hsml = All.MinGasHsml;

	  dt_entr = (time1 - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

	  SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
    
#ifdef VARIABLE_VISC_CONST
    SphP[i].Alpha += SphP[i].DtAlpha * dt_drift;
#endif
	}
    }

  /* if domain-decomp and tree are not going to be reconstructed, update dynamically.  */
  if(All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
    {
      for(i = 0; i < Numnodestree; i++)
	for(j = 0; j < 3; j++)
	  Nodes[All.MaxPart + i].u.d.s[j] += Extnodes[All.MaxPart + i].vs[j] * dt_drift;

      force_update_len();

      force_update_pseudoparticles();
    }

  t1 = second();

  All.CPU_Predict += timediff(t0, t1);
}



/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
  int i, j;
  double boxsize[3];

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

#ifdef LONG_X
  boxsize[0] *= LONG_X;
#endif
#ifdef LONG_Y
  boxsize[1] *= LONG_Y;
#endif
#ifdef LONG_Z
  boxsize[2] *= LONG_Z;
#endif

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
	while(P[i].Pos[j] < 0)
	  P[i].Pos[j] += boxsize[j];

	while(P[i].Pos[j] >= boxsize[j])
	  P[i].Pos[j] -= boxsize[j];
      }
}
#endif




#ifdef SINK_PARTICLES
void identify_doomed_particles(void)
{
  int n, i, j, k, pindex;
  float seperation, relvel, relenergy, little_L, KeplerL2, sinkrad;
  float EffectiveCutoff;
  int num, startnode;
  int numsinks, numsinkstot;
  int notestflag;
  double *local_sink_posx, *local_sink_posy, *local_sink_posz;
  double *local_sink_velx, *local_sink_vely, *local_sink_velz;
  double *local_sink_mass;  
  int *local_sink_ID;  
  int verbose = 0;
  FLOAT *pos, *vel;
  FLOAT Postemp[3], Veltemp[3];
  
  printf("starting accretion, rank %d, %d accretors\n",ThisTask,NumPart-N_gas);
  AccNum = 0;
    
  numsinks = NumPart - N_gas;
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&numsinks, &numsinkstot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
  local_sink_posx = malloc(sizeof(double) * numsinkstot);
  local_sink_posy = malloc(sizeof(double) * numsinkstot);  
  local_sink_posz = malloc(sizeof(double) * numsinkstot); 
  local_sink_velx = malloc(sizeof(double) * numsinkstot);
  local_sink_vely = malloc(sizeof(double) * numsinkstot);  
  local_sink_velz = malloc(sizeof(double) * numsinkstot);  
  local_sink_ID = malloc(sizeof(int) * numsinkstot);   
  local_sink_mass = malloc(sizeof(double) * numsinkstot);     
  list_sink_posx = malloc(sizeof(double) * numsinkstot * NTask);
  list_sink_posy = malloc(sizeof(double) * numsinkstot * NTask);  
  list_sink_posz = malloc(sizeof(double) * numsinkstot * NTask); 
  list_sink_velx = malloc(sizeof(double) * numsinkstot * NTask);
  list_sink_vely = malloc(sizeof(double) * numsinkstot * NTask);  
  list_sink_velz = malloc(sizeof(double) * numsinkstot * NTask); 
  list_sink_ID = malloc(sizeof(int) * numsinkstot * NTask);    
  list_sink_mass = malloc(sizeof(double) * numsinkstot * NTask);  
  
  
  for(i = 0; i < numsinkstot; i++) local_sink_mass[i] = -1;
  
  for(i = 0; i < numsinks; i++){
    local_sink_posx[i] = P[i+N_gas].Pos[0];
    local_sink_posy[i] = P[i+N_gas].Pos[1];    
    local_sink_posz[i] = P[i+N_gas].Pos[2];
    local_sink_velx[i] = P[i+N_gas].Vel[0];
    local_sink_vely[i] = P[i+N_gas].Vel[1];    
    local_sink_velz[i] = P[i+N_gas].Vel[2];
    local_sink_ID[i] = P[i+N_gas].ID;    
    local_sink_mass[i] = P[i+N_gas].Mass;     
  }  
  /* MPI_Barrier(MPI_COMM_WORLD);
   for(i = 0; i < numsinkstot; i++) printf("loc identify: %d %d %d %f %f\n",ThisTask,numsinks,NTask*numsinkstot,local_sink_posx[i],local_sink_mass[i]);
   */    
  MPI_Barrier(MPI_COMM_WORLD);   
  MPI_Allgather(local_sink_posx, numsinkstot, MPI_DOUBLE, list_sink_posx, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_posy, numsinkstot, MPI_DOUBLE, list_sink_posy, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_posz, numsinkstot, MPI_DOUBLE, list_sink_posz, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_velx, numsinkstot, MPI_DOUBLE, list_sink_velx, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_vely, numsinkstot, MPI_DOUBLE, list_sink_vely, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_velz, numsinkstot, MPI_DOUBLE, list_sink_velz, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_mass, numsinkstot, MPI_DOUBLE, list_sink_mass, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_ID, numsinkstot, MPI_INT, list_sink_ID, numsinkstot, MPI_INT, MPI_COMM_WORLD);           
  MPI_Barrier(MPI_COMM_WORLD); 
  
  for(i = 0; i < numsinkstot * NTask; i++) if(ThisTask == 0 && list_sink_mass[i] > 0) fprintf(FdInfo,"  %d:%e  ",list_sink_ID[i],list_sink_mass[i]);      
  if(ThisTask == 0) fprintf(FdInfo,"\n");
  for(i = 0; i < numsinkstot * NTask; i++) if(ThisTask == 0 && list_sink_mass[i] > 0) fprintf(FdInfo,"  %e:%e  ",list_sink_posx[i],list_sink_velx[i]);      
  if(ThisTask == 0) fprintf(FdInfo,"\n");  
  fflush(FdInfo);
  
  for(i = 0; i < numsinkstot * NTask; i++){ /* go through all the sink particles (From all processors) and find doomed gas */
    notestflag = 0;
    if(list_sink_mass[i] > 0){
      Postemp[0] = list_sink_posx[i];
      Postemp[1] = list_sink_posy[i];
      Postemp[2] = list_sink_posz[i];
      Veltemp[0] = list_sink_velx[i];
      Veltemp[1] = list_sink_vely[i];
      Veltemp[2] = list_sink_velz[i];            
      pos = Postemp;
      vel = Veltemp;
      
      startnode = All.MaxPart;
      sinkrad = All.AccretionRadius;
      
#ifdef NO_ACC_TEST
      notestflag = 1;
#endif             
      
      
      KeplerL2 = All.G * list_sink_mass[i] * sinkrad;
      num = ngb_treefind_variable(&pos[0],sinkrad,&startnode); /* find all particles inside the sink radius */
      
      for(n = 0; n < num; n++){
        k = Ngblist[n];
        
        if(P[k].Type == 0 && k < N_gas){  /* only accrete gas particles! */
          for(seperation = 0,j = 0; j < 3; j++) seperation += (P[k].Pos[j]-pos[j]) * (P[k].Pos[j]-pos[j]);  /* r.r */  
          seperation = sqrt(seperation);   /* r */
          
          if(seperation < sinkrad){
            for(relvel = 0,j = 0; j < 3; j++)
              relvel += (P[k].Vel[j]-vel[j]) * (P[k].Vel[j]-vel[j]);      /* v.v */
            
            
            relenergy = .5*relvel - All.G*list_sink_mass[i]/seperation;
            
            if(notestflag) relenergy = -1;  
            if(relenergy < 0){
              for(little_L = 0, j = 0;j < 3; j++) little_L +=pow( (P[k].Pos[(j+1)%3] - pos[(j+1)%3]) * (P[k].Vel[(j+2)%3] - vel[(j+2)%3]) -  \
                                                                 (P[k].Pos[(j+2)%3] - pos[(j+2)%3]) * (P[k].Vel[(j+1)%3] - vel[(j+1)%3]) ,2); /* L.L */
              if(notestflag) little_L = KeplerL2 /2;								  						
              if(little_L < KeplerL2){             
                if(SphP[k].AccretionTarget == 0){  /* if targeted by another sink, it doesn't get accreted */ 
                  SphP[k].AccretionTarget = list_sink_ID[i];       /* if bound in E and L, accrete it */
                  /*AccreteList[AccNum] = k;*/
                  if(verbose) printf("Particle ID %d provisionally scheduled for destruction onto sink ID %d from %d %f %f\n",P[k].ID,list_sink_ID[i],ThisTask,seperation,sinkrad);
                  /*AccNum++; */ 
                }
/* let it go until it's unambiguously targeted */                
                else{
                  SphP[k].AccretionTarget = -2; /* if targeted by multiple sinks, it doesn't get accreted by any*/
                  if(verbose) printf("%d targeted twice! sink ID %d from %d %f %f\n",P[k].ID,list_sink_ID[i],ThisTask,seperation,sinkrad);
                }
              }  
            }
          }
        }	       	      	     
      } 
    }
  }
  
  if(verbose) printf("Confirmed for accretion: ");
  for(i = 0; i < N_gas; i++){
    if(SphP[i].AccretionTarget > 0){
      AccreteList[AccNum] = i;
      AccNum++;
      if(verbose) printf("%d ",P[i].ID);      	
    }
  } 
  if(verbose) printf("\n");
  
  
#ifdef CUTOFF_RADIUS  
  for(pindex = -1, j = N_gas; j < NumPart; j++) 
    if(P[j].Mass > All.MassSeg) pindex = j;
  for(i = 0; i < N_gas; i++){
    for(seperation = 0, j = 0; j < 3; j++) seperation += P[i].Pos[j] * P[i].Pos[j];
    if(sqrt(seperation) > EffectiveCutoff){
      SphP[i].AccretionTarget = -1; /* don't want to add mass during 'accretion' */
      AccreteList[AccNum] = i;
      AccNum++;
      if(verbose) printf("Particle ID %d (%d) scheduled for ignoring. acctarg %d \n",P[i].ID,i,SphP[i].AccretionTarget); 
    }
  }
#endif
  
  /*  printf("rank %d Accnum: %d\n",ThisTask,AccNum); 
   */
  All.TstepLastAcc = All.NumCurrentTiStep;
  free(list_sink_posx);
  free(list_sink_posy);
  free(list_sink_posz);  
  free(list_sink_velx);
  free(list_sink_vely);
  free(list_sink_velz);    
  free(list_sink_ID);
  free(list_sink_mass);  
  free(local_sink_posx);
  free(local_sink_posy);
  free(local_sink_posz);  
  free(local_sink_velx);
  free(local_sink_vely);
  free(local_sink_velz); 
  free(local_sink_ID);      
  free(local_sink_mass);      
}



void destroy_doomed_particles(void)
{
  int n, i, j, k, s, target, acc_counter, distant_sink_counter, accflag = 0;
  int numsinks, numsinkstot;
  double combined_mass; 
  double *local_sink_posx, *local_sink_posy, *local_sink_posz;
  double *local_sink_velx, *local_sink_vely, *local_sink_velz;
  double *local_sink_mass;
  int *local_sink_ID;
  double dposx, dposy, dposz, dvelx, dvely, dvelz, dmass;
  double dposxtot, dposytot, dposztot, dvelxtot, dvelytot, dvelztot, dmasstot;      
  double Ei, Ef, Einit, Efinal, sep, hinv, u, wp;
  double entropychange,dt_entropy,dt_grav;
  double momx,momy,momz;

  
  for(k = N_gas; k < NumPart; k++) printf("ID %d (%d) init vel, pos, mass: (%e|%e|%e), (%e|%e|%e), %e\n",P[k].ID,k,P[k].Vel[0],P[k].Vel[1],P[k].Vel[2],
                                          P[k].Pos[0],P[k].Pos[1],P[k].Pos[2],P[k].Mass); 
  /* first transfer momentum from all particles that are going to be accreted */
  
  numsinks = NumPart - N_gas;  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&numsinks, &numsinkstot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
  
  local_sink_posx = malloc(sizeof(double) * numsinkstot);
  local_sink_posy = malloc(sizeof(double) * numsinkstot);  
  local_sink_posz = malloc(sizeof(double) * numsinkstot); 
  local_sink_velx = malloc(sizeof(double) * numsinkstot);
  local_sink_vely = malloc(sizeof(double) * numsinkstot);  
  local_sink_velz = malloc(sizeof(double) * numsinkstot);  
  local_sink_ID = malloc(sizeof(int) * numsinkstot);   
  local_sink_mass = malloc(sizeof(double) * numsinkstot);     
  list_sink_posx = malloc(sizeof(double) * numsinkstot * NTask);
  list_sink_posy = malloc(sizeof(double) * numsinkstot * NTask);  
  list_sink_posz = malloc(sizeof(double) * numsinkstot * NTask); 
  list_sink_velx = malloc(sizeof(double) * numsinkstot * NTask);
  list_sink_vely = malloc(sizeof(double) * numsinkstot * NTask);  
  list_sink_velz = malloc(sizeof(double) * numsinkstot * NTask); 
  list_sink_ID = malloc(sizeof(int) * numsinkstot * NTask);    
  list_sink_mass = malloc(sizeof(double) * numsinkstot * NTask);
  
  for(i = 0; i < numsinkstot; i++) local_sink_mass[i] = -1;
  
  for(i = 0; i < numsinks; i++){
    dt_grav = All.Timebase_interval * (All.Ti_Current - (P[i+N_gas].Ti_begstep + P[i+N_gas].Ti_endstep) / 2);
    local_sink_posx[i] = P[i+N_gas].Pos[0];
    local_sink_posy[i] = P[i+N_gas].Pos[1];    
    local_sink_posz[i] = P[i+N_gas].Pos[2];
    local_sink_velx[i] = P[i+N_gas].Vel[0] + dt_grav * P[i+N_gas].GravAccel[0];
    local_sink_vely[i] = P[i+N_gas].Vel[1] + dt_grav * P[i+N_gas].GravAccel[1];    
    local_sink_velz[i] = P[i+N_gas].Vel[2] + dt_grav * P[i+N_gas].GravAccel[2];
    local_sink_ID[i] = P[i+N_gas].ID;    
    local_sink_mass[i] = P[i+N_gas].Mass;      
  }  
  
  /* for(i = 0; i < numsinkstot; i++) printf("loc destroy: %d %d %d %f %f\n",ThisTask,numsinks,numsinkstot,local_sink_posx[i],local_sink_mass[i]);
   */ 
  MPI_Barrier(MPI_COMM_WORLD);   
  MPI_Allgather(local_sink_posx, numsinkstot, MPI_DOUBLE, list_sink_posx, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_posy, numsinkstot, MPI_DOUBLE, list_sink_posy, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_posz, numsinkstot, MPI_DOUBLE, list_sink_posz, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_velx, numsinkstot, MPI_DOUBLE, list_sink_velx, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_vely, numsinkstot, MPI_DOUBLE, list_sink_vely, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_velz, numsinkstot, MPI_DOUBLE, list_sink_velz, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);  
  MPI_Allgather(local_sink_mass, numsinkstot, MPI_DOUBLE, list_sink_mass, numsinkstot, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(local_sink_ID, numsinkstot, MPI_INT, list_sink_ID, numsinkstot, MPI_INT, MPI_COMM_WORLD);           
  MPI_Barrier(MPI_COMM_WORLD);
  
  Ei = 0;
  Ef = 0;
  
  /************ 
   for(j = N_gas;j < NumPart; j++)
   Ei += .5 * P[j].Mass * (P[j].Vel[0]*P[j].Vel[0] + P[j].Vel[1]*P[j].Vel[1] + P[j].Vel[2]*P[j].Vel[2]);
   ************/  
  
  for(s = 0; s < numsinkstot*NTask; s++){  
    if(list_sink_mass[s] > 0){
      MPI_Barrier(MPI_COMM_WORLD);
      dvelx = 0; dvely = 0; dvelz = 0; dposx = 0; dposy = 0; dposz = 0; dmass = 0;      
      target = list_sink_ID[s];
      
      for(k = 0;k < AccNum;k++){
        i = AccreteList[k];
        if(SphP[i].AccretionTarget == target){
          accflag = 1;
          
          dvelx += P[i].Mass*P[i].Vel[0];
          dvely += P[i].Mass*P[i].Vel[1];	    
          dvelz += P[i].Mass*P[i].Vel[2];	
          dposx += P[i].Mass*P[i].Pos[0];
          dposy += P[i].Mass*P[i].Pos[1];	    
          dposz += P[i].Mass*P[i].Pos[2];   
          dmass += P[i].Mass;  
        }
      } /* now we have all particles on the local processor that add to sink s */  
      
      /* accumulate the changes from all processors */
      dvelxtot = 0; dvelytot = 0; dvelztot = 0; dposxtot = 0; dposytot = 0; dposztot = 0; dmasstot = 0; 
      MPI_Barrier(MPI_COMM_WORLD); 
      
      MPI_Allreduce(&dvelx, &dvelxtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
      MPI_Allreduce(&dvely, &dvelytot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
      MPI_Allreduce(&dvelz, &dvelztot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
      MPI_Allreduce(&dposx, &dposxtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
      MPI_Allreduce(&dposy, &dposytot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
      MPI_Allreduce(&dposz, &dposztot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);     
      MPI_Allreduce(&dmass, &dmasstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);        
      MPI_Barrier(MPI_COMM_WORLD);     
      
      /* check to see if the sink being considered is on the local processor.  if so, add the changes to it. */
      for(j = N_gas;j < NumPart; j++){
        if(P[j].ID == target){
          printf("jay = %d\n",j);
          
          dposxtot += P[j].Pos[0] * P[j].Mass;
          dposytot += P[j].Pos[1] * P[j].Mass;	
          dposztot += P[j].Pos[2] * P[j].Mass;
          dmasstot += P[j].Mass;
          
          P[j].Pos[0] = dposxtot / dmasstot;
          P[j].Pos[1] = dposytot / dmasstot;	  	  
          P[j].Pos[2] = dposztot / dmasstot;	  
          
          dt_grav = All.Timebase_interval * (All.Ti_Current - (P[j].Ti_begstep + P[j].Ti_endstep) / 2);
          dt_grav=0;
          dvelxtot += P[j].Mass * (P[j].Vel[0] + dt_grav * P[j].GravAccel[0]);	  
          dvelytot += P[j].Mass * (P[j].Vel[1] + dt_grav * P[j].GravAccel[1]);	  
          dvelztot += P[j].Mass * (P[j].Vel[2] + dt_grav * P[j].GravAccel[2]);
          
          P[j].Vel[0] = dvelxtot / dmasstot - dt_grav * P[j].GravAccel[0];
          P[j].Vel[1] = dvelytot / dmasstot - dt_grav * P[j].GravAccel[1];
          P[j].Vel[2] = dvelztot / dmasstot - dt_grav * P[j].GravAccel[2];	  	  
          
          P[j].Mass = dmasstot;
          printf("ID %d task %d accnum %d final vel, pos, mass: (%e|%e|%e), (%e|%e|%e), %e\n", \
                 P[j].ID,ThisTask,AccNum,P[j].Vel[0],P[j].Vel[1],P[j].Vel[2], \
                 P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],P[j].Mass);
          P[j].Ti_endstep = All.Ti_Current;
          
        }   
      }
    }        
  }
  MPI_Barrier(MPI_COMM_WORLD); 

  if(AccNum > 1) qsort(AccreteList, AccNum, sizeof(int), index_compare_key);
  
  acc_counter = 0;
  for(n = 0;n < AccNum;n++){
    
    i = AccreteList[n] - acc_counter;
    if(SphP[i].AccretionTarget > -2){      
      if(P[i].Ti_endstep == All.Ti_Current){
        NumForceUpdate--;
        NumSphUpdate--;
      }
      for(k = i+1; k<=NumPart; k++){ // actually remove the particle here, 
                                     // and shift everything down to fill the gap in the array.
        P[k-1] = P[k];
        if(P[k].Type == 0)
          SphP[k-1] = SphP[k];      
      }
      NumPart--;   // decrement the local countrs of particles and gas particles
      N_gas--;
      acc_counter++; 
    }
  }
    
  free(list_sink_posx);
  free(list_sink_posy);
  free(list_sink_posz);  
  free(list_sink_velx);
  free(list_sink_vely);
  free(list_sink_velz);    
  free(list_sink_ID);
  free(list_sink_mass);  
  free(local_sink_posx);
  free(local_sink_posy);
  free(local_sink_posz);  
  free(local_sink_velx);
  free(local_sink_vely);
  free(local_sink_velz); 
  free(local_sink_ID);      
  free(local_sink_mass); 
  
  AccNum = 0;
  
  printf("done with accretion, rank %d\n",ThisTask); 
}

int index_compare_key(const void *a, const void *b)
{
  return ( *(int*)a - *(int*)b );
}

#endif /* SINK_PARTICLES */

