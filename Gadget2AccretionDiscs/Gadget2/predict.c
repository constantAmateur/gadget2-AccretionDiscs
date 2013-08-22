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
#ifdef DEAD_GAS
      double r;
      r=sqrt(P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+P[i].Pos[2]*P[i].Pos[2]);
      if(r>All.Doom_radius)
      {
        //printf("Damn! Particle at (%g,%g,%g) r = %g is dead.\n",P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],r);
        for(j=0;j<3;j++)
        {
          P[i].Pos[j] -= All.Drift_speed*(P[i].Pos[j]/r)*dt_drift;
        }
      }
      else
      {
        //printf("Hooray! Particle at (%g,%g,%g) r = %g is alive.\n",P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],r);
        for(j=0;j<3;j++)
        {
	       P[i].Pos[j] += P[i].Vel[j] * dt_drift;
        }
      }
#else
      for(j = 0; j < 3; j++)
	P[i].Pos[j] += P[i].Vel[j] * dt_drift;
#endif
      if(P[i].Type==1)
      {
        printf("Drifting particle %d with vel %g,%g,%g.\n",i,P[i].Vel[0],P[i].Vel[1],P[i].Vel[2]);
        printf("New pos is %g,%g,%g.\n",P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);
      }

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
#ifdef SURFACE
     SphP[i].SurDensity *= exp(-SphP[i].SurDivVel * dt_drift);
     SphP[i].SurHsml *= exp(0.333333333333 * SphP[i].SurDivVel * dt_drift);
#endif

	  if(SphP[i].Hsml < All.MinGasHsml)
	    SphP[i].Hsml = All.MinGasHsml;

	  dt_entr = (time1 - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

     //The value of Entropy in the brackets here is "predicted entropy"
	  SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
   
#if defined MMAV || defined CDAV
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
#ifdef CUTOFF_RADIUS
  int pindex;
  float EffectiveCutoff;
#endif
  int n, i, j, k;
  float seperation, relvel, relenergy, little_L, KeplerL2, sinkrad;
  int num, startnode;
  int numsinks, numsinkstot;
  int notestflag;
  double *local_sink_posx, *local_sink_posy, *local_sink_posz;
  double *local_sink_velx, *local_sink_vely, *local_sink_velz;
  double *local_sink_mass;  
  int *local_sink_ID;  
#ifdef VERBOSE
  int verbose = 1;
#else
  int verbose = 0;
#endif
  FLOAT *pos, *vel;
  FLOAT Postemp[3], Veltemp[3];
  
  //printf("starting accretion, rank %d, %d accretors\n",ThisTask,NumPart-N_gas);
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
  
  //for(i = 0; i < numsinkstot * NTask; i++) if(ThisTask == 0 && list_sink_mass[i] > 0) fprintf(FdInfo,"  %d:%e  ",list_sink_ID[i],list_sink_mass[i]);      
  //if(ThisTask == 0) fprintf(FdInfo,"\n");
  //for(i = 0; i < numsinkstot * NTask; i++) if(ThisTask == 0 && list_sink_mass[i] > 0) fprintf(FdInfo,"  %e:%e  ",list_sink_posx[i],list_sink_velx[i]);      
  //if(ThisTask == 0) fprintf(FdInfo,"\n");  
  //fflush(FdInfo);
  
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
      //This needs to be fixed up to handle the case where the ngb buffer gets full
      //i.e., stick inside do while
      //do
      //{
      //}while(startnode>=0);
      do
      {
        num = ngb_treefind_variable(&pos[0],sinkrad,&startnode); /* find all particles inside the sink radius */
        
        for(n = 0; n < num; n++){
          k = Ngblist[n];
          
          //We want to only mark particles for accretion if they haven't been marked previously
          if(P[k].Type == 0 && k < N_gas && (SphP[k].AccretionTarget ==0 || SphP[k].AccretionTarget==-2)){  /* only accrete gas particles! */
            for(seperation = 0,j = 0; j < 3; j++) seperation += (P[k].Pos[j]-pos[j]) * (P[k].Pos[j]-pos[j]);  /* r.r */  
            seperation = sqrt(seperation);   /* r */
            
            if(seperation < sinkrad){

              if(verbose) printf("Particle ID %d is within accretion radius of sink ID %d from %d %f %f\n",P[k].ID,list_sink_ID[i],ThisTask,seperation,sinkrad);
              for(relvel = 0,j = 0; j < 3; j++)
                relvel += (P[k].Vel[j]-vel[j]) * (P[k].Vel[j]-vel[j]);      /* v.v */
              
              
              relenergy = .5*relvel - All.G*list_sink_mass[i]/seperation;
              
              if(notestflag) relenergy = -1;  
              if(relenergy < 0){
                for(little_L = 0, j = 0;j < 3; j++) little_L +=pow( (P[k].Pos[(j+1)%3] - pos[(j+1)%3]) * (P[k].Vel[(j+2)%3] - vel[(j+2)%3]) -  (P[k].Pos[(j+2)%3] - pos[(j+2)%3]) * (P[k].Vel[(j+1)%3] - vel[(j+1)%3]) ,2); /* L.L */
                if(notestflag) little_L = KeplerL2 /2;								  						
                //This test may be a bit too stringent, turn it off if things aren't being accreted
                if(little_L < KeplerL2){             
                  if(SphP[k].AccretionTarget == 0){  /* if targeted by another sink, it doesn't get accreted */ 
                    SphP[k].AccretionTarget = list_sink_ID[i];       /* if bound in E and L, accrete it */
                    /*AccreteList[AccNum] = k;*/
                    if(verbose) printf("Particle ID %d provisionally scheduled for destruction onto sink ID %d from %d %f %f\n",P[k].ID,list_sink_ID[i],ThisTask,seperation,sinkrad);
                    /*AccNum++; */ 
                  }
/* let   it go until it's unambiguously targeted */                
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
      while(startnode>=0);
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
  //Delete's particle that are beyond a certain radius
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
  int n, i, j, k, s, target, acc_counter, accflag = 0;
  int numsinks, numsinkstot;
  double *local_sink_posx, *local_sink_posy, *local_sink_posz;
  double *local_sink_velx, *local_sink_vely, *local_sink_velz;
  double *local_sink_mass;
  int *local_sink_ID;
  double dposx, dposy, dposz, dvelx, dvely, dvelz, dmass;
  double dposxtot, dposytot, dposztot, dvelxtot, dvelytot, dvelztot, dmasstot;      
  double Ei, Ef;
  double dt_grav;
#ifdef TRACK_ACCRETION_LOSSES
  double starR[3],starv[3],starM,starRtot[3],starvtot[3],starMtot;
  double accretion_int,accretion_rad,accretion_kin,accretion_pot;
  double accretion_angmom[3];
  double acc_int_tot,acc_rad_tot,acc_kin_tot,acc_pot_tot,acc_angmom_tot[3];
  double com[7],delR[3],delv[3];
  int *list_acc_num;
  int accnumtot,start,stop,ii,accflagtot,tmpflag;
#ifdef HIGH_PRECISION_POT
  double acc_pot_start,acc_pot_temp,acc_pot_finish;
#endif
#else
#ifdef CUTOFF_BOX
  double starR[3],starRtot[3];
#endif
#endif

  
  //for(k = N_gas; k < NumPart; k++) printf("ID %d (%d) init vel, pos, mass: (%e|%e|%e), (%e|%e|%e), %e\n",P[k].ID,k,P[k].Vel[0],P[k].Vel[1],P[k].Vel[2],
  //                                        P[k].Pos[0],P[k].Pos[1],P[k].Pos[2],P[k].Mass); 
  /* first transfer momentum from all particles that are going to be accreted */
  
  numsinks = NumPart - N_gas;  
  //Could we replace this with NtypeLocal[1]?
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
#ifdef TRACK_ACCRETION_LOSSES
  list_acc_num = malloc(sizeof(int) * NTask);
#endif
  
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
#ifdef TRACK_ACCRETION_LOSSES
      //Need the initial stars position, mass and velocity
      starR[0]=starR[1]=starR[2]=starv[0]=starv[1]=starv[2]=starM=0;
      for(j = N_gas;j < NumPart; j++){
        if(P[j].ID == target){
          starR[0] = P[j].Pos[0];
          starR[1] = P[j].Pos[1];
          starR[2] = P[j].Pos[2];

          dt_grav = All.Timebase_interval * (All.Ti_Current - (P[j].Ti_begstep + P[j].Ti_endstep) / 2);
          dt_grav=0;
          starv[0] = (P[j].Vel[0] + dt_grav * P[j].GravAccel[0]);	  
          starv[1] = (P[j].Vel[1] + dt_grav * P[j].GravAccel[1]);	  
          starv[2] = (P[j].Vel[2] + dt_grav * P[j].GravAccel[2]);	  

          starM= P[j].Mass;
        }
      }
      MPI_Barrier(MPI_COMM_WORLD); 
      MPI_Allreduce(&starR[0], &starRtot[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
      MPI_Allreduce(&starv[0], &starvtot[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
      MPI_Allreduce(&starM, &starMtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
      MPI_Allreduce(&AccNum, &accnumtot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
      MPI_Allgather(&AccNum, 1, MPI_INT, list_acc_num, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);     
     
      //Now loop over all accreting particles
      accretion_int=accretion_rad=accretion_kin=accretion_pot=accretion_angmom[0]=accretion_angmom[1]=accretion_angmom[2]=0;
#ifdef HIGH_PRECISION_POT
      //Tracking the potential is challenging...
      MPI_Barrier(MPI_COMM_WORLD);
      printf("Calculating potential before accretion.\n");
#ifndef ADD_CENTRAL_GRAVITY
      TreeReconstructFlag =1;
#endif
      compute_potential();
      acc_pot_temp =0;
      for(j=0;j<NumPart;j++){
#ifdef ADD_CENTRAL_GRAVITY
        acc_pot_temp += P[j].Mass * P[j].Potential;
#else
        acc_pot_temp += 0.5 * P[j].Mass * P[j].Potential;
#endif
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(&acc_pot_temp,&acc_pot_start,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
      start=0;
      for(k=0;k<ThisTask;k++)
        start+=list_acc_num[k];
      stop=start+list_acc_num[ThisTask];
      for(k=0;k<accnumtot;k++){
        tmpflag=-1;
        com[0]=com[1]=com[2]=com[3]=com[4]=com[5]=com[6]=0;
        //The particle we want is on this processor!
        if(k>=start && k<stop){
          i = AccreteList[k-start];
          //And it's for this star!
          if(SphP[i].AccretionTarget == target){
            tmpflag = ThisTask;

            if(P[i].Type==0){
              accretion_int -= (P[i].Mass*SphP[i].Entropy * pow(SphP[i].Density,(GAMMA_MINUS1)))/(GAMMA_MINUS1);
#if defined BETA_COOLING && defined EXTRA_STATS
              accretion_rad -= P[i].Mass*SphP[i].RadiatedEnergy;
#endif
            }
            accretion_kin -= ((P[i].Mass*starMtot)/(2.0*(starMtot+P[i].Mass)))*((P[i].Vel[0]-starvtot[0])*(P[i].Vel[0]-starvtot[0])+(P[i].Vel[1]-starvtot[1])*(P[i].Vel[1]-starvtot[1])+(P[i].Vel[2]-starvtot[2])*(P[i].Vel[2]-starvtot[2]));
            delR[0]=starRtot[0]-P[i].Pos[0];
            delR[1]=starRtot[1]-P[i].Pos[1];
            delR[2]=starRtot[2]-P[i].Pos[2];
            delv[0]=starvtot[0]-P[i].Vel[0];
            delv[1]=starvtot[1]-P[i].Vel[1];
            delv[2]=starvtot[2]-P[i].Vel[2];
            com[0]=(starRtot[0]*starMtot+P[i].Pos[0]*P[i].Mass)/(P[i].Mass+starMtot);
            com[1]=(starRtot[1]*starMtot+P[i].Pos[2]*P[i].Mass)/(P[i].Mass+starMtot);
            com[2]=(starRtot[2]*starMtot+P[i].Pos[2]*P[i].Mass)/(P[i].Mass+starMtot);
            com[3]=P[i].Mass;
            com[4]=(P[i].Mass*P[i].Vel[0]+starMtot*starvtot[0])/(P[i].Mass+starMtot);
            com[5]=(P[i].Mass*P[i].Vel[1]+starMtot*starvtot[1])/(P[i].Mass+starMtot);
            com[6]=(P[i].Mass*P[i].Vel[2]+starMtot*starvtot[2])/(P[i].Mass+starMtot);
            accretion_angmom[0] -= ((starMtot*P[i].Mass)/(starMtot+P[i].Mass))*(delR[1]*delv[2]-delR[2]*delv[1]);
            accretion_angmom[1] -= ((starMtot*P[i].Mass)/(starMtot+P[i].Mass))*(delR[2]*delv[0]-delR[0]*delv[2]);
            accretion_angmom[2] -= ((starMtot*P[i].Mass)/(starMtot+P[i].Mass))*(delR[0]*delv[1]-delR[1]*delv[0]);
            //This is going to be slightly inaccurate because we potentially(hehe) haven't calculated the potential in a while, so it's out of date
#ifndef HIGH_PRECISION_POT
#ifdef ADD_CENTRAL_GRAVITY
            accretion_pot += (All.G*P[i].Mass*starMtot)/sqrt((P[i].Pos[0]-starRtot[0])*(P[i].Pos[0]-starRtot[0])+(P[i].Pos[1]-starRtot[1])*(P[i].Pos[1]-starRtot[1])+(P[i].Pos[2]-starRtot[2])*(P[i].Pos[2]-starRtot[2]));
#else
            accretion_pot -= P[i].Mass*P[i].Potential;
#endif
            accretion_pot += (All.G*P[i].Mass*P[i].Mass)/sqrt((P[i].Pos[0]-starRtot[0])*(P[i].Pos[0]-starRtot[0])+(P[i].Pos[1]-starRtot[1])*(P[i].Pos[1]-starRtot[1])+(P[i].Pos[2]-starRtot[2])*(P[i].Pos[2]-starRtot[2]));

#endif
            //printf("Accretion potential is %g\n",accretion_pot);
          }
        }
        MPI_Barrier(MPI_COMM_WORLD);     
        MPI_Allreduce(&tmpflag, &accflagtot, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
        //Did anything accrete?
        if(accflagtot>=0)
        {
          MPI_Barrier(MPI_COMM_WORLD);
          MPI_Bcast(&com,7,MPI_DOUBLE,accflagtot,MPI_COMM_WORLD);
          MPI_Barrier(MPI_COMM_WORLD);
#ifndef HIGH_PRECISION_POT
          for(ii=0;ii < NumPart;ii++)
          {
            //We don't want to count the accreting particle, or the star.  The ID condition takes care of the star, the other one ensures that if we're on the same processor as the accreting particle, we won't process it.
#ifndef ADD_CENTRAL_GRAVITY
            if(P[ii].ID==target)
            {
              accretion_pot += P[ii].Potential * com[3];
            }
#else
            //if((ThisTask!=accflagtot || ii!=i) && (P[ii].ID!=target))
            if(P[ii].ID!=target)
            {
              accretion_pot -= (All.G*P[ii].Mass*com[3])/sqrt((P[ii].Pos[0]-starRtot[0])*(P[ii].Pos[0]-starRtot[0])+(P[ii].Pos[1]-starRtot[1])*(P[ii].Pos[1]-starRtot[1])+(P[ii].Pos[2]-starRtot[2])*(P[ii].Pos[2]-starRtot[2]));
            }
#endif
          }
#endif
          //Update the star's position, mass and velocity for the next particle...
          starRtot[0]=com[0];
          starRtot[1]=com[1];
          starRtot[2]=com[2];
          starMtot+=com[3];
          starvtot[0]=com[4];
          starvtot[1]=com[5];
          starvtot[2]=com[6];
        }
      }//End the loop over accreting particles...
      //Gather and store the total changes...
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(&accretion_int,&acc_int_tot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&accretion_rad,&acc_rad_tot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&accretion_kin,&acc_kin_tot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&accretion_pot,&acc_pot_tot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&accretion_angmom[0],&acc_angmom_tot[0],3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      All.Accretion_int += acc_int_tot;
      All.Accretion_rad += acc_rad_tot;
      All.Accretion_kin += acc_kin_tot;
      All.Accretion_pot += acc_pot_tot;
      All.Accretion_angmom[0] += acc_angmom_tot[0];
      All.Accretion_angmom[1] += acc_angmom_tot[1];
      All.Accretion_angmom[2] += acc_angmom_tot[2];
#else
#ifdef CUTOFF_BOX
      //Need to calculate star radius for box limiting
      starR[0]=starR[1]=starR[2]=0;
      for(j = N_gas;j < NumPart; j++){
        if(P[j].ID == target){
          starR[0] = P[j].Pos[0];
          starR[1] = P[j].Pos[1];
          starR[2] = P[j].Pos[2];
        }
      }
      MPI_Barrier(MPI_COMM_WORLD); 
      MPI_Allreduce(&starR[0], &starRtot[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
      MPI_Barrier(MPI_COMM_WORLD);     
#endif
#endif
      
      for(k = 0;k < AccNum;k++){
        i = AccreteList[k];
        if(SphP[i].AccretionTarget == target){
          accflag = 1;
          
          //Should we be using the predicted velocities here?
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
          //printf("jay = %d\n",j);
          
          dposxtot += P[j].Pos[0] * P[j].Mass;
          dposytot += P[j].Pos[1] * P[j].Mass;	
          dposztot += P[j].Pos[2] * P[j].Mass;
          dmasstot += P[j].Mass;
          
          //Move position to centre of mass
          P[j].Pos[0] = dposxtot / dmasstot;
          P[j].Pos[1] = dposytot / dmasstot;	  	  
          P[j].Pos[2] = dposztot / dmasstot;	  
          
          dt_grav = All.Timebase_interval * (All.Ti_Current - (P[j].Ti_begstep + P[j].Ti_endstep) / 2);
          dt_grav=0;
          //Should this be predicted velocity?
          dvelxtot += P[j].Mass * (P[j].Vel[0] + dt_grav * P[j].GravAccel[0]);	  
          dvelytot += P[j].Mass * (P[j].Vel[1] + dt_grav * P[j].GravAccel[1]);	  
          dvelztot += P[j].Mass * (P[j].Vel[2] + dt_grav * P[j].GravAccel[2]);
          
          //Add momentum to the sink
          P[j].Vel[0] = dvelxtot / dmasstot - dt_grav * P[j].GravAccel[0];
          P[j].Vel[1] = dvelytot / dmasstot - dt_grav * P[j].GravAccel[1];
          P[j].Vel[2] = dvelztot / dmasstot - dt_grav * P[j].GravAccel[2];	  	  
          
          //Add the mass to the sink
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
      //Need to decrement the Ntype and Ntypelocal variables too?
      acc_counter++; 
    }
  }
#ifdef CUTOFF_BOX
  //How far to go?
  i=N_gas;
  j=0;
  //Check each particle in array
  while(j<i)
  {
    //Is it outside the box?
    if((P[j].Pos[0]-starRtot[0])*(P[j].Pos[0]-starRtot[0])+(P[j].Pos[1]-starRtot[1])*(P[j].Pos[1]-starRtot[1]) > All.maxR2 || fabs(P[j].Pos[2]-starRtot[2]) > All.maxZ)
    {
      //It is, delete it.
      if(P[j].Ti_endstep == All.Ti_Current){
        NumForceUpdate--;
        NumSphUpdate--;
      }
      for(k = j+1; k<=NumPart; k++){ // actually remove the particle here, 
                                     // and shift everything down to fill the gap in the array.
        P[k-1] = P[k];
        if(P[k].Type == 0)
          SphP[k-1] = SphP[k];      
      }
      NumPart--;   // decrement the local countrs of particles and gas particles
      N_gas--;
      //There's one less particle to consider now, so stop loop a bit earlier
      i--;
    }
    j++;
  }
#endif
    
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
#ifdef TRACK_ACCRETION_LOSSES
  free(list_acc_num);
#ifdef HIGH_PRECISION_POT
  MPI_Barrier(MPI_COMM_WORLD); 
  printf("Calculating potential after accretion.\n");
#ifndef ADD_CENTRAL_GRAVITY
  TreeReconstructFlag =1;
#endif
  compute_potential();
  acc_pot_temp=0;
  for(i=0;i<NumPart;i++){
#ifdef ADD_CENTRAL_GRAVITY
    acc_pot_temp += P[i].Mass * P[i].Potential;
#else
    acc_pot_temp += 0.5 * P[i].Mass * P[i].Potential;
#endif
  }
  MPI_Barrier(MPI_COMM_WORLD); 
  MPI_Allreduce(&acc_pot_temp,&acc_pot_finish,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  All.Accretion_pot += (acc_pot_finish-acc_pot_start);
#endif
#endif
  
  AccNum = 0;
  
  //printf("done with accretion, rank %d\n",ThisTask); 
}

int index_compare_key(const void *a, const void *b)
{
  return ( *(int*)a - *(int*)b );
}

#endif /* SINK_PARTICLES */
