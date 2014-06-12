#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file init.c
 *  \brief Code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree. Various variables of the particle data are initialised and An intial
 *  domain decomposition is performed. If SPH particles are present, the inial
 *  SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;
  double a3;
#if defined BETA_COOLING || defined ADD_CENTRAL_GRAVITY 
  int starID,*list_starID;
  double starMass,*list_starMass;
#endif
#ifdef INJECT_GAS
  double binaryMasses[2];
  double binaryPositions[6];
  double *list_starPos,*list_starMass;
#endif

  All.Time = All.TimeBegin;
#ifdef SINK_PARTICLES
  All.AccreteFlag = 1;
#ifdef TRACK_ACCRETION_LOSSES
  All.Accretion_int = 0;
  All.Accretion_rad = 0;
  All.Accretion_kin = 0;
  All.Accretion_pot = 0;
  All.Accretion_angmom[0]=0;
  All.Accretion_angmom[1]=0;
  All.Accretion_angmom[2]=0;
#endif
  //Square radius for fast comparison
#ifdef CUTOFF_BOX
  All.maxR2 = All.maxR2*All.maxR2;
#endif
#endif

  switch (All.ICFormat)
    {
    case 1:
#if (MAKEGLASS > 1)
      seed_glass();
#else
      read_ic(All.InitCondFile);
#endif
      break;
    case 2:
    case 3:
      read_ic(All.InitCondFile);
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

#ifdef ACCRETED_MASS_ONLY
  for(i=N_gas;i<NumPart;i++)
    P[i].NumAccreted=0.0;
#endif
//#ifdef ACCRETED_MASS_ONLY
//  double m[2]={0};
//  double ms[2*NTask];
//  for(i=N_gas;i<NumPart;i++)
//  {
//    if(P[i].Mass>m[1])
//    {
//      m[0]=m[1];
//      m[1]=P[i].Mass;
//    }
//    else
//    {
//      if(P[i].Mass>m[0])
//        m[0]=P[i].Mass;
//    }
//  }
//  printf("Set m1=%g m2=%g.\n",m[1],m[0]);
//  MPI_Allgather(&m,2,MPI_DOUBLE,ms,2,MPI_DOUBLE,MPI_COMM_WORLD);
//  m[0]=m[1]=0;
//  for(i=0;i<2*NTask;i++)
//  {
//    if(ms[i]>m[1])
//    {
//      m[0]=m[1];
//      m[1]=ms[i];
//    }
//    else
//    {
//      if(ms[i]>m[0])
//        m[0]=ms[i];
//    }
//  }
//  printf("Final Set m1=%g m2=%g.\n",m[1],m[0]);
//  All.M1=m[1];
//  All.M2=m[0];
//#endif

  All.Time = All.TimeBegin;
  All.Ti_Current = 0;

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      a3 = All.Time * All.Time * All.Time;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      a3 = 1;
    }

  set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;

  All.TotNumOfForces = 0;
  All.NumForcesSinceLastDomainDecomp = 0;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }
#if defined BETA_COOLING || defined ADD_CENTRAL_GRAVITY
  starMass=-1.0;
  starID=-1;
#endif
  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;
#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_endstep = 0;
      P[i].Ti_begstep = 0;

      P[i].OldAcc = 0;
      P[i].GravCost = 1;
      P[i].Potential = 0;
#if defined BETA_COOLING || defined ADD_CENTRAL_GRAVITY
      /* All processors will have the ID with the largest mass in them */
      if(i>=N_gas)
      {
        if(P[i].Mass>starMass)
        {
          starMass=P[i].Mass;
          starID=P[i].ID;
        }
      }
#endif
    }
#if defined BETA_COOLING || defined ADD_CENTRAL_GRAVITY
  /* Gather all the sink IDs together */
  list_starMass = malloc(NTask * sizeof(double));
  list_starID = malloc(NTask * sizeof(int));
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allgather(&starMass,1,MPI_DOUBLE,list_starMass,1,MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgather(&starID,1,MPI_INT,list_starID,1,MPI_INT,MPI_COMM_WORLD);
  /* Now identify the final ID that we'll use*/
  for(i=0;i<NTask;i++)
  {
    if(list_starMass[i]>starMass)
    {
      starMass=list_starMass[i];
      starID=list_starID[i];
    }
  }
  /* Now we all have the same starMass and starID */
  All.StarID=starID;
  if(ThisTask==0)
    printf("The star ID is %d\n",All.StarID);
  free(list_starMass);
  free(list_starID);
#endif

#ifdef INJECT_GAS
  All.LastInjectionTime = All.TimeBegin;
  All.MaxID = All.TotNumPart;
  All.Injected = 0;
  //Determine Binary mass, a and j
  binaryMasses[0]=binaryMasses[1]=-1;
  for(i=N_gas;i<NumPart;i++)
  {
    binaryMasses[i-N_gas]=P[i].Mass;
    binaryPositions[(i-N_gas)*3+0] = P[i].Pos[0];
    binaryPositions[(i-N_gas)*3+1] = P[i].Pos[1];
    binaryPositions[(i-N_gas)*3+2] = P[i].Pos[2];
  }
  //Gather them all together
  list_starMass = malloc(NTask * 2 *sizeof(double));
  list_starPos = malloc(NTask * 6 *sizeof(double));
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allgather(&binaryMasses,2,MPI_DOUBLE,list_starMass,2,MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgather(&binaryPositions,6,MPI_DOUBLE,list_starPos,6,MPI_DOUBLE,MPI_COMM_WORLD);
  j=0;
  for(i=0;i<2*NTask;i++)
  {
    if(list_starMass[i]!=-1)
    {
      binaryMasses[j]=list_starMass[i];
      binaryPositions[j*3+0] = list_starPos[i*3+0];
      binaryPositions[j*3+1] = list_starPos[i*3+1];
      binaryPositions[j*3+2] = list_starPos[i*3+2];
      j++;
    }
  }
  //If we only found one star, make up binary properties
  if(j==1)
  {
    All.Binary_M = binaryMasses[0];
    All.Binary_a = 1.0;
    All.Binary_j = sqrt(All.G * All.Binary_M * All.Binary_a);
    All.Binary_q = 0.2;
  }
  else
  {
    //Calculate the mass,a and j and store them
    All.Binary_M = binaryMasses[0]+binaryMasses[1];
    All.Binary_a = sqrt((binaryPositions[0]-binaryPositions[3])*(binaryPositions[0]-binaryPositions[3]) + (binaryPositions[1]-binaryPositions[4])*(binaryPositions[1]-binaryPositions[4]) + (binaryPositions[2]-binaryPositions[5])*(binaryPositions[2]-binaryPositions[5]));
    All.Binary_j = sqrt(All.G * All.Binary_M * All.Binary_a);
    All.Binary_q = binaryMasses[0] > binaryMasses[1] ? binaryMasses[1]/binaryMasses[0] : binaryMasses[0]/binaryMasses[1];
  }
  if(ThisTask==0)
    printf("Binary has M=%g,a=%g,j=%g,q=%g\n",All.Binary_M,All.Binary_a,All.Binary_j,All.Binary_q);
  free(list_starMass);
  free(list_starPos);
#endif


#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

#ifdef FLEXSTEPS
  All.PresentMinStep = TIMEBASE;
  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
    }
#endif

#ifdef CDAV
  for(i=0; i < NumPart; i++)
    P[i].GravAccelOld[0]=P[i].GravAccelOld[1]=P[i].GravAccelOld[2]=0;
#endif

  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].HydroAccel[j] = 0;
	}

      SphP[i].DtEntropy = 0;
#ifdef BETA_COOILNG
#ifdef OUTPUTRADIATEDENERGY
      SphP[i].DtRadiatedEnergy = 0;
      SphP[i].RadiatedEnergy = 0;
#endif
#endif
#ifdef SINK_PARTICLES
      SphP[i].AccretionTarget = 0; 
#endif
#if defined MMAV
//     SphP[i].Alpha=All.VariableViscAlphaMin;
      SphP[i].DtAlpha=0;
#endif
#ifdef CDAV
//      SphP[i].Alpha=0;
      SphP[i].DtAlpha=0;
//      SphP[i].MaxSignalVel=0;
      SphP[i].DivVelOld=0;
//      SphP[i].DivVel=0;
//      SphP[i].GravAccelOld[0]=SphP[i].GravAccelOld[1]=SphP[i].GravAccelOld[2]=0;
#endif
      if(RestartFlag == 0)
	{
	  SphP[i].Hsml = 0;
	  SphP[i].Density = -1;
	}
    }

  ngb_treeallocate(MAX_NGB);

  force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  domain_Decomposition();	/* do initial domain decomposition (gives equal numbers of particles) */

  ngb_treebuild();		/* will build tree */

  setup_smoothinglengths();

  TreeReconstructFlag = 1;

  /* at this point, the entropy variable normally contains the 
   * internal energy, read in from the initial conditions file, unless the file
   * explicitly signals that the initial conditions contain the entropy directly. 
   * Once the density has been computed, we can convert thermal energy to entropy.
   */
#ifndef ISOTHERM_EQS
  if(header.flag_entropy_instead_u == 0)
    for(i = 0; i < N_gas; i++)
      SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].Density / a3, GAMMA_MINUS1);
#endif
}


/*! This routine computes the mass content of the box and compares it to the
 *  specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(fabs(omega - All.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;
#ifdef SURFACE
  double sur_tgt,sur_deviation;
#endif

  if(RestartFlag == 0)
    {

#ifdef SURFACE
      sur_tgt = 1.2089939655123523 * pow(All.DesNumNgb,2.0/3.0);
      sur_deviation = dmax(
          1.2089939655123523 *pow(All.DesNumNgb+All.MaxNumNgbDeviation,2.0/3.0)-sur_tgt,
          sur_tgt - 1.2089939655123523 *pow(All.DesNumNgb-All.MaxNumNgbDeviation,2.0/3.0));
#endif

      for(i = 0; i < N_gas; i++)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }
#ifndef TWODIMS
	  SphP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  SphP[i].Hsml =
	    pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
#ifdef SURFACE
  SphP[i].SurHsml =
    //pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
    pow(1.0 / (M_PI) * sur_tgt * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
	}
    }
  
#ifdef SURFACE
  sur_density();
#endif
  density();
}


/*! If the code is run in glass-making mode, this function populates the
 *  simulation box with a Poisson sample of particles.
 */
#if (MAKEGLASS > 1)
void seed_glass(void)
{
  int i, k, n_for_this_task;
  double Range[3], LowerBound[3];
  double drandom, partmass;
  long long IDstart;

  All.TotNumPart = MAKEGLASS;
  partmass = All.Omega0 * (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G))
    * (All.BoxSize * All.BoxSize * All.BoxSize) / All.TotNumPart;

  All.MaxPart = All.PartAllocFactor * (All.TotNumPart / NTask);	/* sets the maximum number of particles that may */

  allocate_memory();

  header.npartTotal[1] = All.TotNumPart;
  header.mass[1] = partmass;

  if(ThisTask == 0)
    {
      printf("\nGlass initialising\nPartMass= %g\n", partmass);
      printf("TotNumPart= %d%09d\n\n",
	     (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000));
    }

  /* set the number of particles assigned locally to this task */
  n_for_this_task = All.TotNumPart / NTask;

  if(ThisTask == NTask - 1)
    n_for_this_task = All.TotNumPart - (NTask - 1) * n_for_this_task;

  NumPart = 0;
  IDstart = 1 + (All.TotNumPart / NTask) * ThisTask;

  /* split the temporal domain into Ntask slabs in z-direction */

  Range[0] = Range[1] = All.BoxSize;
  Range[2] = All.BoxSize / NTask;
  LowerBound[0] = LowerBound[1] = 0;
  LowerBound[2] = ThisTask * Range[2];

  srand48(ThisTask);

  for(i = 0; i < n_for_this_task; i++)
    {
      for(k = 0; k < 3; k++)
	{
	  drandom = drand48();

	  P[i].Pos[k] = LowerBound[k] + Range[k] * drandom;
	  P[i].Vel[k] = 0;
	}

      P[i].Mass = partmass;
      P[i].Type = 1;
      P[i].ID = IDstart + i;

      NumPart++;
    }
}
#endif
