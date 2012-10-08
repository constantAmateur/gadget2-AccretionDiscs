#############################
# DESCRIPTION ###############
#############################
#This script, together with "EasyIC" (https://bitbucket.org/constantAmateur/easyic) creates a GADGET format binary initial condition file of an accretion disc.  The parameters describing this accretion disc fall into arbitrary and scale free parameters.  The arbitrary parameters do not change the physical outcome or the run time of the simulation, whereas the scale free parameters change both.

#This disc is setup so that Q does not vary with R initially.  Here Q=c_s Omega/pi G Sigma

#############################
# PREAMBLE/PARAMATERS   #####
#############################
source("~/Projects/Common/Rinit.R")

### SCALE FREE PARAMETERS ###

#r=R_{outer}/R_{inner} the ratio of the outer to inner radii
r=5
#q=M_{disc}/M_* the mass ratio
q=.2
#density_power = Power law index for the surface density profile.  i.e. Sigma ~ R^density_power
density_power=-2.0
#Number of particles...
Npart = 1e6
#Fraction of the mass that should be in the exponential decay
epsilon = .01

### ARBITRARY PARAMETERS  ###

#Inner radius (in units given below)
R_i = 1
#Mass of star (in units given below)
M = 1.0
#Initial value of Q (Q is constructed to be constant initally)
Q_init = 2.0
#ratio of specific heats
gamma = 5/3.0
#Moleculare weight in units of hydrogen mass
mu=2.3

#########UNITS###############

#Unit Mass in grams
M_unit=1.989e33 #Solar mass
#Unit Length in cm
L_unit=1.496e13 #An AU
#We want the unit of time to be one outer rotation periods, so we set that and let the velocity unit follow...
#That is t_unit=2pisqrt(R^3/GM), giving t in units of seconds
#t_unit=2*pi*sqrt(((25*1.496e13)^3)/(6.673e-8*1.989e33))
#This unit choice gives G=1 (or at least it should).  It sets the time unit to be 1/Omega at R=1
t_unit=sqrt((L_unit^3)/(6.673e-08*M_unit))
#Infer the velocity unit
v_unit=L_unit/t_unit
##Unit Velocity in cm/sec
#v_unit=1e5
##Infer unit time
#t_unit=L_unit/v_unit
#Define some other units we'll use
#Boltzman constant in grams, cm, seconds, then converted to internal units
k_b=1.38e-16
k_b=k_b*(1/v_unit)^2*(1/M_unit)
#Mass of hydrogen in grams and internal units
m_H = 1.67e-24
m_H = m_H * (1/M_unit)
#Newton's Gravitational constant in grams,cm,seconds
ginternal=6.673e-08
#Convert to internal units...
ginternal = ginternal*(1/L_unit)^3*(1/M_unit)^-1*(1/t_unit)^-2
#These are the units that need to be set in Gadget...
print("Set the GADGET units to the following...")
print(paste("The unit velocity in cm/s is:",v_unit))
print(paste("The unit length in cm is:",L_unit))
print(paste("The unit mass in g is:",M_unit))
print(paste("The unit time s is:",t_unit))
print(paste("The value of G in this system is:",ginternal))

### Dependent variables #####
M_disc = M*q
R_o = R_i * r
#This ensures Q is constant with R
temp_power= 3 + 2*density_power
#First calculate the Sigma normalization constant
if(density_power==-2.0){
  Sigma_norm=M_disc/(2*pi*log(r))
}else{
  Sigma_norm=(M_disc * (density_power+2))/(2*pi*(R_o^(density_power+2)-R_i^(density_power+2)))
}
Temp_norm = (Q_init*Q_init*mu*m_H*ginternal*pi*pi*Sigma_norm*Sigma_norm)/(gamma*k_b*M)

############I/O################

outputdir=paste(getwd(),"Build/",sep='/')


#############################
# BUILD INITIAL CONDITIONS  #
#############################

#The error function
erf=function(x)
{
  return(2*pnorm(x*sqrt(2))-1)
}

#We need to pick the exponential decay in such a way that the surface density remains continuous, its 1st derivative remains continuous and the mass enclosed is equal to epsilon*M.  The root to the following function are what are needed to meet these requirements
outpeach = function(b,eps,Ri,r,alpha=density_power)
{
  return(Ri*(Ri-b)*(1/alpha)*(exp((alpha*(Ri-b))/(2*Ri))-exp(-(alpha*b^2)/(2*Ri*(b-Ri))))+b*sqrt((pi*Ri*(b-Ri))/(2*alpha))*(erf(sqrt(((b-Ri)*alpha*.5)/Ri))-erf(sqrt((alpha*b^2)/(2*Ri*(b-Ri)))))+(eps*Ri^(-alpha)*log(r)*exp((alpha*.5*(Ri-b))/Ri))/(1-eps))
}
#Determine where to look for root
samp=R_i*(0:1000)/1000
samp=samp[which(samp!=1)]
tmp=outpeach(samp,eps=epsilon,Ri=R_i,r=r,alpha=density_power)
lower=samp[max(which(tmp<0))]
print(paste("searching for root above R=",lower))
#The mean of the distribution
b=uniroot(outpeach,eps=epsilon,Ri=R_i,r=r,alpha=density_power,lower=lower,upper=R_i,f.upper=(R_i^(-density_power)*epsilon*log(r))/(1-epsilon))
b=b$root
#it's standard deviation
c=sqrt(R_i*(b-R_i)/density_power)
print(paste("Root finder settled on standard deviation =",c,"mean=",b))
#Not really needed, but for completeness...
#needs updating...
a=((1-epsilon)*M_disc*R_i^density_power*exp((b-R_i)^2/(2*c*c)))/(2*pi*log(r))
pdf=function(R)
{
  ret=(a*R*exp(-(R-b)^2/(2*c^2)))
  ret[R>=R_i]=(((1-epsilon)*M_disc*R^(density_power+1))/(2*pi*log(r)))[R>=R_i]
  return(ret)
}
rejectionSample = function(no,pdf,xmin=0,xmax=1,ymin=NULL,ymax=NULL,sample=100,safteyFact=2)
{
  if(is.null(ymin) || is.null(ymax))
  {
    samp=(xmax-xmin)*(0:sample)/sample+xmin
    tmp=pdf(samp)
    ymin=min(tmp)/safteyFact
    ymax=max(tmp)*safteyFact
  }
  res=vector()
  #Keep sampling till we have enough points
  while(length(res)<no)
  {
    #We'll make too many this way, but who cares...
    tmp_x=(xmax-xmin)*runif(no)+xmin
    tmp_y=(ymax-ymin)*runif(no)+ymin
    true_y=pdf(tmp_x)
    res=c(res,tmp_x[tmp_y<true_y])
  }
  #Only keep the first no
  res=res[1:no]
  return(res)
}
   



#Rejection sample the distribution...
radi = rejectionSample(Npart,pdf,xmax=R_o)

#First pick a random set of seeds which determines if the particle should go in the disc or in the exponential decay
#seed=runif(Npart)
#epart=which(seed<epsilon)
#print(paste("Assigning",length(epart),"particles to the exponential decay."))
#
##We want Sigma (the surface density) to go like R^density_power, so we have to pick the density of R values to go like R^density_power+1
#tmp=runif(Npart)
#if(density_power==-2.0){
#  radi = R_i * exp(tmp*log(r))
#}else{
#  radi = ((R_o^(2+density_power)-R_i^(2+density_power))*tmp+R_i^(2+density_power))^(1/(2+density_power))
#}
##Need to select radii for the exponential decay particles from the gaussian a*exp(-(R-b)^2/2c^2)
#tmp=rnorm(length(epart),mean=b,sd=c)
#rejectionSample(length(epart),pdf,xmax=R_i)
#tmp[tmp>R_i]=2*R_i-tmp[tmp>R_i]
##tmp[tmp>=R_i]=-1
##Shouldn't really be any of these I would hope...
#if(sum(tmp<=0))
#  print(paste("There were",sum(tmp<=0),"particles that were too close to the star and were accreted."))
#tmp[tmp<=0]=-1
#radi[epart]=tmp
##If there are any bad ones, accrete them
#M=M+(M_disc/Npart)*sum(radi<0)
#radi=radi[radi>0]
#Npart=length(radi)

#Theta probably should never be anything but uniform...
theta= runif(Npart,0,2*pi)
#The temperature
temp = Temp_norm*radi^temp_power
#Calculate the speed of sound for each particle
cs = sqrt(gamma*k_b*temp/(mu*m_H))
#Calculate the velocities (just the usual 1/sqrt(r) keplerian)
vphi=sqrt((ginternal * M)/(radi^3))
#Add the pressure modification? This saves time when settling into marginal stability
vphi=sqrt((ginternal * M)/(radi))
vphi=vphi*sqrt(1-density_power*cs*cs/(vphi*vphi))
vphi=vphi/radi
#This determines the vertical scale of the disk, it goes as r^(3/2-temppower/2)
scale_H = (cs /vphi)
#The z coordinate
z= rnorm(Npart)*scale_H
#An array of masses...
mass=rep(M_disc/Npart,Npart)
#We've calculated everything, now put it in the final variables which are cartesian...
pos=cyl2cart(cbind(radi,theta,z))
vel=cbind(-pos[,2]*vphi,pos[,1]*vphi,0)
#We need to convert the temperature to internal energy per unit mass
u=temp*k_b/(mu*m_H*(gamma-1))
#Add on the star...  It doesn't need an internal energy as it's a different "type"
pos=rbind(pos,c(0,0,0))
vel=rbind(vel,c(0,0,0))
mass=c(mass,M)

#############################
# WRITE TO FILE             #
#############################

#These are the filenames in "order"
noms=c(filefile="Filenames.txt",Header="Header.txt",Pos="Position.txt",Vel="Velocity.txt",ID="ID.txt",Mass="Masses.txt",InternalEnergy="Energy.txt")
tmp=names(noms)
noms=paste(outputdir,noms,sep='')
names(noms)=tmp
#Create the output directory (if exists, nothing is done)
dir.create(outputdir,recursive=TRUE,showWarnings=FALSE)
#Write the Master list
write.table(cbind(names(noms)[-1],as.character(noms)[-1]),file=noms[1],quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
#First write the header
header=list(particles=c(Npart,1,0,0,0,0),MassArray=rep(0,6),time=0,redshift=0,starFormFlag=0,feedbackFlag=0,totalParticles=c(Npart,1,0,0,0,0),coolingFlag=0,numFiles=1,boxSize=0,Omega0=0,OmegaLambda=0,Hubble=1,starFormTimeFlag=0,metalFlag=0,Npart64=rep(0,6),entropyNotEnergyFlag=0)
cat("#This is the header file!\n",file=noms["Header"])
for(i in 1:length(header)){
	cat(paste("#",names(header)[i],"\n",sep=''),file=noms["Header"],append=TRUE)
	#If it's a large integer, must not use scientific notation or things will break...
	tmp=options("scipen")[[1]]
	if(names(header)[i]%in%c("particles","totalParticles","Npart64"))
		options(scipen=100000)
	cat(header[[i]],file=noms["Header"],append=TRUE,sep='\t')
	options(scipen=tmp)
	cat("\n",file=noms["Header"],append=TRUE)
}
cat("#END",file=noms["Header"],append=TRUE)
#Now write out all the data and things
write.table(pos,sep='\t',file=noms["Pos"],quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(vel,sep='\t',file=noms["Vel"],quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(0:(Npart),file=noms["ID"],quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(mass,file=noms["Mass"],quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(u,file=noms["InternalEnergy"],quote=FALSE,row.names=FALSE,col.names=FALSE)
