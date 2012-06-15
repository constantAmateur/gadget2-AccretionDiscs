#Makes a Plummer Sphere
#############################
# PREAMBLE/PARAMATERS   #####
#############################
source("~/Source/Common/Rinit.R")

#########UNITS###############

#Unit Mass in grams
M_unit=1.989e33 #Solar mass
#Unit Length in cm
L_unit=1.496e13 #An AU
#Unit Velocity in cm/sec
v_unit=1e5
#Infer unit time
t_unit=L_unit/v_unit
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

#######SPHERE VARIABLES#########

#Plummer radius
rs=1.0
#Number of sphere particles
Npart = 1e4
#Total mass
M = 1.0
#How does the temperature scale with radius? 0 gives isothermal
temppower = 0.0
#What is the temperature for the disc...
baseTemp = 16.700
#Adiabatic index
gamma = 5/3
#specific heat capacity at constant volume (for converting between temperature and internal energy)
cv=3/2
#Moleculare weight in units of hydrogen mass
mu=2.3
#Power law index, density (and column density) goes like R^plIndex
plIndex=-.5
#Should we distribute velocities according to Aarseth, Henon & Wielen?
distVel=TRUE

############I/O################

outputdir=paste(getwd(),"Build/",sep='/')

#############################
# BUILD INITIAL CONDITIONS  #
#############################

x1=runif(Npart)
x2=runif(Npart)
x3=runif(Npart)
radi = sqrt((x1^(2/3) * rs^2)/(ginternal^(2/3)-x1^(2/3)))
phi = pi*(2*x2-1)
theta=acos(2*x3-1)
pos=sph2cart(cbind(radi,phi,theta))
mass=rep(M/Npart,Npart)
#This is really just some nonesense that should be ignored
temp = baseTemp * ((radi)^temppower)
u=temp*cv*k_b/(mu*m_H)
#Routine for sampling velocities
if(distVel){
  q=rep(-1,Npart)
  while(any(q<0)){
    needed=which(q<0)
    x4=runif(length(needed))
    x5=runif(length(needed))
    good=which(.1*x5 < (1-x4^2)*x4^2)
    q[needed[good]]=x4[good]
  }
  vesc = sqrt((2*ginternal*M)/((radi^2 +rs^2)^.25))
  v = vesc*q
  x6 = runif(Npart)
  x7 = runif(Npart)
  vx=(1-2*x6)*v
  vy=sqrt(v^2-vx^2)*cos(2*pi*x7)
  vz=sqrt(v^2-vx^2)*sin(2*pi*x7)
  vel=cbind(vx,vy,vz)
}else{
  tmp=rep(0,Npart)
  vel=cbind(tmp,tmp,tmp)
}


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
header=list(particles=c(Npart,0,0,0,0,0),MassArray=rep(0,6),time=0,redshift=0,starFormFlag=0,feedbackFlag=0,totalParticles=c(Npart,0,0,0,0,0),coolingFlag=0,numFiles=1,boxSize=0,Omega0=0,OmegaLambda=0,Hubble=1,starFormTimeFlag=0,metalFlag=0,Npart64=rep(0,6),entropyNotEnergyFlag=0)
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
write.table(0:(Npart-1),file=noms["ID"],quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(mass,file=noms["Mass"],quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(u,file=noms["InternalEnergy"],quote=FALSE,row.names=FALSE,col.names=FALSE)
