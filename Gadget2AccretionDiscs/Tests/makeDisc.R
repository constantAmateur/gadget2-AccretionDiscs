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

#######DISK VARIABLES#########

#Number of disk particles
Npart = 1e4
#Mass of the central star
M_star = 1.
#Mass of the disk (all disk particles given same mass)
M_disk = .1
#Inner radius (vacant area around star)
r_inner = 10
#Outer radius (edge of disk)
r_outer = 100
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

############I/O################

outputdir=paste(getwd(),"Build/",sep='/')

#############################
# BUILD INITIAL CONDITIONS  #
#############################

#We want the density to go like R^plIndex, so we have to pick the density of R values to go like R^plIndex+1
tmp=runif(Npart)
radi= ((r_outer^(2+plIndex)-r_inner^(2+plIndex))*tmp+r_inner^(2+plIndex))^(1/(2+plIndex))
#Theta probably should never be anything but uniform...
theta= runif(Npart,0,2*pi)
#Initialize the partcile temperatures
temp = baseTemp * ((radi/r_inner)^temppower)
#Calculate the speed of sound for each particle
cs = sqrt(gamma*k_b*temp/(mu*m_H))
#Calculate the velocities (just the usual 1/sqrt(r) keplerian)
vphi=sqrt(ginternal * M_star/radi)
#This determines the vertical scale of the disk, it goes as r^(3/2-temppower/2)
scale_H = (radi*cs /vphi)
#The z coordinate
z= rnorm(Npart)*scale_H
#An array of masses...
mass=rep(M_disk/Npart,Npart)
#We've calculated everything, now put it in the final variables which are cartesian...
pos=cyl2cart(cbind(radi,theta,z))
#The 1/r is to convert vphi to an angular velocity (i.e. dtheta/dt)
vel=cbind(-pos[,2]*vphi/radi,pos[,1]*vphi/radi,0)
#We need to convert the temperature to internal energy per unit mass
u=temp*cv*k_b/(mu*m_H)
#Add on the star...  It doesn't need an internal energy as it's a different "type"
pos=rbind(pos,c(0,0,0))
vel=rbind(vel,c(0,0,0))
mass=c(mass,M_star)
#Factor...
Hfact=sqrt((gamma*k_b*baseTemp*r_inner^-temppower)/(mu*m_H*ginternal*M_star))
Havg=Hfact*(2/((5+temppower)*(r_outer-r_inner)))*(r_outer^(.5*(5+temppower))-r_inner^(.5*(5+temppower)))

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
