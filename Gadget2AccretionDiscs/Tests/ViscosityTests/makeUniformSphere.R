#Makes a Plummer Sphere
#############################
# PREAMBLE/PARAMATERS   #####
#############################
source("/home/myoung/Projects/Common/Rinit.R")

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

#radius
rs=100.0
#Number of sphere particles
Npart = 2e4
#Total mass
M = 1.0
#How does the temperature scale with radius? 0 gives isothermal
temppower = 0.0
#What is the temperature for the sphere..
baseTemp = 167.00
#We need this to convert a temperature to an internal energy....
gamma = 6/5.0
#Moleculare weight in units of hydrogen mass
mu=2.3

############I/O################

outputdir=paste(getwd(),"Build/",sep='/')

#############################
# BUILD INITIAL CONDITIONS  #
#############################

rho=M/((4/3)*pi*rs^3)
numbox=(6*Npart)/pi
#This creates a grid...
#x=seq(-rs,rs,((4*pi*rs^3)/(3*Npart))^(1/3))
#y=seq(-rs,rs,((4*pi*rs^3)/(3*Npart))^(1/3))
#z=seq(-rs,rs,((4*pi*rs^3)/(3*Npart))^(1/3))
#aa=rep(x,each=length(y)*length(z))
#bb=rep(y,each=length(z),times=length(x))
#cc=rep(z,times=length(x)*length(y))
#pos=data.frame(x=aa,y=bb,z=cc)
#pos=pos[aa^2+bb^2+cc^2<=rs^2,]
#This creates a random sampling
x=runif(numbox,min=-rs,max=rs)
y=runif(numbox,min=-rs,max=rs)
z=runif(numbox,min=-rs,max=rs)
pos=data.frame(x=x,y=y,z=z)
pos=pos[x^2+y^2+z^2<rs^2,]
Npart=dim(pos)[1]
vel=data.frame(vx=rep(0,Npart),vy=rep(0,Npart),vz=rep(0,Npart))
temp=rep(baseTemp,Npart)
mass=rep(M/Npart,Npart)
u=temp*k_b/(mu*m_H*(gamma-1))



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
