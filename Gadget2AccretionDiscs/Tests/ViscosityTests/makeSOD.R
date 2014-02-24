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

#Y-Z size
l=100.0
#Main dimension
x=100.0
#Number of particles
Npart = 1e5
#Adiabatic index
gamma = 1.4
#high density/pressure
rho1=1.0
p1=1.0
#low density/pressure
rho2=.25
p2=.1795
#Total mass
M=(rho1*l*l*x*.5)+(rho2*l*l*x*.5)
#Moleculare weight in units of hydrogen mass
mu=2.3

############I/O################

outputdir=paste(getwd(),"Build/",sep='/')

#############################
# BUILD INITIAL CONDITIONS  #
#############################

#Tile the high half
Lhigh=(M/(rho1*Npart))^(1/3)
a=seq(0,l,Lhigh)
b=seq(0,l,Lhigh)
c=seq(-(x*.5),0,Lhigh)
aa=rep(a,each=length(b)*length(c))
bb=rep(b,each=length(c),times=length(a))
cc=rep(c,times=length(a)*length(b))
highHalf=data.frame(x=cc,y=bb,z=aa)
#Tile the low half
Llow=(M/(rho2*Npart))^(1/3)
a=seq(0,l,Llow)
b=seq(0,l,Llow)
c=seq(max(c)+Llow,(x*.5),Llow)
aa=rep(a,each=length(b)*length(c))
bb=rep(b,each=length(c),times=length(a))
cc=rep(c,times=length(a)*length(b))
lowHalf=data.frame(x=cc,y=bb,z=aa)
#Combine them
pos=rbind(highHalf,lowHalf)
Npart=dim(pos)[1]
#Now calculate the internal energies
u1=(rho1/p1)*(1/(gamma-1))
u2=(rho2/p2)*(1/(gamma-1))
u=c(rep(u1,dim(highHalf)[1]),rep(u2,dim(lowHalf)[1]))
#Other things
vel=pos
vel[,1]=vel[,2]=vel[,3]=0
mass=rep(M/Npart,dim(pos)[1])

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
