#Makes a Evrard Sphere
#This is designed to reproduce the "Evrard test" from the Cullend & Dullen (2010) artificial viscosity paper.  This in turn is based on the 1993 Steinmetz & Muller paper which describes a series of SPH tests and gives "gold standard" grid based answers to compare against.  Note that gamma should be 5/3.
#############################
# PREAMBLE/PARAMATERS   #####
#############################
source("/home/myoung/Projects/Common/Rinit.R")

#########UNITS###############

#Unit Mass in grams
M_unit=1.989e33 #Solar mass
#Unit Length in cm
L_unit=1.496e13 #An AU
#Unit Velocity in cm/sec # Needed for G=1
v_unit=2.979e6
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

#radius should be 1...
rs=1.0
#Number of sphere particles
Npart = 100280
#Total mass.  Should be 1...
M = 1.0

############I/O################

outputdir=paste(getwd(),"Build/",sep='/')

#############################
# BUILD INITIAL CONDITIONS  #
#############################

#First we have to setup a Face centered cubic latice
rho=M/((4/3)*pi*rs^3)
numbox=(6*Npart)/pi
#Make one cell
base=matrix(c(0,0,0,0,0,1,0,1,0,0,1,1,.5,.5,0,0,.5,.5,.5,0,.5),c(7,3),byrow=T)
a=which(base==0)
b=which(base==1)
copy=base
copy[a]=1
copy[b]=0
base=rbind(base,copy)
#For each box, define its top-left corner
x=0:ceiling((4*Npart)^(1/3))
y=0:ceiling((4*Npart)^(1/3))
z=0:ceiling((4*Npart)^(1/3))
xx=rep(x,each=length(y)*length(z))
yy=rep(y,each=length(z),times=length(x))
zz=rep(z,times=length(x)*length(y))
#14 is the number of points per box
allBoxes=as.matrix(data.frame(x=rep(xx,each=14),y=rep(yy,each=14),z=rep(zz,each=14)))
stop=dim(allBoxes)[1]/14
a=matrix(rep(c(t(base)),stop),ncol=3,byrow=TRUE)
final=a+allBoxes
lat=final[!duplicated(final),]
#Scale up so all points are >1 or 0
lat=lat*10
#Get the point closest to the middle
mid=c(max(x)/2,max(y)/2,max(z)/2)
tmp=lat-mid
dist=(tmp[,1]^2+tmp[,2]^2+tmp[,3]^2)
mid=lat[which(dist==min(dist))[1],]
#Shift everything relative to this
lat=lat-mid
#Now we can finally stretch the fucker out
r=sqrt(lat[,1]^2+lat[,2]^2+lat[,3]^2)^(1/3)
lat=lat*matrix(c(r,r,r),ncol=3)
#Keep only as many particles as we need
r=sqrt(lat[,1]^2+lat[,2]^2+lat[,3]^2)
maxr=r[order(r)[Npart]]
pos=lat[which(r<=maxr),]
#Now we need to rescale everything so that rmax=1
lat=lat*(rs/maxr)
#All done...
Npart=dim(pos)[1]
vel=data.frame(vx=rep(0,Npart),vy=rep(0,Npart),vz=rep(0,Npart))
u=rep(.05*(ginternal*M/rs),Npart)
mass=rep(M/Npart,Npart)


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

#######################
#Useful diagnostics   #
#######################

print(paste("Unit velocity in cm/s:",v_unit))
print(paste("Unit length in cm:",L_unit))
print(paste("Unit mass in g:",M_unit))
print(paste("Unit time in s:",t_unit))
print(paste("G in this unit system:",ginternal))
print(paste("Number of particles used:",Npart))
