
limit=301
import numpy as np
import os
from matplotlib import pyplot as plt

#Calculate potential from ascii
def potential(base,stop=301):
  start=0
  energy=np.zeros((stop-start+1,2))
  for i in xrange(start,stop+1):
    nom=base % i
    tmp=np.genfromtxt(nom)
    x=tmp[:,0]
    y=tmp[:,1]
    z=tmp[:,2]
    pe=0
    for j in xrange(0,len(x)):
      pe = pe - np.sum(1.0/np.sqrt((x[j]-x[j+1:])**2 + (y[j]-y[j+1:])**2 + (z[j]-z[j+1:])**2))
    energy[i,0]=i
    energy[i,1]=pe
  return energy[:,1]


#Gravity test
#Original
first="OrigAV_FixedGrav/energy.txt"
second="Pristene_OrigAV_VarGrav/energy.txt"
third="OrigAV_PriceGrav/energy.txt"
#New
first="NewAV_FixedGrav/energy.txt"
second="NewAV_VarGrav/energy.txt"
third="NewAV_PriceGrav/energy.txt"
#AV test
first="OrigAV_FixedGrav/energy.txt"
second="NewAV_FixedGrav/energy.txt"
third="NewAV_FixedGrav/energy.txt"
#potential returns sum(1/r), multiply it by this factor...
potFac=5e-7
calcPot=False

fdat=np.genfromtxt(first)
fdat=fdat[:limit,]
sdat=np.genfromtxt(second)
sdat=sdat[:limit,]
tdat=np.genfromtxt(third)
tdat=tdat[:limit,]
#Calculate the potential?
if calcPot:
  fpot=potential(os.path.dirname(first)+"/gassphere_%03d.ascii")*potFac
  fdat[:,2]=fpot[:limit]
  spot=potential(os.path.dirname(second)+"/gassphere_%03d.ascii")*potFac
  sdat[:,2]=spot[:limit]
  tpot=potential(os.path.dirname(third)+"/gassphere_%03d.ascii")*potFac
  tdat[:,2]=tpot[:limit]


#First column is time, second is thermal, third is potential, fourth is kinetic
t=(fdat[:,0]+sdat[:,0]+tdat[:,0])/3.0
#Thermal
plt.plot(t,fdat[:,1],'r-')
plt.plot(t,sdat[:,1],'r--')
plt.plot(t,tdat[:,1],'r-.')
#Kinetic
plt.plot(t,fdat[:,3],'g-')
plt.plot(t,sdat[:,3],'g--')
plt.plot(t,tdat[:,3],'g-.')
#potential
plt.plot(t,fdat[:,2],'b-')
plt.plot(t,sdat[:,2],'b--')
plt.plot(t,tdat[:,2],'b-.')
#Total
plt.plot(t,fdat[:,1]+fdat[:,2]+fdat[:,3],'-',color='black')
plt.plot(t,sdat[:,1]+sdat[:,2]+sdat[:,3],'--',color='black')
plt.plot(t,tdat[:,1]+tdat[:,2]+tdat[:,3],'-.',color='black')
#Pretty things up
plt.xlabel("Free Fall Times")
plt.ylabel("Energy")
plt.title("Energy evolution for different code setups")
#plt.legend(('$E_{th}^{def}$','$E_{th}^{vis}$','$E_{th}^{grav+vis}$',
#'$E_{kin}^{def}$','$E_{kin}^{vis}$','$E_{kin}^{grav+vis}$',
#'$E_{pot}^{def}$','$E_{pot}^{vis}$','$E_{pot}^{grav+vis}$',
#'$E_{tot}^{def}$','$E_{tot}^{vis}$','$E_{tot}^{grav+vis}$'),ncol=4)

