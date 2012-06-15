
limit=301
import numpy as np
from matplotlib import pyplot as plt

#first="GadgetOriginal/energy.txt"
first="NewViscPriceGrav/energy.txt"
second="NewVisc/energy.txt"
third="energy.txt"

fdat=np.genfromtxt(first)
fdat=fdat[:limit,]
sdat=np.genfromtxt(second)
sdat=sdat[:limit,]
tdat=np.genfromtxt(third)
tdat=tdat[:limit,]
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

