import numpy as np
fixed=np.genfromtxt("fixedGravSmooth/energy.txt")
price=np.genfromtxt("priceGravSmooth/energy.txt")
var=np.genfromtxt("varGravSmooth/energy.txt")

plt.plot(price[:,0],price[:,1]+price[:,2]+price[:,3])
plt.plot(var[:,0],var[:,1]+var[:,2]+var[:,3])
plt.plot(fixed[:,0],fixed[:,1]+fixed[:,2]+fixed[:,3])
plt.xlabel("Time")
plt.ylabel("Total Energy")
plt.title("Total energy for N=1.5 perturbed polytrope")
plt.legend(("Adaptive gravitational softening with correction terms","Adaptive gravitational softening without correction terms"))


