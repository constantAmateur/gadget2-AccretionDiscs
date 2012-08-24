import os
import re
#Parameters
base="/home/myoung/Output/Tests/Scaling"
timename="cpu.txt"
folders=["1CPU","2CPU","4CPU","8CPU","12CPU","24CPU","48CPU","96CPU"]
folders=["1CPUHalfParticles","2CPUHalfParticles","4CPUHalfParticles","8CPUHalfParticles","12CPUHalfParticles","24CPUHalfParticles","48CPUHalfParticles"]
folders=["1CPUTenthParticles","2CPUTenthParticles","4CPUTenthParticles","8CPUTenthParticles","12CPUTenthParticles","24CPUTenthParticles","48CPUTenthParticles","96CPUTenthParticles"]
folders=["1CPU_5e6","2CPU_5e6","4CPU_5e6","8CPU_5e6","12CPU_5e6","24CPU_5e6","48CPU_5e6","96CPU_5e6","192CPU_5e6"]
#Load information
r=re.compile("Step ([0-9]+), Time: (.+), CPUs: ([0-9]+)")
timeData=[]
gold=inf
goldData=[]
for folder in folders:
  os.chdir(base)
  os.chdir(folder)
  nn=np.genfromtxt(timename,comments="S")
  f=open(timename)
  tmp=[r.findall(x) for x in f.readlines() if x[0]=="S"]
  ncpus=np.array([float(x[0][2]) for x in tmp])
  time=np.array([float(x[0][1]) for x in tmp])
  nstep=np.array([float(x[0][0]) for x in tmp])
  if not np.all(ncpus==ncpus[0]):
    raise ValueError("Number of CPUs used not constant!")
  nn=np.column_stack((nstep,time,nn))
  timeData.append([nn,ncpus[0]])
  if ncpus[0]<gold:
    gold=ncpus[0]
    goldData=[nn,ncpus[0]]
  
#Relative times for each part of code
for x in timeData:
  #Decomposition of time at a given step...
  labels=["Grav","Hydro","Domain","Potential","Drift","Kick","IO","SPH","PH"]
  ind=1+np.array([2,3,4,5,6,7,8,13,18])
  #Normalize all times by total
  a=x[0][:,ind]
  b=x[0][:,2]
  bb=np.repeat(b,len(ind))
  bb.shape=a.shape
  normed=a/bb
  #Now plot the median fraction of the time for each part of code...
  plt.figure()
  pie(np.median(normed,axis=0),labels=labels)
  plt.title(str(int(x[1]))+"CPUs")

#Scaling for different bits of the code
#set the number after the plus to the number from page 37 of user guide...
what=1+2
for what in xrange(2,12):
  datx=[]
  daty=[]
  for x in timeData:
    mm=min(x[0].shape[0],goldData[0].shape[0])
    t=x[0][:mm,1]
    first=(x[0][:mm,what])
    gs=(goldData[0][:mm,what])
    speedup=np.median(gs/first)
    datx.append(x[1])
    daty.append(speedup)
    style='-'
    if what>8:
      style='--'
  plt.plot(np.log2(datx),np.log2(daty),style)
plt.legend(["Total","Grav","Hydro","Domain","Pot energy","Drift","Kick","Snapshot","Tree walks","Tree build"],loc='upper left')
plt.title("Scaling for 50k particles.")
plt.xlabel("log2(Number of processors)")
plt.ylabel("log2(Speedup)")
    #tmp=plt.plot(t[:-1],np.log2(diff(first)/diff(gs)))
    #tmp=plt.plot(x[0][:mm,1],((x[0][:mm,what]/goldData[0][:mm,what])/(goldData[1]/x[1]))**(-1),'--')
    #tmp=plt.plot(x[0][:mm,1],((x[0][:mm,what]/goldData[0][:mm,what])/(1.))**(-1),'--')
    #plt.plot(x[0][:mm,1],np.repeat(goldData[1]/x[1],mm),tmp[0].get_color()+'--')
