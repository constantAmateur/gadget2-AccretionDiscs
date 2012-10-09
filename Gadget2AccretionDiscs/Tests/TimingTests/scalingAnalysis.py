import os
import re
#Parameters
base="/home/myoung/Output/Tests/Scaling/r5noPH"
timename="cpu.txt"
#Load information for all runs...
dat=[]
##This is the old runs with r=100
#folders=["1CPUTenthParticles","2CPUTenthParticles","4CPUTenthParticles","8CPUTenthParticles","12CPUTenthParticles","24CPUTenthParticles","48CPUTenthParticles","96CPUTenthParticles"]
#dat.append([50000,folders])
##Too finely grained to be useful
##folders=["1CPUHalfParticles","2CPUHalfParticles","4CPUHalfParticles","8CPUHalfParticles","12CPUHalfParticles","24CPUHalfParticles","48CPUHalfParticles"]
##dat.append([250000,folders])
#folders=["1CPU","2CPU","4CPU","8CPU","12CPU","24CPU","48CPU","96CPU"]
#dat.append([500000,folders])
#folders=["1CPU_5e6","2CPU_5e6","4CPU_5e6","8CPU_5e6","12CPU_5e6","24CPU_5e6","48CPU_5e6","96CPU_5e6","192CPU_5e6"]
#dat.append([5000000,folders])
##Lack of baseline confounds this analysis
#folders=["60CPU_5e7","120CPU_5e7","240CPU_5e7","480CPU_5e7"]
#dat.append([50000000,folders])
#These are the new runs with r=5...
folders=["1CPU_1e5","2CPU_1e5","4CPU_1e5","12CPU_1e5","24CPU_1e5","48CPU_1e5"]#,"96CPU_1e5"]#,"192CPU_1e5"]
dat.append([100000,folders])
folders=["1CPU_5e5","2CPU_5e5","4CPU_5e5","12CPU_5e5","24CPU_5e5","48CPU_5e5","96CPU_5e5"]#,"192CPU_5e5"]
dat.append([500000,folders])
folders=["1CPU_1e6","2CPU_1e6","4CPU_1e6","12CPU_1e6","24CPU_1e6","48CPU_1e6","96CPU_1e6","192CPU_1e6"]
dat.append([1000000,folders])
folders=["1CPU_5e6","2CPU_5e6","4CPU_5e6","12CPU_5e6","24CPU_5e6","48CPU_5e6","96CPU_5e6","192CPU_5e6"]
dat.append([5000000,folders])
folders=["1CPU_1e7","2CPU_1e7","4CPU_1e7","12CPU_1e7","24CPU_1e7","48CPU_1e7","96CPU_1e7","192CPU_1e7","384CPU_1e7"]
dat.append([10000000,folders])
#Load all information at once...
for i in xrange(len(dat)):
  folders=dat[i][1]
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
  dat[i].append(goldData)
  dat[i].append(timeData)

#For multiple N
plotDat=[]
for i in xrange(len(dat)):
  what=1+1
  datx=[]
  daty=[]
  speedup=1
  for c,x in enumerate(dat[i][3]):
    cind=max(0,c-1)
    mm=min(x[0].shape[0],dat[i][3][cind][0].shape[0])
    t=x[0][:mm,1]
    first=(x[0][:mm,what])
    gs=(dat[i][3][cind][0][:mm,what])
    speedup=speedup*np.min(gs/first)
    datx.append(x[1])
    daty.append(speedup)
  plotDat.append([datx,daty,dat[i][0]])
#Make some plots...
speedup=plt.figure()
speedup=speedup.add_subplot(1,1,1)
speedup.set_title("Total speedup for different numbers of particles")
speedup.set_xlabel("log2(Number of processors)")
speedup.set_ylabel("log2(Speedup)")
efficiency=plt.figure()
efficiency=efficiency.add_subplot(1,1,1)
efficiency.set_title("Scaling efficiency for different numbers of particles")
efficiency.set_xlabel("log2(Number of processors)")
efficiency.set_ylabel("speedup/No. proc")
speedup.bigx=[]
efficiency.bigx=[]
for tmp in plotDat:
  speedup.plot(np.log2(tmp[0]),np.log2(tmp[1]),label="No. pcles="+str(tmp[2]))
  if len(tmp[0])>len(speedup.bigx):
    speedup.bigx=tmp[0]
  efficiency.plot(np.log2(tmp[0]),np.array(tmp[1])/(np.array(tmp[0])/tmp[0][0]),label="No. pcles="+str(tmp[2]))
  if len(tmp[0])>len(efficiency.bigx):
    efficiency.bigx=tmp[0]

#Add guiding lines
speedup.plot(np.log2(speedup.bigx),np.log2(speedup.bigx),'k--',label="Perfect scaling")
efficiency.plot(np.log2(efficiency.bigx),np.repeat(1,len(efficiency.bigx)),'k--',label="Perfect scaling")
#Add legends
speedup.legend(loc=2)
efficiency.legend(loc=3)
#Set range
speedup.set_ylim(bottom=0)
efficiency.set_ylim(0,1.1)
#Display
speedup.get_figure().show()
efficiency.get_figure().show()





#For single N
i=0
goldData=dat[i][2]
timeData=dat[i][3]
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
  speedup=1
  #Assumes that timeData is ordered in increasing No. CPUs
  for c,x in enumerate(timeData):
    cind=max(0,c-1)
    mm=min(x[0].shape[0],timeData[cind][0].shape[0])
    t=x[0][:mm,1]
    first=(x[0][:mm,what])
    gs=(timeData[cind][0][:mm,what])
    speedup=speedup*np.median(gs/first)
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
