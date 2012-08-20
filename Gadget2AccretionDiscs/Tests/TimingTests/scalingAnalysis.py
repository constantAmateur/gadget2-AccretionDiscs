import os
import re
#Parameters
base="/home/myoung/Output/Tests/TimingTests"
timename="cpu.txt"
folders=["1CPU","2CPU","4CPU"]
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
  
for x in timeData:
  tmp=plt.plot(x[0][:,1],x[0][:,2]/goldData[0][:len(x[0][:,2]),2])
  plt.plot(x[0][:,1],np.repeat(goldData[1]/x[1],len(x[0][:,1])),tmp[0].get_color()+'--')
