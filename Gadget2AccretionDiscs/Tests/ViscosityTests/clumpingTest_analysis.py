from scipy.stats import pearsonr
from scipy import spatial as sp
import numpy as np
import os
import re
import progressbar

dir="/home/myoung/Output/Tests/ViscosityTests/clumping/noVisc_reverseSpin"
ff=re.compile("clumping_([0-9]+).ascii")
nsamp=10

os.chdir(dir)
files=os.listdir(".")
files=[x for x in files if ff.match(x) is not None]
tpoints = np.array([int(ff.match(x).groups()[0]) for x in files])
widgets = [progressbar.Percentage(), progressbar.Bar(),' ',progressbar.FormatLabel('Time elapsed: %(elapsed)s'),' ',progressbar.ETA()]
samp=np.floor(np.linspace(0,len(files)-1,nsamp))
processed=np.zeros((nsamp,5))
output=[]
j=0
for t in samp:
  file=files[np.where(t==tpoints)[0]]
  dat=np.genfromtxt(file)
  pos=dat[:,(0,1,2)]
  r=np.sqrt(pos[:,0]**2+pos[:,1]**2)
  theta=np.arctan2(pos[:,1],pos[:,0])
  h=dat[:,9]
  #make a tree
  tree=sp.KDTree(pos)
  #Do each of them
  ccs=np.zeros(len(h))
  pbar2 = progressbar.ProgressBar(widgets=widgets, maxval=len(h)).start()
  for i in xrange(len(h)):
    pbar2.update(i)
    pts=tree.query_ball_point(pos[i,],h[i]*1.0)
    if len(pts)>3:
      #Damned modulo arithmatic
      tt=theta[pts]-theta[i]
      o=np.where(np.abs(tt)>np.pi)[0]
      if len(o)>0:
        tt[o]=tt[o]-sign(tt[o])*2*pi
      cc=pearsonr(r[pts]-r[i],tt)[0]
      if cc<-.5:
        ccs[i]=-1
      elif cc>.5:
        ccs[i]=1
      else:
        ccs[i]=2
  pbar2.finish()
  #Store and summarise
  processed[j,0]=t
  processed[j,1]=np.sum(ccs==1)
  processed[j,2]=np.sum(ccs==-1)
  processed[j,3]=np.sum(ccs==0)
  processed[j,4]=len(ccs)
  j=j+1
  output.append(ccs)

#Plot the progression
plt.plot(processed[:,0],processed[:,1]/processed[:,4])
plt.plot(processed[:,0],processed[:,2]/processed[:,4])

#Plot the neighbours of i
for i in prograde:
  pts=tree.query_ball_point(pos[i,],h[i])
  tt=theta[pts]-theta[i]
  o=np.where(np.abs(tt)>np.pi)[0]
  if len(o)>0:
   tt[o]=tt[o]-sign(tt[o])*2*pi
  plt.plot(tt,r[pts]-r[i],'.')
