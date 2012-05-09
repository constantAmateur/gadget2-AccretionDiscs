

frun="/home/vspringe/Examples/galaxy/"

f=frun+"/energy.txt"


spawn,"wc "+f,result
lines=long(result)
lines=lines(0)


en=dblarr(4+6*3+6,LINES)

openr,1,f
readf,1,en
close,1


ti=dblarr(LINES)
ke=dblarr(LINES)
po=dblarr(LINES)
th=dblarr(LINES)


ti(*)=en(0,*)
ke(*)=en(3,*)
po(*)=en(2,*)
th(*)=en(1,*)


tot=th+ke+po


window,xsize=500,ysize=700

!p.multi=[0,1,2]

plot,ti,tot, yrange=[min(po),max(ke)]	

oplot,ti,ke,linestyle=2
oplot,ti,po,linestyle=3


plot,ti,(tot-tot(0))/abs(tot(0)) 


end




