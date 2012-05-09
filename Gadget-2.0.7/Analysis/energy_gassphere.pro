

frun="/home/vspringe/Examples/gassphere/"


f=frun+"energy.txt"


spawn,"wc "+f,result
lines=long(result)
lines=lines(0)


en=fltarr(4+6*3+6,LINES)


openr,1,f
readf,1,en
close,1


ti=fltarr(LINES)
ke=fltarr(LINES)
po=fltarr(LINES)
th=fltarr(LINES)



ti(*)=en(0,*)
ke(*)=en(3,*)
po(*)=en(2,*)
th(*)=en(1,*)


tot=th+ke+po

window,xsize=500,ysize=800

!p.multi=[0,1,2]

plot,ti,tot,yrange=[-3,2],xrange=[0,3.0]

oplot,ti,th,linestyle=4,thick=3.0

oplot,ti,ke,linestyle=2

oplot,ti,po,linestyle=3



plot,ti,tot/tot(0) -1 


end




