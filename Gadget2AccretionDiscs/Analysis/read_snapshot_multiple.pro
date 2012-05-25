
base="/afs/mpa/data/sun-3/volker/iap8/"


window,xsize=900,ysize=900
!x.margin=[6,6]
!y.margin=[6,6]


FILES=2

num=0

Etot=0
Ntot=0


for fi=0,FILES-1 do begin	

    exts='000'
    exts=exts+strcompress(string(num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)

    f=base+"snapshot_"+exts+"."+string(fi)
    f=strcompress(f,/remove_all)
    

    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=50.0D
    bytesleft=256-6*4 - 6*8 - 8 - 8
    
    la=intarr(bytesleft/2)


    

    openr,1,f,/f77_unformatted

    readu,1,npart,massarr,time,redshift,la

    print,npart,massarr
    print
    print,"Time=",time
    print,"Redshift=",1/time -1, redshift
    print


    N=total(npart)

    pos=fltarr(3,N)
    vel=fltarr(3,N)
    id=lonarr(N)

    readu,1,pos
    readu,1,vel
    readu,1,id
    close,1



    NN=npart(1)
    NN=N

    x=fltarr(NN) 
    y=fltarr(NN)
    z=fltarr(NN)
    x(*)=pos(0,0:NN-1)
    y(*)=pos(1,0:NN-1)
    z(*)=pos(2,0:NN-1)

;    vx=fltarr(NN) 
;    vy=fltarr(NN)
;    vz=fltarr(NN)
;    vx(*)=vel(0,0:NN-1)
;    vy(*)=vel(1,0:NN-1)
;    vz(*)=vel(2,0:NN-1)

;    E1=total(vx^2+vy^2+vz^2)
;    print,"E=",E1
;    Etot=Etot+E1
    Ntot=Ntot+NN
    
    print,"total: id=", total(double(id)),max(id)
    print,"total: N=", total(npart)
    

    ind = where(abs(z) le 5000) 
    if fi eq 0 then begin
;        plot,x(ind),y(ind),psym=3 ,xrange=[-15000,15000],yrange=[-15000,15000]

        plot,x(ind),y(ind),psym=3 ,xrange=1.0*[-15000,15000],yrange=1.0*[-15000,15000],xstyle=1,ystyle=1
    endif else begin
        oplot,x(ind),y(ind),psym=3
    endelse	
    
endfor


print,"<V^2>=",Etot/(Ntot)


end








