;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Example of how file format of Gadget can be read-in  ;;
;; in IDL                                               ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    fname="/home/vspringe/ICs/galaxy.dat.intel"  ; Filename

    
    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedback=0L
    npartTotal=lonarr(6)	
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
    la=intarr(bytesleft/2)

    openr,1,fname,/f77_unformatted

    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,la

    N=total(npart)
    pos=fltarr(3,N)
    vel=fltarr(3,N)
    id=lonarr(N)
 
    ind=where((npart gt 0) and (massarr eq 0)) 
    if ind(0) ne -1 then begin
	Nwithmass= total(npart(ind))
        mass=fltarr(Nwithmass)
    endif else begin	
        Nwithmass= 0
    endelse

    readu,1,pos
    readu,1,vel
    readu,1,id

    if Nwithmass gt 0 then begin
      readu,1,mass
    endif

    NGas=npart(0)

    if Ngas gt 0 then begin

        u=fltarr(Ngas)
        readu,1,u

        rho=fltarr(Ngas)
        readu,1,rho
    endif
    close,1


end








