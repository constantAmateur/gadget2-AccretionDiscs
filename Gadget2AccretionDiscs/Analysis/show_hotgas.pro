


    base="/home/vspringe/Examples/lcdm_gas/"

    num=5

    Boxsize= 50000.0


    exts='000'
    exts=exts+strcompress(string(num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)
    f=base+"snapshot_"+exts
    f=strcompress(f,/remove_all)
    fname=f 

    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedback=0L
    npartTotal=lonarr(6)	
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
    la=intarr(bytesleft/2)


    print,fname
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
    NHalo=npart(1)
    NDisk=npart(2)
    NBulge=npart(3)
    NStars=npart(4)

    if Ngas gt 0 then begin
        u=fltarr(Ngas)
        readu,1,u

        rho=fltarr(Ngas)
        readu,1,rho
    endif
    close,1



    if Ngas gt 0 then begin
        xgas=fltarr(Ngas) &  ygas=fltarr(Ngas)  & zgas=fltarr(Ngas) & mgas=fltarr(Ngas)
        xgas(*)=pos(0,0:Ngas-1)
        ygas(*)=pos(1,0:Ngas-1)
        zgas(*)=pos(2,0:Ngas-1)
        if massarr(0) eq 0 then begin
            mgas(*)=mass(0:Ngas-1)	
        endif else begin
            mgas(*)= massarr(0)
	endelse
        print,"Mass gas: ",total(total(double(mgas)))
    endif



    GRAVITY   = 6.672d-8
    BOLTZMANN = 1.3806d-16
    UnitLength_in_cm =        3.085678d21        ;  1.0 kpc 
    UnitMass_in_g    =        1.989d43 ;  1.0e10 solar masses 
    UnitVelocity_in_cm_per_s = 1d5 ;  1 km/sec 
    
    UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s 
    G=GRAVITY/ UnitLength_in_cm^3 * UnitMass_in_g * UnitTime_in_s^2
    UnitDensity_in_cgs= UnitMass_in_g/ UnitLength_in_cm^3
    UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ UnitTime_in_s^2
    UnitEnergy_in_cgs= UnitMass_in_g * UnitLength_in_cm^2 / UnitTime_in_s^2


    PROTONMASS = 1.6726e-24

    ; convert to physical units 

    rho= rho* UnitDensity_in_cgs
    u  = u* UnitEnergy_in_cgs/ UnitMass_in_g

    MeanWeight= 1.0 * PROTONMASS


    gamma= 5.0/3

    temp = MeanWeight/BOLTZMANN * (gamma-1) * u


    ind=where(temp gt 5.0e6) ; select particles above 5.0*10^6 Kelvin 


    window, xsize=500, ysize=500

    plot, ygas(ind), zgas(ind), psym=3, yrange=[0, BoxSize], xrange=[0, BoxSize], xstyle=1, ystyle=1


end









