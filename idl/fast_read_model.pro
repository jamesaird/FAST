;........................................................................
;Description:
; Reads stellar population binary files (ised) and returns synthetic
; spectra for given ages
;
;Parameters
; input:     ised        ised file (bc03, cb07,ma05)
;            log_age     Ages (log) for which spectra are returned [yr]
; output:    wave        wavelength [Angstrom]
;            sfr         star formation rate M_sun / yr 
;            mass
; return:    fluxes      spectrum at each age [L_sun / Angstrom]
;
; optional keywords: 
;  id   = description of SED
;  imf  = [ ml, mu, nseg, [xx, lm, um, baux, cn, cc] ]
;  info = [totm, totn, avs, jo, tau, tau1, tau2, tau3, tau4, iop, stelib]
;
; Notes
; + Based on script Ivo Labbe
; + if log_age is supplied, the code interpolates between the nearest 
;   age steps in the SED, otherwise it returns the full age grid
;........................................................................


FUNCTION fast_read_model, ised, log_age, wave, sfr, mass, imf=imf, $
                          id=id, info=info, extras=extras, specs=specs

;... NOTE: we are reading a mixed type raw file so we need to tell idl
;    beforehand whether to read 4 or 8 byte words
tsteps = 0L ; 8 byte long integers
wsteps = 0L ; 8 byte long integers
junk   = 0e ; 4 byte
junk1  = bytarr(4)
junk2  = bytarr(12)

;... test endian by reading second word (time steps)
;    if it doesnt make sense swap endian
openr, lun, ised, /get_lun
readu, lun, junk, tsteps & close, lun 

;... open again with correct endian
;    NOTE: the swap_endian test is not super robust, could go wrong.
openr, lun,  ised, swap_endian=(tsteps lt 0 or tsteps gt 1e4)
readu, lun, junk, tsteps ; time steps
time = fltarr(tsteps)    ; read array of tsteps x 4 byte words
readu, lun, time  

;... the info section. what a freaking pain... 
;    NOTE: the order of the info block is *different* in the ASCII file 
;    and ised file
;    NOTE: in cb07 csp_galaxev code nextra has changed to 12 always.
imf_lim = fltarr(2)      ; lower, upper limit imf
imf_seg = 0l             ; number of imf segments: long integer
info1   = fltarr(5)    
info2   = fltarr(4)
id      = bytarr(80)
iop     = 1s             ; integer 2b
stelib  = 1s             ; "logical" 2b
readu, lun, imf_lim, imf_seg
if imf_seg gt 0 then imf_par = fltarr(imf_seg*6) else imf_par = 0
for i=0,imf_seg-1 do begin 
    tmp = fltarr(6) 
    readu, lun, tmp
    imf_par[i*6:(i+1)*6-1] = tmp
endfor 
readu, lun, info1, id, info2, id, id, iop, stelib, junk, junk, junk 
imf  = [imf_lim, imf_seg, imf_par]
info = [info1,info2,iop,stelib] 
id   = string(id)
if info1[3] eq 0 then nextra=12 else nextra=10 ; 2 more extra columns if SSP 

;... read seds at all time steps + extra spectral indices info
readu, lun, wsteps ; wave steps
wave   = fltarr(wsteps) 
readu, lun, wave 
fluxes = fltarr(tsteps,wsteps)
ssteps = 0L
for i=0,tsteps-1 do begin
    readu, lun, junk, junk, wsteps 
    flux = fltarr(wsteps)
    readu, lun, flux, ssteps ; get number of spectral index points  
    spec = fltarr(ssteps)
    readu, lun, spec         ; get spectral index points
    if i eq 0 then specs = fltarr(tsteps,ssteps) 
    fluxes[i,0:wsteps-1] = flux
    specs[i,0:ssteps-1]  = spec
endfor  
    
;... read extra sfr, mass, bolflux etc vs time at end of ised
;    bolflux, M*, sfr are the first three columns
extras = fltarr(tsteps,nextra)
for i=0,nextra-1 do begin
    readu, lun, junk, junk, tsteps 
    extra = fltarr(tsteps)
    readu, lun, extra ;   
    extras[0:tsteps-1,i] =  extra
endfor
close, lun & free_lun, lun

;... if log_age is set: interpolate fluxes etc onto supplied age grid
;    NOTE: extrapolations will select nearest neighbour value
;    NOTE: BC03: M* = Mgal - Mgas, MA05 = M*T = M*_alive+Mremn
log_age_model = alog10(time) > 0  
if not keyword_set(log_age) then begin
    log_age = alog10(time) > 0
    mass    = extras[*,1]
    sfr     = extras[*,2]    
    return, fluxes > 0.
endif else begin 
    age_ip = INTERPOL(INDGEN(tsteps), log_age_model, log_age)  
    extras = INTERPOLATE(extras, age_ip, INDGEN(nextra), /grid) > 0.
    specs  = INTERPOLATE(specs, age_ip, INDGEN(ssteps), /grid) > 0.
    mass   = extras[*,1] 
    sfr    = extras[*,2]    
    return, INTERPOLATE(fluxes, age_ip, INDGEN(n_elements(wave)), /grid) > 0. 
endelse

end

; the extras block contains 12 columns (10 if not an SSP):
; 	  bflx     bol flux
;	  strm     stellar mass (M* = Mtot-Mgas) 
;	  sfr      SFR
;	  evfl     spec. evol. flux
;	  snbr     SN rate/yr/Lsun
;	  pnbr     PN rate
;	  bhtn     # BH
;         sntn     # NS
;	  wdtn     # WD
;	  rmtm     mass in stellar remnants
;	  toff     turn off mass
;	  bolms    ratio of MS / PMS bol flux
; 
;  bin_ised/csp_galaxev decides if ised is an SSP by looking at [jo]
;  in info block (=info[3]) from bin_ised.f 
;     open (1,file=name,form='unformatted',status='old',err=3)
;     read  (1) nsteps,(tb(i),i=1,nsteps),ml,mu,iseg,(xx(i),lm(i),
;     um(i), baux(i),cn(i),cc(i),i=1,iseg),totm,totn,avs,jo,tau,id,
;     tau1,tau2,tau3,tau4,id2,id3,iop,stelib
;
;  where jo is the code:
;   Choose SFR: 0 = SSP (Delta Burst = zero length burst)'
;               1 = Exponential (enter Tau)'
;              -1 = Exponential (enter mu_SFR parameter)'
;               2 = Single Burst of finite length'
;               3 = Constant'
;               4 = Delayed'
;               5 = Linearly decreasing'
;               6 = Read SFR(t) from ASCII file'

