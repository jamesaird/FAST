;+
; NAME:
;   fast
; PURPOSE:
;   Performs SED fitting to synthetic galaxy templates (+AGN if
;   required). Most settings are defined in fast.param file.
; INPUT KEYWORDS: (can be used from within IDL only)
;   subset - gives a list of indices, FAST will only run for these
;                sources (useful for testing)
;   subname - name for the output .fout file, overrides param file
;   param - use this param file (otherwise defaults to fast.param).
; MODIFICATION HISTORY:
;  Jan 2016, J. Aird - added calculation of means (mass, SFR)
;  Mar 2016, J. Aird - added two-template fitting to include AGN
;  July 2017, J. Aird - fixed output for AB_ZEROPOINT!=25 - now output
;                       SED files always in 1e-19 erg/s/cm2/Ang regardless of
;                       input
;- 
PRO fast,subset=subset,subname=subname,param=param

version = '1.1'
date    = 'Sept 2017'
;                       '
print,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print,"+                                                              +"
print,"+      Fitting & Assessment of Synthetic Templates (FAST)      +"
print,"+                                                              +"
print,"+                   Mariska Kriek & Ivo Labbé                  +"
print,"+                                                              +"
print,"+                                                              +"
print,"+      Version "+strtrim(version,1)+": "+strtrim(date,1)+$
      "                        +"
print,"+                                                              +"
print,"+      Info:  astro.berkeley.edu/~mariska/FAST.html            +"
print,"+             Kriek, M., et al. 2009, ApJ, 700, 221            +"
print,"+                                                              +"
print,"+      If you use FAST in your publication, please cite        +"
print,"+      this paper!                                             +"
print,"+                                                              +"
print,"+            Updates by J. Aird, 2017 (v1.1).                  +"
print,"+                 see Aird et al. 2017, MNRAS, 465,3390;       +"
print,"+                     Aird et al. 2017b, arXiv:1705.01132      +"
print,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print,""


;--- READ PARAMETER FILE AND GRID IF PRESENT --------------------------
; 
;----------------------------------------------------------------------

t_start = systime(1) 
args=command_line_args(count=count)
if count gt 0 then param=args[0] else if n_elements(param) eq 0 then param='fast.param'
print,'Read parameter file: '+param
READCOL,param,key,f='a,a',COMMENT='#',DELIM='#',/SILENT
for i=0,n_elements(key)-1 do r=execute(key[i])
if not KEYWORD_SET(SPECTRUM) then AUTO_SCALE=0

if n_elements(subname) gt 0 then begin
   splog,'Fitting a subset, OUTPUT_FILE=',subname
   OUTPUT_FILE=subname
endif

filters = FAST_READ_FILTERS(FILTERS_RES,lambda,CATALOG=CATALOG,$
                            SPECTRUM=SPECTRUM,no_filt=no_filt)
z       = fast_z_arr(Z_MIN,Z_MAX,Z_STEP,Z_STEP_TYPE,CATALOG,NAME_ZPHOT=$
                     NAME_ZPHOT,N_SIM=N_SIM)
f_num   = MEAN(lambda)*100.
m_num   = TOTAL([AB_ZEROPOINT+50,Z_MIN,Z_MAX,mean(z),LOG_TAU_MIN,$
                 LOG_TAU_MAX,LOG_TAU_STEP,LOG_AGE_MIN,LOG_AGE_MAX,$
                 LOG_AGE_STEP,A_V_MIN,A_V_MAX,A_V_STEP,METAL] $
                * 1000. * PRIMES(13+n_elements(METAL)))

if DUST_LAW eq 'noll' then PDUST = DUST_LAW+'_'+string(E_B,format='(f5.3)')+$
                                   '_'+string(delta+5,format='(f5.3)') else $
                                      PDUST = DUST_LAW

if KEYWORD_SET(MY_SFH) then SFH = MY_SFH
tmp_name  = LIBRARY+'_'+RESOLUTION+'_'+SFH+'_'+IMF+'_'+PDUST+'_'+$
  strtrim(fix(m_num,type=3),1)+'_'+strtrim(fix(f_num,type=3),1)
if KEYWORD_SET(MIN_AGE) then begin
   ;; JAA - add min_age requirement to tmp_name
   if ZFORM gt 0 then begin
      tmp_name=tmp_name+'_'+string(AGE_FRAC*100,format='(I03)')+string(ZFORM,format='(I02)')
   endif else begin
      tmp_name=tmp_name+'_minage'
   endelse
endif

if not file_test(tmp_name) then file_mkdir, tmp_name


;--- MAKE OR READ MODEL CUBE ------------------------------------------
;
;----------------------------------------------------------------------

    
if not FILE_TEST(tmp_name+'/grid.save') then begin

    print,'Make grid: '+tmp_name+'/grid.save'
    FAST_PRINT_GRID,key,no_filt=no_filt

    ;step sizes cannot be zero
    if LOG_TAU_STEP eq 0. then LOG_TAU_STEP=0.1
    if LOG_AGE_STEP eq 0. then LOG_AGE_STEP=0.1
    if Z_STEP eq 0. then Z_STEP=0.01
    if A_V_STEP eq 0. then A_V_STEP=0.1
    
    n_z     = n_elements(z) 
    n_met   = n_elements(METAL)
    n_Av    = round((A_V_MAX - A_V_MIN) / A_V_STEP) + 1
    A_v     = A_V_MIN + A_V_STEP * findgen(n_Av)
    n_age   = round((LOG_AGE_MAX - LOG_AGE_MIN) / LOG_AGE_STEP) + 1
    log_age = LOG_AGE_MIN+LOG_AGE_STEP*findgen(n_age)
 
    if KEYWORD_SET(MY_SFH) then begin
       n_tau   = 1.
       log_tau = -99.
       print,'    Read in your custom SFH '+MY_SFH
    endif else begin
       n_tau   = round((LOG_TAU_MAX - LOG_TAU_MIN) / LOG_TAU_STEP) + 1
       log_tau = LOG_TAU_MIN+LOG_TAU_STEP*findgen(n_tau)
    endelse

    lum2fl = FAST_LUM2FLUX(z,H0=H0,lambda0=OMEGA_L,omega_m=$
                           OMEGA_M,AB_ZEROPOINT=AB_ZEROPOINT)
    
    if not KEYWORD_SET(NO_MAX_AGE) then begin
       la_univ = Alog10(galage(z,1.e4,H0=H0,lambda0=OMEGA_L,$
                               omega_m=OMEGA_M,/silent))
    endif else begin
        la_univ = REPLICATE(LOG_AGE_MAX,n_z)
     endelse

    if KEYWORD_SET(MIN_AGE) then begin
       ;; AGE_FRAC * time since z=zform
       if zform gt 0 then begin
          la_zmin=alog10( (galage(z,ZFORM,H0=H0,lambda0=OMEGA_L,omega_m=OMEGA_M,/silent)*AGE_FRAC)>0.1)
       endif else begin
          la_zmin=8.2-0.4*z   ; JAA update
       endelse
    endif else begin       
       la_zmin=fltarr(size(z,/dim))-1.0 ; young ages will not be exclued
    endelse
    

    grid      = fltarr(n_z,n_tau,n_met,n_age,n_Av,n_elements(lambda))
    mass_grid = fltarr(n_tau,n_met,n_age)
    sfr_grid  = fltarr(n_tau,n_met,n_age)
    galfl1m   = fltarr(n_z,n_tau,n_met,n_age,n_av) ; galaxy template flux denstiy at rest 1micron
    galf2800   = fltarr(n_z,n_tau,n_met,n_age,n_av) ; galaxy template flux density at rest 2800Ang
    galf5000   = fltarr(n_z,n_tau,n_met,n_age,n_av) ; galaxy template flux density at rest 5000Ang

    for i=0,n_met-1 do for k=0,n_tau-1 do begin

        if not KEYWORD_SET(MY_SFH) then begin
           name    = LIBRARY_DIR+'/ised_'+SFH+'.'+RESOLUTION+'/'+$
                     LIBRARY+'_'+RESOLUTION+'_'+IMF+$
                     repstr(string(METAL[i],f='(g0)'),'0.','_z')+'_ltau'+$
                     strtrim(string(log_tau(k),f='(f5.1)'),1)+'.ised'
        endif else begin
           name    = LIBRARY_DIR+'/'+LIBRARY+'_'+RESOLUTION+'_'+IMF+$
                     repstr(string(METAL[i],f='(g0)'),'0.','_z')+'_'+MY_SFH+$
                     '.ised'
        endelse

        print,'    Read: '+file_basename(name)
        
        if not FILE_TEST(name) then begin
           print,'ERROR:" Cannot find '+name
           print,'        Check path name, resolution and whether file exist' $
                 & exit
        endif
        synspec = FAST_READ_MODEL(name,log_age,wl_rest,sfr,mass)
        synspec = FAST_DUST(temporary(synspec),wl_rest,A_v,law=DUST_LAW,$
                            E_B=E_B,delta=delta)

        for r=0,n_z-1 do begin
            g_age = where(log_age le la_univ(r) and log_age ge la_zmin(r),n_g_age)
            synspec2 = FAST_MADAU(TEMPORARY(synspec(*,g_age,*)),wl_rest,z(r))
            grid(r,k,i,g_age,*,*) = lum2fl(r) * $
              FAST_INTEGRATE(wl_rest*(1.+z(r)),synspec2,filters,$
                             FILTER_FORMAT=FILTER_FORMAT)
            i1m = round(findex(wl_rest,1e4))
            galfl1m[r,k,i,g_age,*] = lum2fl[r] * reform(synspec2[i1m,*,*],n_g_age,n_av)
            i2800 = round(findex(wl_rest,2800.))
            galf2800[r,k,i,g_age,*] = lum2fl[r] * reform(synspec2[i2800,*,*],n_g_age,n_av)
            i5000 = round(findex(wl_rest,5000.))
            galf5000[r,k,i,g_age,*] = lum2fl[r] * reform(synspec2[i5000,*,*],n_g_age,n_av)
        endfor
        mass_grid(k,i,*) = mass & sfr_grid(k,i,*) = sfr
    endfor
    synspec  = [0]
    synspec2 = [0]

    print,format='(a11,a'+strtrim(strlen(tmp_name)+10,1)+',f'+$
      strtrim(45-strlen(tmp_name),1)+'.1,a4)','Save grid: ',$
      tmp_name+'/grid.save',(systime(1)-t_start)/60.,'min'
    SAVE,z,log_tau,metal,log_age,A_v,grid,mass_grid,sfr_grid,lum2fl,galfl1m,galf2800,galf5000,$
      filename=tmp_name+'/grid.save'
    
    if not file_test('00README') then begin
        close,2 & openw,2,'00README'
        printf,2,'Decoding of grid names'
        printf,2,'' & printf,2,tmp_name,':'
        FAST_PRINT_GRID,key,no_filt=no_filt,/file,un=2,version=version 
        close,2
    endif else begin
        close,2 & openu,2,'00README',/append
        printf,2,'' & printf,2,tmp_name,':'
        FAST_PRINT_GRID,key,no_filt=no_filt,/file,un=2,version=version
        close,2
    endelse

endif else begin
    
    print,format='(a11,a'+strtrim(strlen(tmp_name)+10,1)+',f'+$
      strtrim(35-strlen(tmp_name),1)+'.2,a4)','Read grid: ',$
      tmp_name+'/grid.save'
    FAST_PRINT_GRID,key,no_filt=no_filt
    RESTORE,tmp_name+'/grid.save'

endelse

;; JAA - calculate SFR100 for grid
if keyword_set(calcsfr100) then begin
   sfr100_grid=fast_calcsfr100(sfr_grid,log_age,log_tau,/grid,sfh=SFH)
endif

;; JAA - prior on dust
if keyword_set(DUSTPRIOR) then begin
   p_av=a_v^2.0*exp(-a_v*2.0/avpeak)
   p_av=transpose(rebin(p_av,n_elements(a_v),1,n_elements(log_tau),$
                        n_elements(metal),n_elements(log_age)),[1,2,3,4,0])
endif

;----- READ AGN MODELS -----------------------------------------------
;
;---------------------------------------------------------------------
if keyword_set(fitagn) then begin
   if n_elements(AGNLIBRARY) eq 0 then $
      message,'Must specify AGNLIBRARY in parameter file if FITAGN=1'
   mnum_agn = string(floor(AB_ZEROPOINT*100),format='(I04)')+$
              string(Z_MIN*100,format='(I04)')+string(Z_MAX*100,format='(I04)')+$
              '_'+strtrim(fix(f_num,type=3),1)
   agnlibfile=AGNLIBRARY+'_'+mnum_agn+'.save'
   if file_test(agnlibfile) then begin
      restore,agnlibfile
      splog,'Restore AGN library: '+agnlibfile
   endif else begin
      splog,'generating the AGN model library'
      lum2fl = FAST_LUM2FLUX(z,H0=H0,lambda0=OMEGA_L,omega_m=$
                             OMEGA_M,AB_ZEROPOINT=AB_ZEROPOINT)
      ;; read the list of AGN templates
      readcol,LIBRARY_DIR+'/AGN/'+AGNLIBRARY+'.list',AGNtempl_num,AGNtempl_name,format='L,A'
      n_agn=n_elements(AGNtempl_name)
      n_z=n_elements(z)
      AGNgrid      = fltarr(n_z,n_agn,n_elements(lambda))
      AGNf2800   = fltarr(n_z,n_agn) ;; agn template flux density at rest 2800Ang
      AGNf5000   = fltarr(n_z,n_agn) ;; agn template flux density at rest 5000Ang
      for iagn=0,n_agn-1 do begin
         splog,'Reading AGN template: '+AGNtempl_name[iagn]
         readcol,LIBRARY_DIR+'/AGN/'+AGNtempl_name[iagn],wl_rest,synspec,$
                 format='F,F',comment='#'
         synspec = reform(synspec,n_elements(wl_rest),1,1)
         for r=0,n_z-1 do begin
            synspec2 = fast_madau(synspec,wl_rest,z[r])
            ;; normalize to 1 at rest-frame 1micron
            synspec2 = synspec2 / interpol(synspec2[*,0,0],wl_rest, 1e4)
            synspec2 = reform(synspec2,n_elements(wl_rest),1,1)
            AGNgrid[r,iagn,*]= FAST_INTEGRATE(wl_rest*(1+z[r]),synspec2,filters,$
                                              FILTER_FORMAT=FILTER_FORMAT)
            AGNf2800[r,iagn] = interpol(synspec2,wl_rest,2800.0)
            AGNf5000[r,iagn] = interpol(synspec2,wl_rest,2800.0)
         endfor
      endfor
      save,AGNgrid,AGNf2800,AGNf5000,filename=agnlibfile
   endelse
endif
;--- READ AND FIT DATA ------------------------------------------------
;
;----------------------------------------------------------------------

data = FAST_READ_OBSERVATIONS(lambda,CATALOG=CATALOG,SPECTRUM=$
                              SPECTRUM,AB_ZEROPOINT=AB_ZEROPOINT,$
                              C_INTERVAL=C_INTERVAL,NAME_ZPHOT=NAME_ZPHOT)
n_gal  = n_elements(data.id)
if n_gal lt 100 then n_step = 10. else n_step=10000.
do_p   = n_gal/n_step

print,format='(f66.1,a4)',(systime(1)-t_start)/60.,'min'
print,'Fit '+strtrim(n_gal,1)+' object(s)'

if KEYWORD_SET(TEMP_ERR_FILE) then if FILE_TEST(TEMP_ERR_FILE) then $
  print,"Apply template error function to photometry"
close,1
if not keyword_set(OUTPUT_DIR) then OUTPUT_DIR='.' 
if not keyword_set(OUTPUT_FILE) then OUTPUT_FILE=data.name
if not file_test(OUTPUT_DIR) then file_mkdir, OUTPUT_DIR
if not file_test(tmp_name+'/'+data.name+'.param') then FILE_COPY,param,$
  tmp_name+'/'+data.name+'.param'

openw,1,OUTPUT_DIR+'/'+OUTPUT_FILE+'.fout'
FAST_PRINT_GRID,key,no_filt=no_filt,/file,un=1,version=version,/header


if n_elements(subset) eq 0 then subset=lindgen(n_gal)
  
   
for j=0l,n_elements(subset)-1 do begin
   i=subset[j]
   if i lt n_gal then begin
   ;for i=0l,199 do begin        ;n_gal-1 do begin

    if i ne 0 and i eq round(do_p) then begin
        print,format='($,a2,i2,6(a1))','',round(do_p/n_gal*100.),'%',$
          string(8b),string(8b),string(8b),string(8b),string(8b)
        do_p = TEMPORARY(do_p) + n_gal/n_step    
     endif
    FAST_FIT,data.id(i),data.flux(*,i),data.eflux(*,i),data.zspec(i),$
             data.zphot(*,i),grid,mass_grid,sfr_grid,Z,LOG_TAU,METAL,LOG_AGE,A_V,$
             sfr100_grid=sfr100_grid,CALCMEAN=CALCMEAN,$
             DUSTPRIOR=DUSTPRIOR,P_AV=p_av,$
             data.name,filters,key,tmp_name,C_INTERVAL=C_INTERVAL,N_SIM=N_SIM,$
             AUTO_SCALE=AUTO_SCALE,TEMP_ERR_FILE=TEMP_ERR_FILE,SAVE_CHI_GRID=$
             SAVE_CHI_GRID,BEST_FIT=BEST_FIT,scale_bands=scale_bands,$
             OUTPUT_DIR=OUTPUT_DIR, LAMBDA=LAMBDA,$
             AB_ZEROPOINT=AB_ZEROPOINT,$
             FITAGN=FITAGN,AGNmodel=AGNgrid,$
             galfl1m=galfl1m,galf2800=galf2800,AGNf2800=AGNf2800,galf5000=galf5000,agnf5000=agnf5000
    if i eq 0 and KEYWORD_SET(AUTO_SCALE) then begin
        print,"Scale spectrum using photometric bands:",no_filt(scale_bands)
        print,"  NOTE: Bands used for scaling are not included in the fit"
     endif
    if i mod 10 eq 0 then print,i
    endif
endfor
close,1
if file_test(tmp_name+'/'+data.name+'.fout') then FILE_DELETE,tmp_name+'/'+$
  data.name+'.fout'
FILE_COPY,OUTPUT_DIR+'/'+OUTPUT_FILE+'.fout',tmp_name+'/'+data.name+'.fout'

print,format='(a4,f62.1,a4)','Done',(systime(1)-t_start)/60.,'min'
print,''
junk = check_math()

END

