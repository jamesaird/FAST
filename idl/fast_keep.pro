PRO fast,E_b=E_b,delta=delta

version = '1.0'
date    = 'April 26 2012  '
;                       '
print,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print,"+                                                              +"
print,"+      Fitting & Assessment of Synthetic Templates (FAST)      +"
print,"+                                                              +"
print,"+                   Mariska Kriek & Ivo Labbé                  +"
print,"+                                                              +"
print,"+                                                              +"
print,"+      Version "+strtrim(version,1)+": "+strtrim(date,1)+$
      "                            +"
print,"+                                                              +"
print,"+      Info:  www.cfa.harvard.edu/~mkriek/FAST.html            +"
print,"+             Kriek, M., et al. 2009, ApJ, 700, 221            +"
print,"+                                                              +"
print,"+      If you use FAST in your publication, please cite        +"
print,"+      this paper!                                             +"
print,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print,""


;--- READ PARAMETER FILE AND GRID IF PRESENT --------------------------
; 
;----------------------------------------------------------------------

t_start = systime(1) 

args=command_line_args(count=count)
if count gt 0 then param=args[0] else param='fast.param'
print,'Read parameter file: '+param
READCOL,param,key,f='a,a',COMMENT='#',DELIM='#',/SILENT
for i=0,n_elements(key)-1 do r=execute(key[i])
if not KEYWORD_SET(SPECTRUM) then AUTO_SCALE=0

filters = FAST_READ_FILTERS(FILTERS_RES,lambda,CATALOG=CATALOG,$
                            SPECTRUM=SPECTRUM,no_filt=no_filt)
z       = fast_z_arr(Z_MIN,Z_MAX,Z_STEP,Z_STEP_TYPE,CATALOG,NAME_ZPHOT=$
                     NAME_ZPHOT,N_SIM=N_SIM)
f_num   = MEAN(lambda)*100.
m_num   = TOTAL([AB_ZEROPOINT+50,Z_MIN,Z_MAX,mean(z),LOG_TAU_MIN,$
                 LOG_TAU_MAX,LOG_TAU_STEP,LOG_AGE_MIN,LOG_AGE_MAX,$
                 LOG_AGE_STEP,A_V_MIN,A_V_MAX,A_V_STEP,METAL] $
                * 1000. * PRIMES(13+n_elements(METAL)))

if DUST eq 'noll' then PDUST = DUST+'_'+string(E_B,format='(f5.3)')+$
                               '_'+string(delta+5,format='(f5.3)') else $
                                  PDUST = DUST

if KEYWORD_SET(MY_SFH) then SFH = MY_SFH
tmp_name  = LIBRARY+'_'+RESOLUTION+'_'+SFH+'_'+IMF+'_'+PDUST+'_'+$
  strtrim(fix(m_num,type=3),1)+'_'+strtrim(fix(f_num,type=3),1)
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

    grid      = fltarr(n_z,n_tau,n_met,n_age,n_Av,n_elements(lambda))
    mass_grid = fltarr(n_tau,n_met,n_age)
    sfr_grid  = fltarr(n_tau,n_met,n_age)

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

        synspec = FAST_DUST(temporary(synspec),wl_rest,A_v,law=DUST,E_b=E_b,delta=delta)
        for r=0,n_z-1 do begin
            g_age = where(log_age le la_univ(r))
            synspec2 = FAST_MADAU(TEMPORARY(synspec(*,g_age,*)),wl_rest,z(r))
            grid(r,k,i,g_age,*,*) = lum2fl(r) * $
              FAST_INTEGRATE(wl_rest*(1.+z(r)),synspec2,filters,$
                             FILTER_FORMAT=FILTER_FORMAT)
        endfor
        mass_grid(k,i,*) = mass & sfr_grid(k,i,*) = sfr
    endfor
    synspec  = [0]
    synspec2 = [0]

    print,format='(a11,a'+strtrim(strlen(tmp_name)+10,1)+',f'+$
      strtrim(35-strlen(tmp_name),1)+'.1,a4)','Save grid: ',$
      tmp_name+'/grid.save',(systime(1)-t_start)/60.,'min'
    SAVE,z,log_tau,metal,log_age,A_v,grid,mass_grid,sfr_grid,lum2fl,$
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


;--- READ AND FIT DATA ------------------------------------------------
;
;----------------------------------------------------------------------

data = FAST_READ_OBSERVATIONS(lambda,CATALOG=CATALOG,SPECTRUM=$
                              SPECTRUM,AB_ZEROPOINT=AB_ZEROPOINT,$
                              C_INTERVAL=C_INTERVAL,NAME_ZPHOT=NAME_ZPHOT)
n_gal  = n_elements(data.id)
if n_gal lt 100 then n_step = 10. else n_step=100.
do_p   = n_gal/n_step

print,format='(f56.1,a4)',(systime(1)-t_start)/60.,'min'
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

for i=0l,n_gal-1 do begin

    if i ne 0 and i eq round(do_p) then begin
        print,format='($,a2,i2,6(a1))','',round(do_p/n_gal*100.),'%',$
          string(8b),string(8b),string(8b),string(8b),string(8b)
        do_p = TEMPORARY(do_p) + n_gal/n_step    
    endif

    FAST_FIT,data.id(i),data.flux(*,i),data.eflux(*,i),data.zspec(i),$
      data.zphot(*,i),grid,mass_grid,sfr_grid,Z,LOG_TAU,METAL,LOG_AGE,A_V,$
      data.name,filters,key,tmp_name,C_INTERVAL=C_INTERVAL,N_SIM=N_SIM,$
      AUTO_SCALE=AUTO_SCALE,TEMP_ERR_FILE=TEMP_ERR_FILE,SAVE_CHI_GRID=$
      SAVE_CHI_GRID,BEST_FIT=BEST_FIT,scale_bands=scale_bands,$
      OUTPUT_DIR=OUTPUT_DIR,DISPERSION=DISPERSION

    if i eq 0 and KEYWORD_SET(AUTO_SCALE) then begin
        print,"Scale spectrum using photometric bands:",no_filt(scale_bands)
        print,"  NOTE: Bands used for scaling are not included in the fit"
    endif

endfor
close,1
if file_test(tmp_name+'/'+data.name+'.fout') then FILE_DELETE,tmp_name+'/'+$
  data.name+'.fout'
FILE_COPY,OUTPUT_DIR+'/'+OUTPUT_FILE+'.fout',tmp_name+'/'+data.name+'.fout'

print,format='(a4,f52.1,a4)','Done',(systime(1)-t_start)/60.,'min'
print,''
junk = check_math()

END

