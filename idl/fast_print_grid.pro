PRO fast_print_grid,key,no_filt=no_filt,file=file,un=un,version=version,$
                    header=header

for i=0,n_elements(key)-1 do r=execute(key[i])

CASE IMF OF
    'ch': print_imf='Chabrier'
    'sa': print_imf='Salpeter'
    'kr': print_imf='Kroupa  '
ENDCASE

CASE LIBRARY OF
    'bc03': print_lib='Bruzual & Charlot (2003)'
    'cb07': print_lib='Charlot & Bruzual (2007)'
    'ma05': print_lib='Maraston (2005)         '
    'ma11': print_lib='Maraston (2011)         '
    'co11': print_lib='FSPS (Conroy et al.)    '
    else: print_lib ='Library unknown         '
 ENDCASE

if keyword_set(MY_SFH) then begin
   print_sfh = 'custom SFH:  '+MY_SFH
endif else begin
   CASE SFH OF
      'exp': print_sfh='Exponentially declining SFH: SFR ~ exp(-t/tau)'
      'del': print_sfh='Delayed exponential SFH: SFR ~ t exp(-t/tau)  '
      'tru': print_sfh='Truncated SFH: single burst of length tau     '
      else: print_sfh = SFH
   ENDCASE 
endelse

if KEYWORD_SET(DUST_LAW) then begin
   if DUST_LAW eq 'mw'       then print_dust='Milky Way dust attenuation law                    '
   if DUST_LAW eq 'calzetti' then print_dust='Calzetti (2000) dust attenuation law              '
   if DUST_LAW eq 'kc'       then print_dust='Kriek & Conroy (2013) average dust attenuation law'
   if DUST_LAW eq 'noll'     then print_dust='Noll et al. dust law, E_B='+$
                                             string(E_B,format='(f5.2)')+', delta='+$
                                             string(delta,format='(f5.2)')
endif else begin
   print,'ERROR: dust law not recognized'
   exit
endelse

if not KEYWORD_SET(file) then begin

    print,format='(a2,a8,a5,a25)','','Library:','',print_lib
    print,format='(a2,a4,a9,a47)','','SFH:','',print_sfh
    print,format='(a2,a12,a1,a9)','','Stellar IMF:','',print_imf
    print,format='(a2,a9,a1,a54)','','Dust law:','',print_dust
    print,format='(a2,a12,a1,'+strcompress(n_elements(METAL))+'(f7.4))','',$
      'metallicity:','',METAL
     print,format='(a2,a12,a1,f5.2,a4,f6.2,a13,f5.2)','','log(tau/yr):','',$
      LOG_TAU_MIN,'-',LOG_TAU_MAX, ', in steps of',LOG_TAU_STEP
    print,format='(a2,a12,a1,f5.2,a4,f6.2,a13,f5.2)','','log(age/yr):','',$
      LOG_AGE_MIN,'-',LOG_AGE_MAX, ', in steps of',LOG_AGE_STEP
    print,format='(a2,a4,a9,f5.2,a4,f6.2,a13,f5.2)','','A_V:','',A_V_MIN,$
      '-',A_V_MAX, ', in steps of',A_V_STEP
    print,format='(a2,a2,a11,f7.4,a2,f8.4,a13,f7.4)','','z:','',Z_MIN,$
      '-',Z_MAX, ', in steps of',Z_STEP

    if KEYWORD_SET(no_filt) then print,$
      format='(a2,a8,a4,'+strmid(strcompress(n_elements(no_filt)),1)+$
      '(i4))','','Filters:','',no_filt
    if KEYWORD_SET(SPECTRUM) then if SPECTRUM ne '' then $
      print,format='(a2,a9,a5,a'+strtrim(strlen(SPECTRUM),1)+')',$
      '','Spectrum:','',SPECTRUM

endif else begin

    printf,un,format='(a1,a1,a13,a5)','#','','FAST version:',version
    if file_test(CATALOG+'.cat') then printf,un,$
      "# Photometric catalog file: "+CATALOG+'.cat'
    if file_test(CATALOG+'.zout') then printf,un,$
      "# Photometric redshift file: "+CATALOG+'.zout'
    if file_test(SPECTRUM+'.spec') then printf,un,"# Spectrum file: "+$
      SPECTRUM+'.spec'
    if KEYWORD_SET(TEMP_ERR_FILE) then $
       printf,un,'# Template error function: '+file_basename(TEMP_ERR_FILE)
    if KEYWORD_SET(AB_ZEROPOINT) then $
       printf,un,format='(a1,a1,a6,a7,a6)','#','','AB ZP:','',$
              string(AB_ZEROPOINT,format='(f6.2)')

    printf,un,format='(a1,a1,a8,a5,a25)','#','','Library:','',print_lib
    printf,un,format='(a1,a1,a4,a9,a47)','#','','SFH:','',print_sfh
    printf,un,format='(a1,a1,a12,a1,a9)','#','','Stellar IMF:','',print_imf
    printf,un,format='(a1,a1,a9,a1,a54)','#','','Dust law:','',print_dust
    printf,un,format='(a1,a1,a12,a1,'+strcompress(n_elements(METAL))+$
      '(f7.4))','#','','metallicity:','',METAL
    printf,un,format='(a1,a1,a12,a1,f4.1,a5,f5.1,a13,f5.2)','#','',$
      'log(tau/yr):','',LOG_TAU_MIN,'-',LOG_TAU_MAX, ', in steps of',$
      LOG_TAU_STEP
    printf,un,format='(a1,a1,a12,a1,f4.1,a5,f5.1,a13,f5.2)','#','',$
      'log(age/yr):','',LOG_AGE_MIN,'-',LOG_AGE_MAX, ', in steps of',$
      LOG_AGE_STEP
    printf,un,format='(a1,a1,a4,a9,f4.1,a5,f5.1,a13,f5.2)','#','','A_V:',$
      '',A_V_MIN,'-',A_V_MAX, ', in steps of',A_V_STEP
    printf,un,format='(a1,a1,a2,a11,f7.4,a2,f8.4,a13,f7.4)','#','','z:','',$
      Z_MIN,'-',Z_MAX, ', in steps of',Z_STEP
    if  KEYWORD_SET(no_filt) then printf,un,$
      format='(a1,a1,a8,a4,'+strmid(strcompress(n_elements(no_filt)),1)+$
      '(i4))','#','','Filters:','',no_filt

endelse

if keyword_set(header) then begin

    ;... print also output information in file
    printf,un,'# ltau: log[tau/yr], lage: log[age/yr], '+$
      'lmass: log[mass/Msol], lsfr: log[sfr/(Msol/yr)], '+$
      'lssfr: log[ssfr*yr], la2t: log[age/tau]'
    printf,un,'# For sfr=0. lsfr is set to -99'
    
    par = ['z','ltau','metal','lage','Av','lmass','lsfr','lssfr','la2t']

    printf,un,format='($,a1,a6)','#','id'

    if N_SIM gt 0 then begin
        conf = ['']
        for i=0,n_elements(C_INTERVAL)-1 do begin
            CASE C_INTERVAL(i) of 
                68: conf=[conf,'l68_','u68_']
                95: conf=[conf,'l95_','u95_']
                99: conf=[conf,'l99_','u99_']
                else: begin
                    print, "ERROR choose confidence interval: 68%, 95% or 99%" 
                    RETURN
                endelse
            ENDCASE
        endfor
        n_conf = n_elements(C_INTERVAL) * 2
        conf   = conf[1:n_conf]
        for i=0,n_elements(par)-1 do $
          printf,un,format='($,'+strtrim(n_conf+1,1)+'(a10))',par(i),$
          conf+par(i)
    endif else begin
        for i=0,n_elements(par)-1 do printf,un,format='($,a10)',par(i)
     endelse
    
    ;; JAA additional outputs (optional)
    if keyword_set(calcsfr100) then begin
       para=['lsfr100']
    endif
    if keyword_set(calcmean) then begin
       push,para,['lmass_mean','lmass_med','lsfr_mean','lsfr_med']
       if keyword_set(calcsfr100) then begin
          push,para,['lsfr100_mean','lsfr100_med']
       endif
    endif

    push,para,'L2800'
    if keyword_set(calcmean) then push,para,'L2800_mean'

    
    ;; AGN outputs
    if keyword_set(fitagn) then begin
       para=[para,'agn_templ','fagn1m','fagn5000','fagn2800']
       if keyword_set(calcmean) then para=[para,'fagn1m_mean','fagn5000_mean','fagn2800_mean']
    endif
    

    if n_elements(para) gt 0 then for i=0,n_elements(para)-1 do printf,un,format='($,a15)',para(i)

    

    printf,un,format='(a10)','chi2'
endif

END
