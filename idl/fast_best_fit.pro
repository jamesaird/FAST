PRO fast_best_fit, DIR, NAME, ID, KEY, BEST, SCALE, LAMBDA, BEST_LFIT,$
                   agn_templ=agn_templ, agn_scale=agn_scale,$
                   agn_lfit=agn_lfit,AB_ZEROPOINT=AB_ZEROPOINT

junk    = CHECK_MATH()
for i=0,n_elements(key)-1 do r=execute(key[i])

if not KEYWORD_SET(MY_SFH) then begin
   ised    = LIBRARY_DIR+'/ised_'+SFH+'.'+RESOLUTION+'/'+$
             LIBRARY+'_'+RESOLUTION+'_'+IMF+$
             repstr(string(BEST[2],f='(g0)'),'0.','_z')+'_ltau'+$
             strtrim(string(BEST[1],f='(f5.1)'),1)+'.ised'
endif else begin
   ised    = LIBRARY_DIR+'/'+LIBRARY+'_'+RESOLUTION+'_'+IMF+$
             repstr(string(BEST[2],f='(g0)'),'0.','_z')+'_'+MY_SFH+$
             '.ised'
endelse
synspec = FAST_READ_MODEL(ised,BEST(3),wl)
intspec = FAST_READ_MODEL(ised,BEST(3),wl)
synspec = FAST_DUST(temporary(synspec),wl,BEST(4),law=DUST_LAW,$
                    E_B=E_B,delta=delta)
junk    = CHECK_MATH()
synspec = FAST_LUM2FLUX(BEST(0)) * FAST_MADAU(TEMPORARY(synspec),wl,BEST(0))
intspec = FAST_LUM2FLUX(BEST(0)) * TEMPORARY(intspec)
synspec = REFORM(SCALE * TEMPORARY(synspec)  / (10^((25.0+48.57)/2.5 - 19.)))
intspec = REFORM(SCALE * TEMPORARY(intspec)  / (10^((25.0+48.57)/2.5 - 19.)))
wl_obs  = (1.+BEST(0))*wl

if not file_test(DIR+'/BEST_FITS') then file_mkdir, DIR+'/BEST_FITS'

Close,5
openw,5,DIR+'/BEST_FITS/'+NAME+'_'+strtrim(ID,1)+'.fit'
printf,5,'# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)'
for i=0,n_elements(wl)-1 do printf,5,wl_obs(i),synspec(i),intspec(i)
close,5

best_lfit = best_lfit  / (10^((AB_ZEROPOINT+48.57)/2.5 - 19.))
openw,6,DIR+'/BEST_FITS/'+NAME+'_'+strtrim(ID,1)+'.input_res.fit'
printf,6,'# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)'
for i=0,n_elements(lambda)-1 do printf,6,lambda[i],best_lfit(i)
close,6

if n_elements(agn_templ) gt 0 then begin
   if agn_templ eq 0 then begin
      ;; best fit is a pure galaxy template - write out 0s for AGN
      Close,7
      openw,7,DIR+'/BEST_FITS/'+NAME+'_'+strtrim(ID,1)+'.AGN.fit'
      printf,7,'# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)'
      for i=0,n_elements(wl_obs)-1 do printf,7,wl_obs[i],0.0,0.0
      close,7

      openw,6,DIR+'/BEST_FITS/'+NAME+'_'+strtrim(ID,1)+'.AGN.input_res.fit'
      printf,6,'# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)'
      for i=0,n_elements(lambda)-1 do printf,6,lambda[i],0.0
      close,6
   endif else begin
      ;; read the list of AGN templates
      readcol,LIBRARY_DIR+'/AGN/'+AGNLIBRARY+'.list',AGNtempl_num,AGNtempl_name,format='L,A',/silent
      iagn=where(AGNtempl_num eq agn_templ,n)
      if n eq 0 then message,'Cannot find the best AGN template in the library'
      ;; read the best fit 
      readcol,LIBRARY_DIR+'/AGN/'+AGNtempl_name[iagn],wl_rest,synspec,$
              format='F,F',comment='#',/silent
      ;; normalize to 1 at rest-frame 1micron
      synspec = synspec / interpol(synspec,wl_rest, 1e4)
	    intspec = synspec / interpol(synspec,wl_rest, 1e4)
      synspec = reform(synspec,n_elements(wl_rest),1,1)
	    intspec = reform(intspec,n_elements(wl_rest),1,1)
      ;; redden
      synspec = FAST_MADAU(TEMPORARY(synspec),wl_rest,BEST(0))
      synspec = REFORM(AGN_SCALE * TEMPORARY(synspec) / (10^((AB_ZEROPOINT+48.57)/2.5 - 19.)))
	    intspec = REFORM(AGN_SCALE * TEMPORARY(intspec) / (10^((AB_ZEROPOINT+48.57)/2.5 - 19.)))
      wl_obs  = (1.+BEST(0))*wl_rest
      
      Close,7
      openw,7,DIR+'/BEST_FITS/'+NAME+'_'+strtrim(ID,1)+'.AGN.fit'
      printf,7,'# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)'
      for i=0,n_elements(wl_obs)-1 do printf,7,wl_obs[i],synspec[i],intspec[i]
      close,7

      agn_lfit = agn_lfit / (10^((AB_ZEROPOINT+48.57)/2.5 - 19.))
      openw,6,DIR+'/BEST_FITS/'+NAME+'_'+strtrim(ID,1)+'.AGN.input_res.fit'
      printf,6,'# wl fl (x 10^-19 ergs s^-1 cm^-2 Angstrom^-1)'
      for i=0,n_elements(lambda)-1 do printf,6,lambda[i],agn_lfit[i]
      close,6

      
   endelse
endif


END
