;+
; NAME:
;  fast_agnout
; PURPOSE:
;  Creates an output catalog providing the predicted AGN fluxes in
;  every observed band. Require FAST to have been run with FITAGN=1
;  and BEST_FIT=1
; CALLING_SEQUENCE:
;  fast_agnout,[param=param]
;  (must be run in the FAST directory)
; INPUT KEYWORDS:
;   param - use this parameter file (defaults to fast.param)
;   outputcat - name of the output catalog (defaults to CATALOG+'.AGNOUT.cat')
; MODIFICATION HISTORY:
;  Aug 2019, J. Aird - created routine
;-
pro fast_agnout,param=param,outputcat=outputcat

  if n_elements(param) eq 0 then param='fast.param'

  print,'Read parameter file: '+param
  READCOL,param,key,f='a,a',COMMENT='#',DELIM='#',/SILENT
  for i=0,n_elements(key)-1 do r=execute(key[i])

  if not keyword_set(OUTPUT_DIR) then OUTPUT_DIR='.' 
  if n_elements(outputcat) eq 0 then outputcat=CATALOG+'.AGNOUT.cat'
    
  ;; Read filters
  filters = FAST_READ_FILTERS(FILTERS_RES,lambda,CATALOG=CATALOG,$
                              SPECTRUM=SPECTRUM,no_filt=no_filt)
  
  ;; read the input catalog
  data = FAST_READ_OBSERVATIONS(lambda,CATALOG=CATALOG,SPECTRUM=$
                              SPECTRUM,AB_ZEROPOINT=AB_ZEROPOINT,$
                              C_INTERVAL=C_INTERVAL,NAME_ZPHOT=NAME_ZPHOT,H_CAT=H_CAT)
  n_gal  = n_elements(data.id)

  ;; read original header
  openr,lun,CATALOG+'.cat',/get_lun
  hline=''
  readf,lun,hline
  free_lun,lun
  h_orig=strsplit(hline,/extract)
  if h_orig[0] eq '#' then h_orig=h_orig[1:-1]

  ;; open output catalog and write header
  openw,outlun,OUTPUT_DIR+'/'+outputcat,/get_lun
  fl_ind=where(strmatch(h_cat,'F*[1234567890]') eq 1,n_filt)
  printf,outlun,'# ID '+strjoin(h_orig[fl_ind],'   ')

  for igal=0,n_gal-1 do begin
     ;; read output for this galaxy
     file=OUTPUT_DIR+'/BEST_FITS/'+data.name+'_'+strtrim(data.ID[igal],1)+'.AGN.input_res.fit'
     if file_test(file) eq 0 then begin
        print,'WARNING: no output for '+strtrim(data.ID[igal],1)+' ... filling in -99'
        printf,outlun,strtrim(data.ID[igal],1)+' '+strjoin(string(fltarr(n_filt)-99.,format='(G15.8)'),' ')
     endif else begin
        readcol,file,wl,fl,comment='#',/silent
        ;; convert fluxes back to input units
        fnu = (lambda*1d-8)^2/299792458d2 * fl *1d-19 *1d8 ; egs/s/cm^2/Hz
        fl_conv = fnu *10^((48.57+ab_zeropoint)/2.5)
        ;; print to file
        printf,outlun,strtrim(data.ID[igal],1)+' '+strjoin(string(fl_conv,format='(G15.8)'),' ')
     endelse

  endfor
  print,'Created catalog of predicted AGN fluxes in '+outputcat
  free_lun,outlun
  

end




  
