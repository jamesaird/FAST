;......................................................................
;description:
; Correct model IGM absorption (Madau prescriptions)
;
;parameters:
; input    spectrum        model spectra for different ages 
;                          [[n_wl],[n_age],[n_Av]]
;          wl              rest-frame wavelength array in Angstrom
;          z               array of input redshifts
; return   out_spec        corrected model
;                          [[n_wl],[n_age],[n_Av]]
;
; Note: check with Subaru verhaal, less than 912?
;......................................................................


FUNCTION fast_madau,spectrum,wl_rest,z

if (size(spectrum))[0] eq 1 then begin
   spectrum = REFORM(spectrum,n_elements(spectrum),1,1)
endif
if (size(spectrum))[0] eq 2 then begin
   spectrum = REFORM(spectrum,(size(spectrum))[1],(size(spectrum))[2],1)
endif

;Calculate deficients
w_st     = 1050.*(1.+z)
w_en     = 1170.*(1.+z)
w_step   = 100.
wtemp    = w_st+(w_en-w_st)/w_step*findgen(w_step)
ptau     = fltarr(w_step)
ptau     = EXP(-3.6E-3*(wtemp/1216.)^3.46)
da       = (1./(120.*(1.+z)))*TOTAL(ptau)*(w_en-w_st)/w_step
w_st     = 920.*(1.+z)
w_en     = 1015.*(1.+z)
wtemp    = w_st+(w_en-w_st)/w_step*findgen(w_step)
ptau     = EXP(-1.7E-3*(wtemp/1026.)^3.46 - 1.2E-3*(wtemp/972.5)^3.46 $
               - 9.3E-4*(wtemp/950.)^3.46)
db       = (1./(95.*(1.+z)))*TOTAL(ptau)*(w_en-w_st)/w_step

;Do corrections
outspec  = spectrum
wl_defa  = where(wl_rest ge 1026. and wl_rest lt 1216.)
wl_defb  = where(wl_rest lt 1026.)  
wl_zero  = where(wl_rest lt 912)        
outspec(wl_defa,*,*) = da*spectrum(wl_defa,*,*)
outspec(wl_defb,*,*) = db*spectrum(wl_defb,*,*)
outspec(wl_zero,*,*) = 0.
outspec  = REFORM(TEMPORARY(outspec),(size(spectrum))[1],$
                  (size(spectrum))[2],(size(spectrum))[3])

RETURN,outspec

END
