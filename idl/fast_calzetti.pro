function fast_calzetti,synspec,wl,A_v

;...................................................................
;Description:
; Apply Calzetti extinction law 
;
;Parameters
; input:     synspec     Input model spectrum [[n_ages],[n_wl]]
;            wl          Wavelength array / Anstrom
;            A_v         Extinctions
; return:    corr_spec   Corrected model spectrum 
;                        [[n_wl],[n_ages],[n_Av]]
;
;...................................................................

n_ages   = (SIZE(synspec))[1]
n_wl     = (SIZE(synspec))[2]
n_Av     = N_ELEMENTS(A_v)
A_V      = REFORM(A_V,n_Av)
R_v      = 4.05
errR_v   = 0.80

;...Determine k_lambda
k_lambda = fltarr(N_ELEMENTS(wl))
wl_tmp   = wl/1.e4              ;wavelength in microns
ff11     = 2.659*(-2.156+1.509/0.11-0.198/(0.11^2)+0.011/(0.11^3))+R_v
ff12     = 2.659*(-2.156+1.509/0.12-0.198/(0.12^2)+0.011/(0.12^3))+R_v
ff99     = 2.659*(-1.857+1.040/2.19)+R_v
ff100    = 2.659*(-1.857+1.040/2.2)+R_v

wl1      = where(wl_tmp gt 0.12 and wl_tmp le 0.63)
wl2      = where(wl_tmp gt 0.63 and wl_tmp le 2.2)
wl3      = where(wl_tmp lt 0.12,n_wl3)
wl4      = where(wl_tmp gt 2.2,n_wl4)

k_lambda      = fltarr(N_ELEMENTS(wl))
k_lambda(wl1) = 2.659 * (-2.156 + 1.509/wl_tmp(wl1) - 0.198/(wl_tmp(wl1)^2) + $
                           0.011/(wl_tmp(wl1)^3)) + R_v
k_lambda(wl2) = 2.659 * (-1.857 + 1.040/wl_tmp(wl2)) + R_v
if n_wl3 ge 1 then k_lambda(wl3) = ff11+(wl(wl3)-1100.)*(ff12-ff11)/100.
if n_wl4 ge 1 then k_lambda(wl4) = ff99+(wl(wl4)-21900.)*(ff100-ff99)/100. 

;...Modify synthetic spectrum
corr     = 10^((-0.4/R_v) * A_v#k_lambda)  < 1.0

if n_ages gt 1 then begin
    corr_spec = TRANSPOSE(REBIN(corr,n_Av,n_wl,n_ages),[1,2,0]) * $
      REBIN(TRANSPOSE(synspec),n_wl,n_ages,n_Av)
endif else begin
    corr_spec = REFORM(TRANSPOSE(REFORM(corr,n_Av,n_wl,n_ages),[1,2,0]),$
                       n_wl,n_ages,n_Av) * $
      REBIN(TRANSPOSE(synspec),n_wl,n_ages,n_Av)
endelse

RETURN,corr_spec 

end

