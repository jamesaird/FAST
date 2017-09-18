function calzetti,lambda

  R_v      = 4.05

  k_lam = fltarr(N_ELEMENTS(lambda))
  x     = lambda/1.e4  ;wavelength in microns
  x1    = where(x le 0.63,n_x1)
  x2    = where(x ge 0.63,n_x2)

  if n_x1 ge 1 then k_lam(x1) = 2.659 * (-2.156 + 1.509/x(x1) - 0.198/(x(x1)^2) $
                                         + 0.011/(x(x1)^3)) + R_v
  if n_x2 ge 1 then k_lam(x2) = 2.659 * (-1.857 + 1.040/x(x2)) + R_v
  
  return,k_lam/R_v

end


function mw,lambda,R_V=R_V

  ;lambda in Angstrom
  if not KEYWORD_SET(R_V) THEN R_V = 3.1
  x     = 1.e4/lambda ;wavelength in microns
  y     = x-1.82
  n_lam = n_elements(lambda)
  a_lam = findgen(n_lam) ; A(lambda) / A(V)

  x1 = where(x le 1.1,n_x1)
  if n_x1 ge 1 then a_lam(x1) = 0.574*x(x1)^1.61 + (-0.527*x(x1)^1.61)/R_V

  x2 = where(x ge 1.1 and x le 3.3,n_x2)
  if n_x2 ge 1 then $
     a_lam(x2) = (1. + 0.17699*y(x2) - 0.50447*y(x2)^2 - 0.02427*y(x2)^3 + $
              0.72085*y(x2)^4 + 0.01979*y(x2)^5 - 0.77530*y(x2)^6 + $
              0.32999*y(x2)^7) + $
             (1.41338*y(x2) + 2.28305*y(x2)^2 + 1.07233*y(x2)^3 - $
              5.38434*y(x2)^4 - 0.62251*y(x2)^5 + 5.3026*y(x2)^6 - $
              2.09002*y(x2)^7)/R_V
  
  x3 = where(x ge 3.3 and x le 5.9,n_x3)
  if n_x3 ge 1 then $
     a_lam(x3) = 1.752-0.316*x(x3)-0.104/((x(x3)-4.67)^2+0.341) + $
                 (-3.09+1.825*x(x3)+1.206/((x(x3)-4.62)^2+0.263))/R_V

  x4 = where(x ge 5.9 and x le 8., n_x4)
  fa = -0.04473*(x(x4)-5.9)^2-0.009779*(x(x4)-5.9)^3
  fb =  0.2130*(x(x4)-5.9)^2+0.1207*(x(x4)-5.9)^3
  if n_x4 ge 1 then $
     a_lam(x4) = 1.752-0.316*x(x4)-0.104/((x(x4)-4.67)^2+0.341)+fa+$
                 (-3.09+1.825*x(x4)+1.206/((x(x4)-4.62)^2+0.263)+fb)$
                 /R_V

  x5 = where(x ge 8.,n_x5)
  if n_x5 ge 1 then a_lam(x5) = -1.073-0.628*(x(x5)-8.)+0.137*(x(x5)-8.)^2 + $
                                (13.67+4.257*(x(x5)-8.)-0.42*(x(x5)-8.)^2 + $
                                 0.374*(x(x5)-8.)^3)/R_V


  return,a_lam

end


function noll,lambda,E_B,delta

  if not keyword_set(E_B) or not keyword_set(delta) then begin
     print,'ERROR:  you have to define "E_B" and "delta" in the '+$
           'parameter file for the Noll et al. dust law'
     exit
  endif

  R_v     = 4.05
  width   = 350.
  clam    = 2175.
  D_lam   = E_B * lambda^2 * width^2 / ( (lambda^2 - clam^2)^2 + lambda^2 * width^2)
  alam2av = (calzetti(lambda)*R_v + D_lam) / R_v * (lambda/5500.)^delta

  return,alam2av

end



function fast_dust,synspec,wl,A_v,law=law,E_B=E_B,delta=delta

;...................................................................
;Description:
; Apply Calzetti extinction law 
;
;Parameters
; input:     synspec     Input model spectrum [[n_ages],[n_wl]]
;            wl          Wavelength array / Anstrom
;            A_v         Extinctions
; keywords:  calzetti    Choose Calzetti (2000) extinction curve
;            mw          Choose Milky Way extinction curve 
;                        (Cardelli et al. 1989)
;            f_mw        Combine two laws with f_mw the fraction of
;                        the Milky way, and (1-f_mw) the Calzetti
;                        fraction
; return:    corr_spec   Corrected model spectrum 
;                        [[n_wl],[n_ages],[n_Av]]
;
;...................................................................

n_ages   = (SIZE(synspec))[1]
n_wl     = (SIZE(synspec))[2]
n_Av     = N_ELEMENTS(A_v)
A_V      = REFORM(A_V,n_Av)

;...Determine attenuation curve A(lam)/A(V) for different dust laws
if keyword_set(law) then begin
   if law eq 'mw' then  Alam2Av = mw(wl)
   if law eq 'calzetti' then Alam2Av = calzetti(wl) 
   if law eq 'noll' then Alam2Av = noll(wl,E_B,delta)
   if law eq 'kc' then Alam2Av = noll(wl,1.,-0.1)
endif else Alam2Av = calzetti(wl)

;...Determine correction factor 
corr = 10^(-0.4 * A_v#Alam2Av) 

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

