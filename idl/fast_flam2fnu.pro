;+
; NAME:
;  fast_flam2fnu
; PURPOSE:
;  Convert flux in flambda(ergs s^-1 cm^-2 Angstrom^-1) to fnu (erg
;  s^-1 cm^-2 Hz^-1)
; CALLING SEQUENCE:
;  fnu=ja_flam2fnu(lambda,flam)
; INPUTS:
;  lambda - in Angstrom
;  flambda - in erg s^-1 cm^-2 Angstrom^-1
; RETURN VALUE:
;  fnu in erg s^-1 cm^-2 Hz^-1
;
;-
function fast_flam2fnu,lambda,flam

  c2=299792458d2                ; c in cm^2
  x=lambda*1d-8                      ; cm
  fnu=x^2/c2*flam*1d8
  return,fnu
end

  
