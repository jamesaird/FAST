;...convert luminosity into flux 
;...3.826e33 / (3.0856776e24)^2 ->  3.826e1 / (3.0856776e8)^2

FUNCTION fast_lum2flux,z,H0=H0,lambda0=omega_l,omega_m=omega_m,$
                       AB_ZEROPOINT=AB_ZEROPOINT

if not KEYWORD_SET(AB_ZEROPOINT) then AB_ZEROPOINT = 25.; (-48.57)

zp_c_fl = 10^((48.57+AB_ZEROPOINT)/2.5)
dist    = LUMDIST(z,H0=H0,lambda0=omega_l,omega_m=omega_m,/SILENT) 
l2f     = 10^(ALOG10(3.826e1) + ALOG10(zp_c_fl) - $
              ALOG10(4.*!PI*(1.+z)*(dist*3.0856776e8)^2))

RETURN,l2f

END
