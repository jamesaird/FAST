;...................................................................
;Description:
; Determine flux for each filter / spectral elements for all input 
; synthetic models
;
;Parameters
; input:     wl_model    wavelength array model [wl/Angstrom]
;            flux_model  flux (in F_lambda) for all models 
;                        [[n_wl],[n_age],[n_Av]]
;            filters     array with spectral and filter information
;                        made by FAST_READ_FILTERS 
;                        [# filt,wl/Angstrom,tr]
; output:    out_mod     flux (in F_lambda) for all filters / 
;                        spectral elements for the full grid
;                        [[n_age],[n_Av],[n_filt]]
;  
;Notes
; + Integrate in lambda * F_lambda (instructions Ivo)
; + In case of broadband photometry, integrate over both the 
;   wavelength arrays of the model and filters
; + For the spectrum this doesn't work as some wavelengths can be
;   masked. So for the spectrum we only integrate the model over the
;   wavelength array set by the spectrum.
;...................................................................

function trapz1, x, y
  n = (size(x))[1]
  return, total((x[1:n-1] - x[0:n-2])*(y[1:n-1]+y[0:n-2]),1)/2.
end

function trapz1_3d, x, y
  n = (size(x))[1]
  return, total((x[1:n-1,*,*] - x[0:n-2,*,*])*(y[1:n-1,*,*] + $
                                               y[0:n-2,*,*]),1)/2.
end


FUNCTION fast_integrate,wl_model,flux_model,filters,FILTER_FORMAT=$
                        FILTER_FORMAT

if not KEYWORD_SET(FILTER_FORMAT) then FILTER_FORMAT = 1

n_dat    = N_ELEMENTS(UNIQ(filters(0,*))) 
n_age    = (SIZE(flux_model))[2]
n_Av     = (SIZE(flux_model))[3]
out_mod  = fltarr(n_age,n_Av,n_dat)

for j=0,n_dat-1 do begin

    ;...Select filtercurve
    filt_ind  = where(filters(0,*) eq j,n_el_f)
    wl_filt   = REFORM(filters(1,filt_ind))
    tr_filt   = REFORM(filters(2,filt_ind))

    ;...Cut relevant part from model spectrum
    l_wl      = where(abs(wl_model-min(wl_filt)) eq $
                      min(abs(wl_model-min(wl_filt))))
    h_wl      = where(abs(wl_model-max(wl_filt)) eq $
                      min(abs(wl_model-max(wl_filt))))
    mod_ind   = where(wl_model ge MIN(wl_model(l_wl)) and wl_model $
                      le MAX(wl_model(h_wl)),n_el_m)
    mod_ind   = [min(mod_ind)-1,mod_ind,max(mod_ind)+1]
    flux_mod  = REFORM(flux_model(mod_ind,*,*),n_elements(mod_ind),$
                       (size(flux_model))[2],(size(flux_model))[3])
    wl_mod    = wl_model(mod_ind)
    
    ;...Determine model fluxes at filter wavelengths 
    wl_int    = INTERPOL(INDGEN(n_el_m+2),wl_mod,wl_filt)
    flux_filt = INTERPOLATE(flux_mod,wl_int,INDGEN(n_age),$
                            INDGEN(n_Av),/grid) 

    if MIN(filters(3,filt_ind)) eq 1 then begin

        ;...Photometric grid points
        wl_tot    = [wl_mod,wl_filt]
        isort     = SORT(wl_tot)
        flux_tot  = ([flux_mod,flux_filt])[isort,*,*]
        wl_tot    = ([wl_mod,wl_filt])[isort]
        tr_tot    = ([INTERPOL(tr_filt,wl_filt,wl_mod),tr_filt])[isort] > 0
        
        if FILTER_FORMAT eq 1 then tr_wl_tot = wl_tot*tr_tot else $
          tr_wl_tot = wl_tot*wl_tot*tr_tot
        
        ;...Interpolate filter curve
        n_wl_tot       = n_elements(wl_tot)
        area           = REPLICATE(trapz1(wl_tot,tr_wl_tot),n_age,n_Av)
        wl_tot         = REBIN(TEMPORARY(wl_tot),n_wl_tot,n_age,n_Av)
        tr_wl_tot      = REBIN(TEMPORARY(tr_wl_tot),n_wl_tot,n_age,n_Av)
        out_mod[0,0,j] = trapz1_3d(wl_tot,tr_wl_tot*flux_tot) / area

    endif else begin

        if TOTAL(tr_filt) eq 0 then begin
            print,"ERROR: Total transmission of a spectral bin cannot be zero"
            exit
        endif

        ;...Spectroscopic grid points
        out_mod[0,0,j] = TOTAL(REBIN(tr_filt,n_el_f,n_age,n_Av)*flux_filt,1) $
          / TOTAL(tr_filt)

    endelse

endfor

RETURN, out_mod

END

