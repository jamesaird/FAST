function fast_auto_scale,flux,eflux,filters,scale_bands=scale_b,$
                         fit_bands=fit_bands,fscale=scale

n_dat   = n_elements(flux)
phot    = filters(0,UNIQ(filters(0,where(filters(3,*) eq 1))))
n_phot  = n_elements(phot)
spec    = max(phot)+1+indgen(n_dat-n_phot)
wl_spec = REFORM(filters(1,where(filters(3,*) eq 0)))
fl_spec = fltarr(n_elements(wl_spec))
j       = 0
for i=0,n_elements(spec)-1 do begin
    bin = where(filters(0,*) eq spec(i),n_bin)
    fl_spec(j:(j+n_bin-1)) = REPLICATE(flux(spec(i)),n_bin)
    j   = TEMPORARY(j)+n_bin
endfor

int_sp  = fltarr(n_dat)
for j=0,n_phot-1 do begin
    filt_ind  = where(filters(0,*) eq phot(j))
    good_tr   = where(filters(2,filt_ind)/max(filters(2,filt_ind)) gt 0.04)
    if min((filters(1,filt_ind))[good_tr]) ge min(wl_spec) and $
      max((filters(1,filt_ind))[good_tr]) le max(wl_spec) and $
      flux(phot(j)) ne -99 then begin
        wl_filt = REFORM(filters(1,filt_ind))
        tr_filt = REFORM(filters(2,filt_ind))
        tr_spec = INTERPOL(tr_filt,wl_filt,wl_spec) > 0
        fl_filt = INTERPOL(fl_spec,wl_spec,wl_filt) > 0
        i_sort  = SORT([wl_filt,wl_spec])
        wl_tot  = ([wl_filt,wl_spec])[i_sort]
        fl_tot  = ([fl_filt,fl_spec])[i_sort]
        tr_tot  = ([tr_filt,tr_spec])[i_sort]
        int_sp(j) = trapz(wl_tot,wl_tot*fl_tot*tr_tot) / $
          trapz(wl_tot,wl_tot*tr_tot)
    endif
endfor

weight   = 1./(eflux*eflux)
scale_b  = where(int_sp ne 0)
scale    = TOTAL(weight(scale_b) * int_sp(scale_b) * flux(scale_b)) / $
  TOTAL(weight(scale_b) * int_sp(scale_b) * int_sp(scale_b))
flux_as  = [REFORM(flux(phot)),REFORM(scale*flux(spec))] 
eflux_as = [REFORM(eflux(phot)),REFORM(scale*eflux(spec))]

rep = where(flux eq -99,n_rep)
if n_rep ge 1 then flux_as(rep) = -99
fit_bands = where(flux ne -99 and int_sp eq 0)

as = CREATE_STRUCT("flux",flux_as,"eflux",eflux_as)

return,as

end
