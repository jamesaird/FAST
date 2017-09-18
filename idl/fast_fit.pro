function fast_conf_int,chi,adi,prop,chi_thr
  n_prop   = n_elements(prop)
  if n_prop gt 1 then begin
     chi_dim  = min(min(min(min(chi,dim=adi[0]),dim=adi[1]),$
                            dim=adi[2]),dim=adi[3])
     chi_reb  = REBIN(chi_dim,n_prop*100.)
     prop_reb = REBIN(prop,n_prop*100.)
     in       = where(chi_reb le chi_thr,n_in)
     if n_in gt 0 then begin
        l_prop = min(prop_reb(in))
        h_prop = max(prop_reb(in))
     endif else begin
        l_prop = prop_reb(where(chi_reb eq min(chi_reb)))
        h_prop = prop_reb(where(chi_reb eq min(chi_reb)))
     endelse
  endif else begin
     l_prop   = prop
     h_prop   = prop
  endelse
  return,[l_prop,h_prop]
end


function trapz1, x, y
  n = (size(x))[1]
  return, total((x[1:n-1] - x[0:n-2])*(y[1:n-1]+y[0:n-2]),1)/2.
end


FUNCTION fast_scale,flux,eflux,model,det,scale=scale,x_err=x_err,id=id
;...Calculate scale and chi^2 (loop is quicker)
n_dat  = n_elements(det)
chi    = 0.
wmm    = 0.
wfm    = 0.
n_z    = (SIZE(model))[1]
n_tau  = (SIZE(model))[2]
n_metal= (SIZE(model))[3]
n_age  = (SIZE(model))[4]
n_Av   = (SIZE(model))[5]
if KEYWORD_SET(x_err) then begin

   n_fl   = n_elements(eflux)
   eflux2 = sqrt( (REBIN(eflux^2,n_fl,n_z)) + $
                  (x_err*(REBIN(flux,n_elements(flux),n_z)))^2 )

   ;weight=1./(eflux2*eflux2)
   weight = TRANSPOSE(REFORM(REBIN(1./(eflux2*eflux2),n_fl,n_z,n_tau,n_metal,$
                                   n_age,n_Av),n_fl,n_z,n_tau,n_metal,n_age,$
                                   n_Av),[1,2,3,4,5,0])
   wmm    = TOTAL((weight*model*model)[*,*,*,*,*,det],6)
   for k=0,n_dat-1 do begin
      wfm = TEMPORARY(wfm) + weight(*,*,*,*,*,det(k)) $
            * TOTAL(flux(det(k))) * model[*,*,*,*,*,det(k)]
   endfor
   scale  = wfm / wmm
   

   wfm    = [0] & wmm = [0] & weight = [0]

   ;; - redundant with weight
   eflux2 = TRANSPOSE(REFORM(REBIN(eflux2,n_fl,n_z,n_tau,n_metal,n_age,n_Av),$
                             n_fl,n_z,n_tau,n_metal,n_age,n_Av),[1,2,3,4,5,0])

   for k=0,n_dat-1 do begin
      tmp_chi = (TOTAL(flux(det(k))) - scale * model[*,*,*,*,*,det(k)] )/ $
                eflux2(*,*,*,*,*,det(k))
      chi     = TEMPORARY(chi)+(tmp_chi*tmp_chi);*weight(*,*,*,*,*,det(k)))
      tmp_chi = [0]
   endfor
   
   weight=[0]
   eflux2 = [0]
   
endif else begin
   
   weight = REFORM(1. / (eflux*eflux))
   for k=0,n_dat-1 do wmm = TEMPORARY(wmm) + TOTAL(weight(det(k))) * $
      model[*,*,*,*,*,det(k)] * model[*,*,*,*,*,det(k)]
   for k=0,n_dat-1 do wfm = TEMPORARY(wfm) + TOTAL(weight(det(k))) * $
      TOTAL(flux(det(k))) * model[*,*,*,*,*,det(k)]
   scale = wfm / wmm 
   wfm   = [0] & wmm = [0] 
   
   for k=0,n_dat-1 do begin
      tmp_chi = TOTAL(flux(det(k))/eflux(det(k))) - $
                1./TOTAL(eflux(det(k))) * scale * model[*,*,*,*,*,det(k)]
      chi     = TEMPORARY(chi)+tmp_chi*tmp_chi
      tmp_chi = [0]
   endfor
   
endelse
return,chi

END


FUNCTION fast_2scale,flux,eflux,model1,model2,det,scale=scale,x_err=x_err,id=id,agnscale=scale_agn
;...Calculate scale and chi^2 for two grids of model components
n_dat  = n_elements(det)
chi    = 0.
wmm1    = 0.
wfm    = 0.

n_z    = (SIZE(model1))[1]
n_tau  = (SIZE(model1))[2]
n_metal= (SIZE(model1))[3]
n_age  = (SIZE(model1))[4]
n_Av   = (SIZE(model1))[5]

n_agn  = (size(model2))[2]

scale     = fltarr(n_z,n_tau,n_metal,n_age,n_av,n_agn+1) ; normalization of galaxy template
scale_agn = fltarr(n_z,n_tau,n_metal,n_age,n_av,n_agn+1) ; normalization of AGN template
chi  = fltarr(n_z,n_tau,n_metal,n_age,n_av,n_agn+1)      ; +1 as include AGN=0, pure galaxy option

n_fl   = n_elements(eflux)
if KEYWORD_SET(x_err) then begin
   eflux2 = sqrt((REBIN(eflux,n_fl,n_z))^2 + $
                 (x_err*(REBIN(flux,n_elements(flux),n_z)))^2 )
   weight = TRANSPOSE(REFORM(REBIN(1./(eflux2*eflux2),n_fl,n_z,n_tau,n_metal,n_age,n_Av),$
                             n_fl,n_z,n_tau,n_metal,n_age,n_Av),[1,2,3,4,5,0])
endif else begin
   weight = transpose(REFORM(REBIN(1. / (eflux*eflux),n_fl,n_z,n_tau,n_metal,n_age,n_Av),$
                             n_fl,n_z,n_tau,n_metal,n_age,n_Av),[1,2,3,4,5,0])
endelse

fl_rebin = transpose(rebin(flux,n_fl,n_z,n_tau,n_metal,n_age,n_Av),[1,2,3,4,5,0])

for iz=0,n_z-1 do begin         ; should only have one redshift.
   ;; including pure galaxy best fit
   wmm = 0.0D & wmf=0.0D
   for k=0,n_dat-1 do begin
      wmm += weight[iz,*,*,*,*,det[k]] * model1[iz,*,*,*,*,det[k]] * model1[iz,*,*,*,*,det[k]]
      wmf += weight[iz,*,*,*,*,det[k]] * model1[iz,*,*,*,*,det[k]] * fl_rebin[iz,*,*,*,*,det[k]]
   endfor
   a = wmf / wmm
   scale[iz,*,*,*,*,0] = a
   chi_i=0.0D
   for k=0,n_dat-1 do begin
      chi_i += weight[iz,*,*,*,*,det[k]] * ( a*model1[iz,*,*,*,*,det[k]] - flux[det[k]] )^2
   endfor
   chi[iz,*,*,*,*,0] = chi_i

   
   for iagn=1,n_agn do begin
      wml=0.0D
      wlf=0.0D
      wll=0.0D
      for k=0,n_dat-1 do begin
         wml += weight[iz,*,*,*,*,det[k]] * model1[iz,*,*,*,*,det[k]] * model2[iz,iagn-1,det[k]]
         wlf += weight[iz,*,*,*,*,det[k]] * model2[iz,iagn-1,det[k]] * fl_rebin[iz,*,*,*,*,det[k]]
         wll += weight[iz,*,*,*,*,det[k]] * model2[iz,iagn-1,det[k]] * model2[iz,iagn-1,det[k]]
      endfor
      
      a = ( (wmf/wml) - (wlf/wll) )/( (wmm/wml) - (wml/wll) )
      b = (wlf- a*wml)/wll
      ar = rebin(a, n_z,n_tau,n_metal,n_age,n_Av, n_fl)
      br = rebin(b, n_z,n_tau,n_metal,n_age,n_Av, n_fl)
      chi_i=0.0D
      for k=0,n_dat-1 do begin
         chi_i += weight[iz,*,*,*,*,det[k]] * ( a*model1[iz,*,*,*,*,det[k]] + b*model2[iz,iagn-1,det[k]] - $
                                                         flux[det[k]] )^2
      endfor
      chi[iz,*,*,*,*,iagn] = chi_i

      scale[iz,*,*,*,*,iagn] = a
      scale_agn[iz,*,*,*,*,iagn] = b
   endfor

endfor
return,chi

END



PRO fast_fit,id,flux,eflux,zspec,zphot,model,mass_model,sfr_model,z,$
             log_tau,metal,log_age,A_v,name_out,filters,key,tmp_name,$
             AB_ZEROPOINT=AB_ZEROPOINT,$
             SFR100_GRID=sfr100_model,CALCMEAN=CALCMEAN,$
             DUSTPRIOR=DUSTPRIOR,P_AV=p_av,$
             C_INTERVAL=C_INTERVAL,N_SIM=N_SIM,AUTO_SCALE=AUTO_SCALE,$
             TEMP_ERR_FILE=TEMP_ERR_FILE,SAVE_CHI_GRID=SAVE_CHI_GRID,$
             BEST_FIT=BEST_FIT,scale_bands=scale_bands,OUTPUT_DIR=$
             OUTPUT_DIR,lambda=lambda,$
             FITAGN=FITAGN,AGNmodel=AGNmodel,$
             galfl1m=galfl1m,galf2800=galf2800,agnf2800=agnf2800,galf5000=galf5000,agnf5000=agnf5000
  if N_SIM ne 0 then n_int = N_ELEMENTS(C_INTERVAL) else n_int=0
  format  = '(i7,'+strtrim(1+2*n_int,1)+'(f10.4),'+strtrim(1+2*n_int,1)+$
            '(f10.2),'+strtrim(1+2*n_int,1)+'(f10.4),'+strtrim(2+4*n_int,1)+$
            '(f10.2),'+strtrim(4+8*n_int,1)+'(f10.2)'
  noutcol=9+n_int*9*2+1
  if n_elements(sfr100_model) gt 0 then begin
     format=format+',(f15.2)'
     noutcol++
  endif
  if keyword_set(calcmean) then begin
     format=format+',4(f15.2)'
     noutcol+=4
     if n_elements(sfr100_model) gt 0 then begin
        format=format+',2(f15.2)'
        noutcol+=2
     endif
  endif

  format=format+',(e15.2)'      ; L2800
  noutcol++
  if keyword_set(calcmean) then begin
     format=format+',(e15.2)'
     noutcol++
  endif

  if keyword_set(fitagn) then begin
     format=format+',(i15),(f15.2),(f15.2),(f15.2)' ; best template, fagn1m, fagn5000,fagn2800
     noutcol+=4
     if keyword_set(calcmean) then begin
        format=format+',(f15.2),(f15.2),(f15.2)' ; fagn1m_mean, fagn5000_mean, fagn2800_mean
        noutcol+=3
     endif
  endif
  
  format=format+',(e10.2))'
  
  format2 = '(i7,'+strtrim(5+10*n_int,1)+'(i10),'+strtrim(4+8*n_int,1)+$
            '(i10),i10)'
  det     = where(flux ne -99,n_det)
  ;skip fitting if zphot=NAN: if .zout is given and both zph and zsp are 
  ;not defined (this is set in fast_read_observations)
  ;or if all bands have no detection
  if not FINITE(zphot(0), /NAN) and n_det gt 0 then begin
     
     if not KEYWORD_SET(N_SIM) then N_SIM=0 
     if N_SIM ne 0 and not KEYWORD_SET(C_INTERVAL) then C_INTERVAL=68
                              
     ;...Exclude bands that have no coverage (don't forget filters!)
     n_tau   = (SIZE(model))[2]
     n_metal = (SIZE(model))[3]
     n_age   = (SIZE(model))[4]
     n_Av    = (SIZE(model))[5]
     n_dat   = (SIZE(model))[6]
     ;...Reduce model grid to zspec or zphot if N_SIM=0
     if zspec ne -1 or (zphot(0) ge 0 and N_SIM eq 0) then begin
        ;;*** DON'T COPY FULL MODEL GRID IN MEMORY -- SLOW!! INSTEAD rename the slice
        ;cp_model = model & model = [0] &
        cp_z = z
        if zspec ne -1 then zslice = zspec
        if zspec eq -1 then zslice = zphot(0)
        ind_z = where(abs(z-TOTAL(zslice)) eq min(abs(z-TOTAL(zslice))),$
                      n_ind_z)
        ind_z = ROUND(TOTAL(ind_z)/n_ind_z)
        ;model = cp_model(ind_z,*,*,*,*,*)
        z     = cp_z(ind_z)
        modeli = model[ind_z,*,*,*,*,*]

        ;*** AGN fit - reduce AGN model grid
        if keyword_set(fitagn) then begin
           AGNmodeli = AGNmodel[ind_z,*,*]
        endif
     endif
     ;...Reduce model grid between lzphot and hzphot for zphot & N_SIM!=0
     ;   if two confidence intervals are given, take outer!
     if zspec eq -1 and (zphot(0) ge 0 and N_SIM ne 0) then begin
        cp_model = model & model = [0] & cp_z = z
        good_z   = where(z ge zphot(1+2*(n_int-1)) and z le $
                         zphot(2+2*(n_int-1)),tmp_n_z)
        if tmp_n_z eq 0 then good_z = $
           where(abs(z-TOTAL(zphot(0))) eq min(abs(z-TOTAL(zphot(0)))))
        model    = cp_model[good_z,*,*,*,*,*]
        z        = cp_z(good_z)
     endif
     n_z = n_elements(z)

     ;...Calculate extra error from template error function
     if KEYWORD_SET(TEMP_ERR_FILE) then if FILE_TEST(TEMP_ERR_FILE) then $
        begin
        f_err  = FAST_READ(TEMP_ERR_FILE,'flt',comment='#')
        fl_err = f_err[1,*]
        phot   = filters(0,UNIQ(filters(0,where(filters[3,*] eq 1))))
        n_phot = n_elements(phot)
        x_err  = fltarr(n_dat,n_z)
        
        for i=0,n_z-1 do begin
           wl_err = f_err[0,*] * (1.+z(i)) 
           for j=0,n_phot-1 do begin
              filt_ind = where(filters(0,*) eq phot(j))
              good_tr  = where(filters(2,filt_ind) gt 0.005)
              wl_filt  = REFORM(filters(1,filt_ind))
              tr_filt  = REFORM(filters(2,filt_ind))
              good_err = where(wl_err ge min(wl_filt) and $
                               wl_err le max(wl_filt))
              good_err = [min(good_err)-1,good_err,max(good_err)+1]
              wl_err2  = wl_err(good_err)
              fl_err2  = fl_err(good_err)
              tr_err   = INTERPOL(tr_filt,wl_filt,wl_err2) > 0
              fl_filt  = INTERPOL(fl_err,wl_err,wl_filt) > 0
              wl_tot   = [wl_filt,wl_err2]
              i_sort   = SORT(wl_tot)
              wl_tot   = TEMPORARY(wl_tot(i_sort))
              tr_tot   = ([tr_filt,tr_err])[i_sort]
              fl_tot   = ([fl_filt,fl_err2])[i_sort]
              x_err(j,i) = trapz1(wl_tot,tr_tot*wl_tot*fl_tot) / $
                           trapz1(wl_tot,tr_tot*wl_tot)
           endfor
        endfor
     endif
     ;;...Fit real observations
     if keyword_set(FITAGN) then begin
        chi = FAST_2SCALE(flux,eflux,modeli,AGNmodeli,det,scale=galscale,x_err=x_err,id=id,agnscale=agnscale)
        bad = where(galscale lt 0 or agnscale lt 0,nbad) ; remove solutions with negative scalings
        if nbad gt 0 then chi[bad]=!values.f_nan
        n_agn  = (size(AGNmodel))[2]
        scale=galscale
        fagn1m = agnscale/(agnscale + rebin(galfl1m[ind_z,*,*,*,*],n_z,n_tau,n_metal,n_age,n_av,n_agn+1)*galscale)
        agnf2800_r = transpose(rebin([0.0,reform(agnf2800[ind_z,*],n_agn)],n_agn+1,n_z,n_tau,n_metal,n_age,n_av),[1,2,3,4,5,0])
        fagn2800 = agnscale*agnf2800_r/(agnscale*agnf2800_r + rebin(galf2800[ind_z,*,*,*,*],n_z,n_tau,n_metal,n_age,n_av,n_agn+1)*galscale)
        agnf5000_r = transpose(rebin([0.0,reform(agnf5000[ind_z,*],n_agn)],n_agn+1,n_z,n_tau,n_metal,n_age,n_av),[1,2,3,4,5,0])
        fagn5000 = agnscale*agnf5000_r/(agnscale*agnf5000_r + rebin(galf5000[ind_z,*,*,*,*],n_z,n_tau,n_metal,n_age,n_av,n_agn+1)*galscale)
     endif else begin
        if KEYWORD_SET(AUTO_SCALE) then begin
           as  = FAST_AUTO_SCALE(flux,eflux,filters,scale_bands=scale_bands,$
                                 fit_bands=fit_bands,fscale=spec_scale)
           chi = FAST_SCALE(as.flux,as.eflux,modeli,fit_bands,$
                            scale=scale,x_err=x_err)
        endif else begin
           chi = FAST_SCALE(flux,eflux,modeli,det,scale=scale,x_err=x_err,id=id)
        endelse
        chi  = REFORM(chi,n_z,n_tau,n_metal,n_age,n_Av)
     endelse
     i_zb = where(abs(z-zphot(0)) eq min(abs(z-zphot(0))),n_i_zb)
     i_zb = ROUND(TOTAL(TEMPORARY(i_zb))/n_i_zb)


     ;...Make full mass and sfr grid
     if keyword_set(fitagn) then begin
        n_all  = [n_tau,n_metal,n_age,n_Av,n_agn+1,n_z]
        n_all2 = [n_z,n_tau,n_metal,n_age,n_Av,n_agn+1]
        transp = [5,0,1,2,3,4]
     endif else begin
        n_all  = [n_tau,n_metal,n_age,n_Av,n_z]
        n_all2 = [n_z,n_tau,n_metal,n_age,n_Av]
        transp = [4,0,1,2,3]
     endelse
     mass   = REFORM(TRANSPOSE(REFORM(REBIN(mass_model,n_all),n_all),$
                               transp) * scale,n_all2)
     sfr    = REFORM(TRANSPOSE(REFORM(REBIN(sfr_model,n_all),n_all),$
                               transp) * scale,n_all2)


     f2800= rebin(galf2800[ind_z,*,*,*,*],n_all2)*scale
     nfn2800=fast_flam2fnu(2800.,f2800) * 299792458d2/(2800.*1d-8)
     L2800 = nfn2800 / fast_lum2flux(z)

     
     if n_elements(sfr100_model) gt 0 then $
        sfr100    = REFORM(TRANSPOSE(REFORM(REBIN(sfr100_model,n_all),n_all),$
                                     [4,0,1,2,3]) * scale,n_all2)    
     
                                
     ;...Derive best-fit values, can be fast min.val=min(array,i_min)
     if keyword_set(fitagn) then begin
        mchi=min(chi,wsol,/nan)
        ind=array_indices(chi,wsol)
        agn_sol = ind[5]
        tmp_sol = where(chi[*,*,*,*,*,agn_sol] eq mchi,n_sol)
        scale = galscale
        ;; storing best AGN template (0=pure galaxy)
        b_agntempl = agn_sol
        b_fagn1m   = fagn1m[wsol]
        b_fagn2800 = fagn2800[wsol]
        b_fagn5000 = fagn5000[wsol]
     endif else begin
        if zphot(0) ne -1 and zspec eq -1 and N_SIM ne 0 then begin
           mchi=min(chi,wsol,/nan)
           tmp_sol = where(chi(i_zb,*,*,*,*) eq mchi,n_sol)
        endif else begin
           mchi=min(chi,wsol,/nan)
           tmp_sol = where(chi eq mchi,n_sol)
        endelse
     endelse


           
        
     
     if n_sol eq 0 then begin
        print,"    WARNING: no solution for object "+strtrim(id,1)
        print,"             check whether all photometric erros are non-zero"
        printf,1,format=format,id,REPLICATE(-1.0,noutcol)
     endif else begin
        
        tmp_sol  = ROUND(TOTAL(TEMPORARY(tmp_sol))/n_sol)    
        if zphot(0) ne -1 and zspec eq -1 and N_SIM ne 0 then begin
           b_val    = array_indices(chi(i_zb,*,*,*,*),tmp_sol)
           b_val(0) = i_zb
           b_z      = zphot(0)
           min_chi  = min(chi(i_zb,*,*,*,*),/NAN)
        endif else begin
           b_val    = array_indices(chi,tmp_sol)
           b_z      = z(b_val(0))
           min_chi  = min(chi,/NAN)
        endelse
        if keyword_set(fitagn) then b_val[-1]=agn_sol
        b_ltau  = log_tau(b_val(1))
        b_metal = metal(b_val(2))
        b_lage  = log_age(b_val(3))
        b_Av    = A_v(b_val(4))
        b_mass  = mass[wsol]
        ;b_mass  = mass(b_val(0),b_val(1),b_val(2),b_val(3),b_val(4))
        b_lmass = ALOG10(b_mass)
        b_sfr = sfr[wsol]
        ;b_sfr   = sfr(b_val(0),b_val(1),b_val(2),b_val(3),b_val(4))
        if n_elements(sfr100) gt 0 then $
           b_sfr100= sfr100[wsol]
        b_ssfr  = b_sfr / b_mass
        b_efold = b_lage - b_ltau
        if b_sfr eq 0 then b_lsfr = -99 else b_lsfr = Alog10(b_sfr)
        if b_ssfr eq 0 then b_lssfr = -99 else b_lssfr = Alog10(b_ssfr)
        if n_elements(sfr100) gt 0 then $
           if b_sfr100 eq 0 then b_lsfr100 = -99 else b_lsfr100 = Alog10(b_sfr100)
        best    = [b_z,b_ltau,b_metal,b_lage,b_Av,b_lmass,b_lsfr,b_lssfr,$
                   b_efold]
        b_scale = scale[wsol]
        scale   = [0]

        b_L2800 = L2800[wsol]

        ;...fit Monte Carlo simulations
        if N_SIM gt 0 then begin
           sol_sim = fltarr(N_SIM,5) ;z,tau,metal,age,Av
           chi_sim = fltarr(N_SIM)
           fl_sim  = fltarr(n_dat)
           
           ;...Reduce model grid when zphot is know
           if zphot(0) ne -1 and zspec eq -1 then model = $
              TEMPORARY(model(i_zb,*,*,*,*,*))
           
           for i=0,N_SIM-1 do begin

               ;...Make and fit monte carlo simulations
               ;   If zphot are given, MC are all performed at best z_phot
              if n_elements(x_err) ne 0 then begin
                 exflux = sqrt(eflux^2+(x_err(*,i_zb)*flux)^2) 
              endif else exflux = eflux
              for j=0,n_dat-1 do fl_sim(j) = flux(j)+randomn(seed)*exflux(j)
              
              if KEYWORD_SET(AUTO_SCALE) then begin
                 as      = FAST_AUTO_SCALE(fl_sim,eflux,filters)
                 tmp_chi = FAST_SCALE(as.flux,as.eflux,model,fit_bands,$
                                      x_err=x_err)
              endif else begin
                 tmp_chi = FAST_SCALE(fl_sim,eflux,model,det,x_err=x_err)
              endelse
              tmp_min = where(tmp_chi eq min(tmp_chi,/NAN),n_min)
              if n_min eq 0 then begin
                 chi_sim(i) = 1.e9 
              endif else begin
                 sol_sim(i,*) = array_indices(tmp_chi,tmp_min[0])
                 if zphot(0) ne -1 then sol_sim(i,0) = i_zb
                 chi_sim(i) = chi(sol_sim(i,0),sol_sim(i,1),sol_sim(i,2),$
                                  sol_sim(i,3),sol_sim(i,4))
              endelse
           endfor
           
           ;...Derive confidence intervals
           low    = fltarr(n_int,9)
           high   = fltarr(n_int,9)
           efold  = TRANSPOSE(REBIN(log_age,n_age,n_tau)) - $
                    REBIN(log_tau,n_tau,n_age)
           i_sort  = SORT(chi_sim)
           chi_sim = chi_sim(i_sort)
           chi_thr = fltarr(n_int)
           
           for k=0,n_int-1 do begin
                
              chi_thr(k) = INTERPOL(chi_sim,findgen(n_sim),$
                                    C_INTERVAL(k)/100.*n_sim-1.)

              ;...Reduce grid for zphot and two or more confidence intervals
              if zphot(0) ne -1 and zspec eq -1 and k lt n_int-1 then begin
                 good_z  = where(z ge zphot(1+2*(n_int-1)) and z $
                                 le zphot(2+2*(n_int-1)),n_z)
                 cp_chi  = chi
                 cp_mass = mass
                 cp_sfr  = sfr
                 chi     = cp_chi(good_z,*,*,*,*)
                 sfr     = cp_sfr(good_z,*,*,*,*)
                 mass    = cp_mass(good_z,*,*,*,*)
              endif
              ssfr = sfr / mass
              
              if zphot(0) ne -1 and zspec eq -1 then begin
                 diff       = min(chi(i_zb,*,*,*,*),/NAN) - min(chi,/NAN)
                 chi_thr(k) = chi_thr(k) - diff
              endif
              

              in_int  = where(chi le chi_thr(k) or chi eq min(chi,/NAN), $
                              n_grid)
              grid_1s = array_indices(chi,in_int)
              if zphot(0) ne -1 and zspec eq -1 then begin
                 l_z = zphot(1+2*k)
                 h_z = zphot(2+2*k)
              endif else begin
                 l_z = min(z(grid_1s(0,*)))
                 h_z = max(z(grid_1s(0,*)))
              endelse

              l_lmass = Alog10(min(mass(grid_1s(0,*),grid_1s(1,*),$
                                        grid_1s(2,*),grid_1s(3,*),$
                                        grid_1s(4,*)))) < b_lmass
              h_lmass = Alog10(max(mass(grid_1s(0,*),grid_1s(1,*),$
                                        grid_1s(2,*),grid_1s(3,*),$
                                        grid_1s(4,*)))) > b_lmass
              l_sfr   = min(sfr(grid_1s(0,*),grid_1s(1,*),grid_1s(2,*),$
                                grid_1s(3,*),grid_1s(4,*))) < b_sfr
              h_sfr   = max(sfr(grid_1s(0,*),grid_1s(1,*),grid_1s(2,*),$
                                grid_1s(3,*),grid_1s(4,*))) > b_sfr
              l_ssfr  = min(ssfr(grid_1s(0,*),grid_1s(1,*),grid_1s(2,*),$
                                 grid_1s(3,*),grid_1s(4,*))) < b_ssfr
              h_ssfr  = max(ssfr(grid_1s(0,*),grid_1s(1,*),grid_1s(2,*),$
                                 grid_1s(3,*),grid_1s(4,*))) > b_ssfr
              l_efold = min(efold(grid_1s(1,*),grid_1s(3,*))) < b_efold
              h_efold = max(efold(grid_1s(1,*),grid_1s(3,*))) > b_efold
              
              ci_av   = fast_conf_int(chi,[1,1,1,1],A_v,chi_thr(k))
              l_Av    = ci_av[0]
              h_Av    = ci_av[1]
              ci_met  = fast_conf_int(chi,[1,1,2,2],metal,chi_thr(k))
              l_metal = ci_met[0]
              h_metal = ci_met[1]
              ci_lage = fast_conf_int(chi,[1,1,1,2],log_age,chi_thr(k))
              l_lage  = ci_lage[0]
              h_lage  = ci_lage[1]
              ci_ltau = fast_conf_int(chi,[1,2,2,2],log_tau,chi_thr(k))
              l_ltau  = ci_ltau[0]
              h_ltau  = ci_ltau[1]

              if l_sfr eq 0 then l_lsfr = -99 else l_lsfr = Alog10(l_sfr)
              if h_sfr eq 0 then h_lsfr = -99 else h_lsfr = Alog10(h_sfr)
              if l_ssfr eq 0 then l_lssfr = -99 else l_lssfr = Alog10(l_ssfr)
              if h_ssfr eq 0 then h_lssfr = -99 else h_lssfr = Alog10(h_ssfr)
              
              low(k,*)  = [l_z,l_ltau,l_metal,l_lage,l_Av,l_lmass,l_lsfr,$
                           l_lssfr,l_efold]
              high(k,*) = [h_z,h_ltau,h_metal,h_lage,h_Av,h_lmass,h_lsfr,$
                           h_lssfr,h_efold]
              
              if zphot(0) ne -1 and zspec eq -1 and k lt n_int-1 then begin
                 chi  = cp_chi  & cp_chi  = [0]
                 mass = cp_mass & cp_mass = [0]
                 sfr  = cp_sfr  & cp_sfr  = [0]
              endif 
              
           endfor
           
        endif else begin
           chi_thr = -1
           low     = -1
           high    = -1
        endelse



        ;; JAA addition - calculate mean+median from chi^2
        if keyword_set(calcmean) then begin
           if keyword_set(fitagn) then begin
              p=exp(-(chi-min(chi,/nan))/2.)
              n_agn=n_elements(agnmodeli[0,*,0])
              sfr = rebin(sfr,n_z,n_tau,n_metal,n_age,n_av,n_agn+1)
              if n_elements(sfr100) gt 0 then sfr100=rebin(sfr100,n_z,n_tau,n_metal,n_age,n_av,n_agn+1)
              mass = rebin(mass,n_z,n_tau,n_metal,n_age,n_av,n_agn+1)
           endif else p=exp(-(chi-min(chi,/nan))/2.)

           if keyword_set(dustprior) then begin
              ;p_av=a_v^2.0*exp(-a_v*2.0/avpeak)
                                ;p_av=transpose(rebin(p_av,n_elements(a_v),1,n_elements(log_tau),$
                                ;                    n_elements(metal),n_elements(log_age)),[1,2,3,4,0])
              p = p*p_av
              p = p/max(p,/nan)
           endif

           ;; mean
           tot_p=total(p,/nan)>1.
           sfr_mean = total(sfr*p,/nan)/tot_p
           if sfr_mean eq 0 then lsfr_mean=-99. else lsfr_mean=alog10(sfr_mean)
           if n_elements(sfr100) gt 0 then begin
              sfr100_mean = total(sfr100*p,/nan)/tot_p
              if sfr100_mean eq 0 then lsfr100_mean=-99. else lsfr100_mean=alog10(sfr100_mean)
           endif
           mass_mean = total(mass*p,/nan)/tot_p
           if mass_mean eq 0 then lmass_mean=-99. else lmass_mean=alog10(mass_mean)

           L2800_mean = total(L2800*p,/nan)/tot_p

           if keyword_set(fitagn) then begin
              fagn1m_mean = total(fagn1m*p,/nan)/tot_p
              fagn2800_mean = total(fagn2800*p,/nan)/tot_p
              fagn5000_mean = total(fagn5000*p,/nan)/tot_p
           endif

           ;; median
           s=sort(sfr)
           pcum=total(p[s],/cum,/nan)
           sfr_med = interpol(sfr[s],pcum/max(pcum),0.5)
           if sfr_med eq 0 then lsfr_med=-99. else lsfr_med=alog10(sfr_med)

           if n_elements(sfr100) gt 0 then begin
              s=sort(sfr100)
              pcum=total(p[s],/cum,/nan)
              sfr100_med = interpol(sfr100[s],pcum/max(pcum),0.5)
              if sfr100_med eq 0 then lsfr100_med=-99. else lsfr100_med=alog10(sfr100_med)
           endif

           s=sort(mass)
           pcum=total(p[s],/cum,/nan)
           mass_med = interpol(mass[s],pcum/max(pcum),0.5)
           if mass_med eq 0 then lmass_med=-99. else lmass_med=alog10(mass_med)
        endif

        
        
        ;...Print values
        CASE n_int OF
           0: all = best
           1: all = REFORM(TRANSPOSE([[best],[REFORM(low)],[REFORM(high)]]),$
                           9*(2*n_int+1))
           2: all = REFORM(TRANSPOSE([[best],[REFORM(low(0,*))],$
                                      [REFORM(high(0,*))],[REFORM(low(1,*))],$
                                      [REFORM(high(1,*))]]),9*(2*n_int+1))
        ENDCASE

        if n_elements(sfr100) gt 0 then all=[all,b_lsfr100]
        if keyword_set(calcmean) then begin
           all=[all,lmass_mean,lmass_med,lsfr_mean,lsfr_med]
           if n_elements(sfr100) gt 0 then all=[all,lsfr100_mean,lsfr100_med]
        endif

        all=[all,b_L2800]
        if keyword_set(calcmean) then all=[all,L2800_mean]

        if keyword_set(fitagn) then begin
           all=[all,b_agntempl,b_fagn1m,b_fagn5000,b_fagn2800]
           if keyword_set(calcmean) then all=[all,fagn1m_mean,fagn5000_mean,fagn2800_mean]
        endif


        
        deg_fix = 0.
        if n_Av gt 1 then deg_fix = deg_fix+1.
        if n_age gt 1 then deg_fix = deg_fix+1.
        if n_metal gt 1 then deg_fix = deg_fix+1.
        if n_tau gt 1 then deg_fix = deg_fix+1.
        if n_z gt 1 then deg_fix = deg_fix+1

        n_degree = n_det - deg_fix 

        printf,1,format=format,id,all,min_chi / n_degree
        
        if keyword_set(SAVE_CHI_GRID) then begin
           if n_elements(cp_z) eq 0 then i_best = b_val else $
              i_best = [where(cp_z eq best(0)),b_val(1:4)]
           if keyword_set(fitagn) then begin
              SAVE,chi,n_degree,mass,sfr,z,chi_thr,best,i_best,low,high,b_scale,key,agnscale,fagn1m,$
                   FILENAME=tmp_name+'/chiagn_'+name_out+'.'+strtrim(id,1)+$
                   '.save'
           endif else if n_elements(spec_scale) ne 0 then begin
              SAVE,chi,n_degree,mass,sfr,z,chi_thr,best,i_best,low,high,b_scale,key,$
                   spec_scale,FILENAME=tmp_name+'/chi_'+name_out+'.'+$
                   strtrim(id,1)+'.save'
           endif else begin
              SAVE,chi,n_degree,mass,sfr,z,chi_thr,best,i_best,low,high,b_scale,key,$
                   FILENAME=tmp_name+'/chi_'+name_out+'.'+strtrim(id,1)+$
                   '.save'
           endelse
        endif
        
        ;best fit at the input resolution
        best_lfit = b_scale*reform(modeli(where(z eq best(0)),where(log_tau eq best(1)),$
                                         where(metal eq best(2)),where(log_age eq best(3)),$
                                         where(A_V eq best(4)),*))

        if KEYWORD_SET(BEST_FIT) then begin
           if keyword_set(fitagn) then begin
              agn_lfit = agnscale[wsol] * agnmodeli[b_val[0],b_agntempl-1,*]
              FAST_BEST_FIT, OUTPUT_DIR, name_out, $
                             id, key, best, b_scale, lambda, best_lfit,$
                             AB_ZEROPOINT=AB_ZEROPOINT,$
                             AGN_TEMPL=b_agntempl,agn_scale=agnscale[wsol],agn_lfit=agn_lfit
           endif else begin
              FAST_BEST_FIT, OUTPUT_DIR, name_out, $                             
                             id, key, best, b_scale, lambda, best_lfit,$
                             AB_ZEROPOINT=AB_ZEROPOINT
           endelse
        endif
        
     endelse
     
     
     if zspec ne -1 or zphot(0) ge 0 then begin
        ;;** DON'T DO THIS - takes most of the time
        ;model    = cp_model
        cp_model = [0]
        z        = cp_z
     endif
     
  endif else begin
     
     printf,1,format=format,id,REPLICATE(-1.0,noutcol)
     
  endelse
END
