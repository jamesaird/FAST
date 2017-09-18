;+
; NAME:
;   fast_calcsfr100
; PURPOSE:
;   Calculate SFR100 (SFR averaged over 100Myr) for a given input SFR,
;   age, tau, and SFH type
; INPUTS:
;  sfr, log_age,log_tau
; KEYWORDS:
;  /grid - if set assumes sfr is a ntau x nmetal x nage
;          grid, with log_age and log_tau of grid given by
;          inputs. Otherwise dimensions of log_age and log_tau must
;          match dimensions of sfr
;  sfh - select star formation history (del or exp only for now)
; RETURN VALUE:
;   sfr100 - sfr average over 100Myr of galaxy age (if shorter)
; MODIFICATION HISTORY:
;  Created 11 Nov 2015, J. Aird, Cambridge
;-
function fast_calcsfr100,sfr,log_age,log_tau,grid=grid,sfh=sfh

  ;; SFR100 - grid
  if keyword_set(grid) then begin
     dim=size(sfr,/dim)
     sfr100=dblarr(dim)
     t=transpose(rebin(10d^(log_age-6.0),dim[2],dim[1],dim[0]),[2,1,0])
     tau=rebin(10d^(log_tau -6.0),dim[0],dim[1],dim[2])
  endif else begin
     sfr100=sfr*0.0
     t=10d^(log_age-6.0)
     tau=10d^(log_tau-6.0)
  endelse 

  if sfh eq 'del' then begin
     w=where(t gt 100,n,complement=q,ncomp=nq)
     if n gt 0 then begin
        fctw = (tau[w]/t[w]) * (tau[w]/100.0) * $
               ( exp(100d/tau[w])*((t[w]/tau[w])-(100d/tau[w])+1) - 1.0 - (t[w]/tau[w])  )
        sfr100[w]=sfr[w]*fctw
     endif
     if nq gt 0 then begin
        fctq =  (tau[q]/t[q])^2 * ( exp(t[q]/tau[q]) - 1.0  - (t[q]/tau[q]) )
        sfr100[q] = sfr[q]*fctq
     endif
  endif else if sfh eq 'exp' then begin
     w=where(t gt 100,n,complement=q,ncomp=nq)
     if n gt 0 then begin
        fctw = (tau[w]/100.)*(exp(100.0/tau[w]) -1.0)
        sfr100[w]=sfr[w]*fctw
     endif
     if nq gt 0 then begin
        fctq =  (tau[q]/t[q])*(exp(t[q]/tau[q])-1.0)
        sfr100[q] = sfr[q]*fctq
     endif
  endif else message,'Must specify SFH as del or exp to calculate SFR100'

  return,sfr100
end





