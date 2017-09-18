FUNCTION fast_z_arr,Z_MIN,Z_MAX,Z_STEP,Z_STEP_TYPE,CATALOG,$
                    NAME_ZPHOT=NAME_ZPHOT,N_SIM=N_SIM

  if Z_MAX eq Z_MIN then z = Z_MIN else begin
     if Z_STEP_TYPE eq 0 then begin
        n_z   = round((Z_MAX - Z_MIN) / Z_STEP) + 1
        z     = Z_MIN + Z_STEP * findgen(n_z)
     endif else begin
        z     = Z_MIN
        while max(z) lt Z_MAX do z = [z, max(z)+Z_STEP*(1+max(z))]
     endelse
  endelse
  
  if not KEYWORD_SET(SPECTRUM) then begin
     if FILE_TEST(CATALOG+'.zout') and N_SIM eq 0 then begin
        zcat = FAST_READ(CATALOG+'.zout','flt',comment='#',header=header)
        n_ob = (size(zcat))[2]
        if n_ob lt n_elements(z) then begin
           
           ;take all photometric redshift
           col_zph = where(header eq NAME_ZPHOT,n_head) ;add z_spec
           if n_head eq 1 then begin
              zphot   = REFORM(zcat(col_zph,*))
              z_out    = zphot
           endif else z_out = replicate(0,n_ob)

          ;replace galaxies for which spec zs are available from both catalogs
          ;ignore zspec in [].zout if available in [].cat
           pcat    = FAST_READ(CATALOG+'.cat','flt',comment='#',header=header2)
           col_zs2 = where(header2 eq 'z_spec',n_head2)
           if n_head2 eq 1 then begin
              zspec2  = REFORM(pcat(col_zs2,*))
              rep_z2  = where(zspec2 gt 0,n_rep2)
              if n_rep2 ge 1 then z_out(rep_z2) = zspec2(rep_z2)
           endif else begin
              col_zs1 = where(header eq 'z_spec',n_head)
              if n_head eq 1 then begin
                 zspec1  = REFORM(zcat(col_zs1,*))
                 rep_z1  = where(zspec1 gt 0,n_rep1)
                 if n_rep1 ge 1 then z_out(rep_z1) = zspec1(rep_z1)
              endif
           endelse
           
           z       = z_out(UNIQ(z_out, SORT(z_out))) > 0.00001
        endif
     endif
  endif  

  RETURN,z

END
