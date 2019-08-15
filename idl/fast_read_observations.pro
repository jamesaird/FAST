FUNCTION fast_read_observations,lambda,CATALOG=CATALOG,SPECTRUM=SPECTRUM,$
                                AB_ZEROPOINT=AB_ZEROPOINT,C_INTERVAL=$
                                C_INTERVAL,NAME_ZPHOT=NAME_ZPHOT,$
                                H_CAT=H_CAT


if not KEYWORD_SET(NAME_ZPHOT) then name_zphot = 'z_phot'
if not KEYWORD_SET(AB_ZEROPOINT) then AB_ZEROPOINT = 25.

if KEYWORD_SET(CATALOG) then begin

    if not KEYWORD_SET(C_INTERVAL) then C_INTERVAL = 68

    print,"Read photometric catalog: "+CATALOG+".cat"
    cat   = FAST_READ(CATALOG+".cat",'flt',comment='#',header=h_cat)
    n_gal = (SIZE(cat))[2]

    ;...change header file in case FILT_TRANSLATE is provided
    if FILE_TEST(CATALOG+'.translate') then begin
        print,"Read catalog translate file: "+CATALOG+".translate"
        readcol,CATALOG+'.translate',tr1,tr2,format='(A,A)',/silent

        for i=0,n_elements(tr1)-1 do begin
            rep_h = where(h_cat eq tr1(i),n_rep_h)
            if n_rep_h eq 1 then begin
                if (size(h_cat))[0] eq 1 then begin
                    h_cat(rep_h) = tr2(i)
                endif else begin
                    rep_h2 = array_indices(h_cat,rep_h)
                    h_cat(rep_h2(0),rep_h2(1)) = tr2(i)
                endelse
            endif
        endfor
    endif

    ;...read and select columns
    if (size(h_cat))[0] eq 1 then n_line_h = 1
    if (size(h_cat))[0] eq 2 then n_line_h = (size(h_cat))[2]
    h_cat   = REFORM(h_cat,(size(h_cat))[1],n_line_h)
    i_line  = (array_indices(h_cat,where(strmatch(h_cat,'F*[1234567890]') eq 1)))[1,0]
    id_ind  = where(h_cat(*,i_line) eq 'id' or h_cat(*,i_line) eq 'ID',n_id)
    fl_ind  = where(strmatch(h_cat(*,i_line),'F*[1234567890]') eq 1,n_filt)
    err_ind = where(strmatch(h_cat(*,i_line),'E*[1234567890]') eq 1)
    tot_ind = where(strmatch(h_cat(*,i_line),'TOT*[1234567890]') eq 1,n_tot)
    zsp_ind = where(strmatch(h_cat(*,i_line),'z_spec',/fold_case),n_zsp)

    if n_filt eq 0 then begin
       print,'ERROR: no flux column identifications found '+$
             '(F[n], E[n], etc.) in '+CATALOG+'.cat'
       exit
    endif
    if n_id eq 0 then begin
       print,'ERROR: no ID column found in '+CATALOG+'.cat'
       exit
    endif

    ;...correct aperture fluxes for total fluxes
    if n_tot eq 1 then begin
        band_tot = TOTAL(fix(strmid(h_cat(tot_ind,i_line),3)))
        band_fl  = fix(strmid(h_cat(fl_ind,i_line),1))
        band_ap  = fl_ind(where(band_fl eq band_tot))
        corr     = REFORM(cat(tot_ind,*) / cat(band_ap,*) )
        corr     = TRANSPOSE(REBIN(TEMPORARY(corr),n_gal,n_filt),[1,0])
    endif else begin
        corr     = 1.
    endelse
    id    = fix(cat(id_ind,*),type=3)
    efnu  = corr*cat(err_ind,*)
    fnu   = corr*cat(fl_ind,*)
    if n_zsp eq 1 then zspec=cat(zsp_ind,*) else zspec=REPLICATE(-1,n_gal)
    rep_zsp = where(zspec lt 0,n_rep_zsp)
    if n_rep_zsp ge 1 then zspec(rep_zsp) = -1

    ;...convert fnu into flambda
    eflux = 3.e8*1.e10*efnu / REBIN(lambda[0:n_filt-1] * $
                                    lambda[0:n_filt-1],n_filt,n_gal)
    flux  = 3.e8*1.e10*fnu / REBIN(lambda[0:n_filt-1]* $
                                   lambda[0:n_filt-1],n_filt,n_gal)
    rep   = where(cat(fl_ind,*) le -99,n_rep) ;can be 0 ; **J Aird - bug fix - use original catalogue to avoid rounding errors
    if n_rep ge 1 then begin
        i_rep = array_indices(flux,rep)
        eflux(i_rep(0,*),i_rep(1,*)) = -99
        flux(i_rep(0,*),i_rep(1,*)) = -99
     endif
    cat   = [0]
    fnu   = 0
    efnu  = 0

    ;...read photometric redshift file when available
    if FILE_TEST(CATALOG+'.zout') then begin
        print,"Read photometric redshifts: "+CATALOG+".zout"
        zcat = FAST_READ(CATALOG+'.zout','flt',comment='#',$
                         header=h_zcat)
        if (size(zcat))[2] ne n_gal then begin
            print,"ERROR: number of objects in "+CATALOG+".cat and "+$
              CATALOG+".zout differ"
            exit
        endif
        if (size(h_zcat))[0] eq 1 then begin
            zcol  = where(h_zcat eq name_zphot,n_zcol)
            zpcol = where(h_zcat eq 'z_spec',n_zsp2)
            n_int   = N_ELEMENTS(C_INTERVAL)
            int_col = fltarr(2*n_int)
            for j=0,n_int-1 do begin
                int_col(2*j)   = where(h_zcat eq 'l'+STRING(C_INTERVAL(j),$
                                                            f='(i2)'))
                int_col(2*j+1) = where(h_zcat eq 'u'+STRING(C_INTERVAL(j),$
                                                            f='(i2)'))
            endfor
        endif else begin
            i_line  = (array_indices(h_zcat,where(h_zcat eq name_zphot)))[1]
            zcol    = where(h_zcat(*,i_line) eq name_zphot,n_zcol)
            zpcol   = where(h_zcat(*,i_line) eq 'z_spec',n_zsp2)
            n_int   = N_ELEMENTS(C_INTERVAL)
            int_col = fltarr(2*n_int)
            for j=0,n_int-1 do begin
                int_col(2*j)   = where(h_zcat(*,i_line) eq 'l'+$
                                       STRING(C_INTERVAL(j),f='(i2)'))
                int_col(2*j+1) = where(h_zcat(*,i_line) eq 'h'+$
                                       STRING(C_INTERVAL(j),f='(i2)'))
            endfor
         endelse
        if n_zcol eq 0 then begin
           if not KEYWORD_SET(SPECTRUM) then begin
              print,"WARNING: column '"+name_zphot+"' not found in "+CATALOG+".zout"
              print,"         If you want to continue fitting with redshift free type "+$
                    "'.cont'"
              stop
           endif
           zphot = REPLICATE(-1,1+2*n_int,n_gal)
        endif else zphot = zcat([zcol,int_col],*)
        if n_zsp eq 0 and n_zsp2 eq 1 then begin
            zspec    = REFORM(zcat(zpcol,*))
            rep_zsp2 = where(zspec lt 0,n_rep_zsp2)
            if n_rep_zsp2 ge 1 then zspec(rep_zsp2) = -1
        endif
        
        ;...put zphot to NAN if .zout: zspec and zphot are not defined
        rep = where((zphot(0,*) lt 0 or finite(zphot(0,*),/nan)) and $
                    zspec lt 0,n_rep)
        if n_rep ge 1 then zphot(0,rep) = alog10(-1)
        rep = where(zspec gt 0 and finite(zphot(0),/nan),n_rep)
        if n_rep ge 1 then zphot(0,rep) = -1

        zcat  = [0]

    endif else begin
        print,"No photometric redshift input"
        zphot = REPLICATE(-1,1+2*N_ELEMENTS(C_INTERVAL),n_gal)
    endelse

    name   = file_basename(CATALOG)

endif

if KEYWORD_SET(SPECTRUM) then begin

    print,"Read spectrum: "+SPECTRUM+'.spec'

    spec    = FAST_READ(SPECTRUM+'.spec','flt',comment='#',header=h_spec)
    fl_ind  = where(strmid(h_spec,0,1) eq 'F',n_gal2)
    err_ind = where(strmid(h_spec,0,1) eq 'E')
    wl_ind  = where(h_spec eq 'wl')
    bin_ind = where(h_spec eq 'bin')
    tr_ind  = where(h_spec eq 'tr')
    id_spec = strmid(h_spec(fl_ind),1)

    good    = where(REFORM(spec(tr_ind,*)) eq 1)
    st_bin  = UNIQ(spec(bin_ind,*))
    n_bin   = n_elements(st_bin)
    corr    = 10^((AB_ZEROPOINT+48.57)/2.5 - 19.)
    fl_sp   = (spec(fl_ind,*))[*,st_bin] * corr
    efl_sp  = (spec(err_ind,*))[*,st_bin] * corr

    if KEYWORD_SET(CATALOG) then if CATALOG ne '' then begin

        if FILE_TEST(CATALOG+'.zout') then print,"  NOTE: Redshift "+$
          "is fixed to z_spec when available, but not to z_phot"

        match = intarr(n_elements(id_spec))
        for i=0,n_elements(id_spec)-1 do begin
            match(i) = where(id eq id_spec(i),n_match)
            if n_match eq 0 then begin
                print,"ERROR: no photometric data for object",id_spec(i)
                exit
            endif
        endfor

        flux  = [TEMPORARY(flux(*,match)),TRANSPOSE(fl_sp,[1,0])]
        eflux = [TEMPORARY(eflux(*,match)),TRANSPOSE(efl_sp,[1,0])] 
        name  = file_basename(CATALOG)+'.'+file_basename(SPECTRUM)
        zphot = REPLICATE(-1,3,n_gal2) ;assume no zphot 
        zspec = TEMPORARY(zspec(match))
        id    = id_spec

    endif
    if n_elements(id) eq 0 then begin
        id    = id_spec
        flux  = TRANSPOSE(fl_sp,[1,0])
        eflux = TRANSPOSE(efl_sp,[1,0])
        zspec = REPLICATE(-1,n_gal2)
        zphot = REPLICATE(-1,3,n_gal2)
        name  = file_basename(SPECTRUM)
    endif


endif

data = CREATE_STRUCT("id",id,"flux",flux,"eflux",eflux,"zphot",zphot,$
                     "zspec",zspec,"name",name)

RETURN,data

END
