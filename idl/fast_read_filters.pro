function fast_no_filters,CATALOG

comment = '#'
l_com   = strlen(comment)
line    = ''
w_head  = 0
w_data  = 0

;...determine the number of data columns
OPENR, unit, CATALOG+'.cat', /GET_LUN
while w_data eq 0 do begin
    READF, unit, line
    col   = strsplit(line, /extract)
    if strmid(col(0),0,l_com) ne comment then begin
        n_col  = n_elements(col)
        w_data = 1
    endif
endwhile
CLOSE, unit
FREE_LUN, unit

;...read header with same number of columns as the data
OPENR, unit, CATALOG+'.cat', /GET_LUN
while w_head eq 0 do begin
    READF, unit, line
    col = strsplit(line, /extract)
    if strmid(col(0),0,l_com) eq comment and n_elements(col) gt 1 $
      then begin
        ;...first entry may be attached to '#'
        if col(0) eq comment then tmp_head = col(1:(n_elements(col)-1)) $
        else tmp_head = [strmid(col(0),1),col(1:(n_elements(col)-1))]
        if n_elements(tmp_head) eq n_col then begin
            if n_elements(header) eq 0 then header = tmp_head else $
              header = [[header],[tmp_head]]
        endif
    endif else begin
        w_head = 1
    endelse
endwhile
CLOSE, unit
FREE_LUN, unit

;...change header file in case CATALOG.TRANSLATE is provided
if FILE_TEST(CATALOG+'.translate') then begin
   readcol,CATALOG+'.translate',tr1,tr2,format='(A,A)',/silent
   for i=0,n_elements(tr1)-1 do begin
      rep_h = where(strcmp(header,tr1(i),/fold_case) eq 1,n_rep_h)
      if n_rep_h eq 1 then begin
         if (size(header))[0] eq 1 then begin
            header(rep_h) = tr2(i)
         endif else begin
            rep_h2 = array_indices(header,rep_h)
            header(rep_h2(0),rep_h2(1)) = tr2(i)
         endelse
      endif
   endfor
endif


;...determine filters from the header info
if (size(header))[0] eq 0 then begin
   print,"ERROR: no header found in "+CATALOG
   print,"       Check that number of data columns is the same as header columns"
   exit
endif
if (size(header))[0] eq 1 then n_line_h = 1
if (size(header))[0] eq 2 then n_line_h = (size(header))[2]
header = REFORM(header,(size(header))[1],n_line_h)
i_line = (array_indices(header,where(strmatch(header,'F*[1234567890]') eq 1)))[1,0]
fl_ind = where(strmatch(header,'F*[1234567890]') eq 1,n_filt)
filt   = strmid(header(fl_ind,i_line),1)

return,filt

end



FUNCTION fast_read_filters,FILTERS_RES,lambda,CATALOG=CATALOG,$
                           SPECTRUM=SPECTRUM,no_filt=no_filt

;...read photometric filters
;   RECENT EDIT: added "or Keyowrd_set(no_filt)
if KEYWORD_SET(CATALOG) or KEYWORD_SET(no_filt) then begin

   if not KEYWORD_SET(no_filt) then no_filt = fast_no_filters(CATALOG)
   n_filt     = n_elements(no_filt)

    if not file_test(FILTERS_RES) then begin
        print,"ERROR: defined filter file is not available" & exit
    endif
    READCOL,FILTERS_RES,no,wl,tr,FORMAT='I,F,F',/silent
    st_filt    = where(no eq 1,n_filt_res)
    en_filt    = [st_filt(1:n_filt_res-1)-1,n_elements(no)-1]
    nel_filt   = en_filt-st_filt+1
    lambda     = fltarr(n_elements(no_filt))
    c_no_filt  = no_filt-1

    for i=0,n_filt-1 do begin
       tmp_filt = [[REPLICATE(i,nel_filt(c_no_filt(i)))],$
                   [wl(st_filt(c_no_filt(i)):en_filt(c_no_filt(i)))],$
                   [tr(st_filt(c_no_filt(i)):en_filt(c_no_filt(i)))],$
                   [REPLICATE(1,nel_filt(c_no_filt(i)))]] ;photo code
       if i eq 0 then begin
          filters = tmp_filt
       endif else begin
          filters = [filters,tmp_filt]
       endelse
       lambda[i] = TOTAL(tmp_filt(*,1)*tmp_filt(*,1)*tmp_filt(*,2)) / $
                   TOTAL(tmp_filt(*,1)*tmp_filt(*,2))
     endfor
    filters = TRANSPOSE(TEMPORARY(filters),[1,0])
endif 

if n_elements(n_filt) eq 0 then n_filt  = 0
    
;...read spectroscopic binning info
if KEYWORD_SET(SPECTRUM) then begin
    spec        = FAST_READ(SPECTRUM+'.spec','flt',comment='#')
    n_spec      = (size(spec))[2]
    n_bin       = n_elements(UNIQ(spec(0,*)))
    bins        = fltarr(4,n_spec)
    bins(0,*)   = spec(0,*)+n_filt
    bins(1:2,*) = spec(1:2,*)
    bins(3,*)   = REPLICATE(0,n_spec) ;spectral code
    if n_elements(filters) eq 0 then filters = bins else $
      filters   = [[temporary(filters)],[bins]]
    tmp_lam     = fltarr(n_bin)
    for i=0l,n_bin-1 do begin
        bin        = where(spec(0,*) eq i)
        tmp_lam(i) = TOTAL(spec(1,bin) * spec(2,bin)) / TOTAL(spec(2,bin))
    endfor
    if n_filt eq 0 then lambda = tmp_lam else $
      lambda = [TEMPORARY(lambda),tmp_lam]
endif

RETURN,filters

END


