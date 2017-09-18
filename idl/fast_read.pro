;......................................................................
;description:
; Reads unformatted data file and header
;
;parameters:
; input    input_file  name input file
;          format      'flt','str','dbl' or 'int'
; keyword  comment
;          header
; return   data        data array
;......................................................................

FUNCTION fast_read, input_file, format, comment=comment, header=header

header = ''
if not KEYWORD_SET(comment) then comment='#'

CASE format OF
    'str': type = 7
    'int': type = 3 ; long 
    'flt': type = 4 
    'dbl': type = 5
    ELSE: type = 4
ENDCASE

; read entire file in string buffer
n = file_lines(input_file)
lines = strarr(n)
openr, lun, input_file, /get_lun
readf, lun, lines
close, lun & free_lun, lun

; locate comment and data lines from first character
first_char = strmid(lines,0,1)
ihdr = where(first_char eq comment, nhdr, complement=idata, ncomplement=ndata)
if ndata eq 0 then message, 'ERROR: '+fname+' contains no data'
ncol = n_elements(strsplit(lines[idata[0]]))

; parse header
if nhdr gt 0 then begin
  h = strmid(lines[ihdr],1)   ; strip comment symbol
  for i=0,n_elements(h)-1 do begin
      head_line = strsplit(h[i],/extract)
      nhcol     = n_elements(head_line)
      if nhcol eq ncol then begin
          if n_elements(header) eq 0 or n_elements(header) eq 1 then $
            header = head_line else header = [[header],[head_line]]
      endif
  endfor
endif

; split lines into words, cast to correct type, and fill array
; 45% of the total time is spent in strsplit, 45% in casting of 
; the data type (fix)
data = make_array(ncol, ndata, type=type)
for i=0L, ndata-1 do begin
    if n_elements(strsplit(lines[idata[i]],/extract)) ne ncol then begin
       print,'ERROR: when reading '+input_file
       print,'       number of columns is not constant throughout file'
       exit
    endif
    data[*,i] = fix( strsplit(lines[idata[i]],/extract), type=type)
endfor
 
return, data

END
