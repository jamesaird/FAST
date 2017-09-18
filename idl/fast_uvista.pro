;+
; Wrapper to run FAST for UVISTA. Runs multiple instances and
; concatenates output catalogues
;-
pro fast_uvista,outname=outname,num_threads=num_threads

  if n_elements(outname) eq 0 then outname='UVISTA_final_v4.1.fout'

  if n_elements(num_threads) eq 0 then num_threads=3

  cat=plxm_read_uvista()

  nset=n_elements(cat)/num_threads +1


  
  CD, CURRENT=dir

  
  for iproc=0,num_threads-1 do begin
     if iproc eq 0 then parproc=obj_new('IDL_IDLBridge',output='sub'+strtrim(string(iproc+1),2)+'.log') $
     else parproc=[parproc,obj_new('IDL_IDLBridge',output='sub'+strtrim(string(iproc+1),2)+'.log')]

     parproc[iproc]->execute,'cd,"'+dir+'"'
     command='fast,subset=lindgen('+string(nset)+')+'+string(iproc*nset)+'L,subname="sub'+strtrim(string(iproc+1),2)+'"'
     print,'Executing: ',command
     print,"...using process: ",iproc     
     parproc[iproc]->execute,command,/nowait
     wait,30.0
  endfor

  for iproc=0,num_threads-1 do  begin
     print,'Waiting for process',iproc,' to end..',format='((A20,I2,A10),$)'
     while parproc[iproc]->status() eq 1 do begin
        print,'.',format='(A1,$)'
        wait,2.
        endwhile
     print,'ended'
     if iproc eq 0 then spawn,'cp sub'+strtrim(string(iproc+1),2)+'.fout '+outname $
     else spawn,'grep -v "#" sub'+strtrim(string(iproc+1),2)+'.fout >>'+outname
     
  endfor
  obj_destroy,parproc
  delvarx,parproc

end
