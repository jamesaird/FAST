function fast_convolve,wl_rest,synspec,dispersion

;apply for a certain dz/(1+z)
;apply for a certain resolution of an instrument


;dispersion: in resolution elements (dz * resolution)
n_ages   = (SIZE(synspec))[1]
synspec2 = synspec
for i=0,n_ages-1 do synspec2(i,*) = SMOOTH(synspec(i,*),dispersion)
synspec  = synspec2

RETURN,synspec

end
