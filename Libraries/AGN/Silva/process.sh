# process the Silva templates to convert from microns,nufnu  to Angstrom,flambda
for file in *.sed_orig;
do
    awk '{print $1*10000" "$2*299792458e2/($1*$1)}' $file > "${file%%.sed_orig}.sed"
done


