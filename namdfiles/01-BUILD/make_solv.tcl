set sol_infile solkcsa.pdb
set kcsa_inbase kcsav
set outbase kcsav_solv_raw
package require psfgen
resetpsf
topology ../top_all27_prot_lipid.rtf
segment SOLV {
auto none
first NONE
last NONE
pdb $sol_infile
}
coordpdb $sol_infile SOLV
readpsf ${kcsa_inbase}.psf
coordpdb ${kcsa_inbase}.pdb
writepdb ${outbase}.pdb
writepsf ${outbase}.psf
exit

