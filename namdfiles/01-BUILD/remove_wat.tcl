mol load psf kcsav_solv_raw.psf pdb kcsav_solv_raw.pdb
set all [atomselect top all]
$all moveby [vecinvert [measure center $all]]
display resetview

set solv [atomselect top "segname SOLV"]
$solv set beta 1
set seltext "segname SOLV and same residue as ((z<-20) or (z>10))"
set sel [atomselect top $seltext]
$sel set beta 0

set badwater [atomselect top "name OH2 and beta >0"]
set seglist [$badwater get segid]
set reslist [$badwater get resid]

mol delete all
package require psfgen
resetpsf
topology ../top_all27_prot_lipid.rtf
readpsf kcsav_solv_raw.psf
coordpdb kcsav_solv_raw.pdb
foreach segid $seglist resid $reslist {
delatom $segid $resid
}
writepdb kcsa_solv.pdb
writepsf kcsa_solv.psf
exit
