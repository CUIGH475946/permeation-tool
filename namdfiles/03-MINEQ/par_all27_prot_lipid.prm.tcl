
set pi 3.1415926

proc calculate_usage { } {

  vmdcon -info "Usage: calculate <psffile> <dcdfile> -topatm <selection> -lowatm <selection> \[options...\]" 
  vmdcon -info "Calculate small molecule permeation"
  vmdcon -info "Required Parameters:"
  vmdcon -info "  <psffile>       Structure file"
  vmdcon -info "  <dcdfile>       Trajectory file"
  vmdcon -info "  -topatm <sel>   Top atom selection (required)"
  vmdcon -info "  -lowatm <sel>   Low atom selection (required)"
  vmdcon -info "Options:"
  vmdcon -info "  -first <int>    First frame (default:  0)"
  vmdcon -info "  -last  <int>    Last  frame (default: -1)"
  vmdcon -info "  -n <int>        Number of segments (default: 4)"
  vmdcon -info "  -s <int>        Stride value (default: 1)"
  vmdcon -info "  -cal <sel>      Calculation atoms (default: type OT)"
  error ""

}

proc calculate {args} {

  global errorInfo errorCode

  set oldcontext [psfcontext new]
  set errflag [catch { eval calculate_core $args } errMsg]
  set savedInfo $errorInfo
  set savedCode $errorCode
  psfcontext $oldcontext delete
  if {$errflag} { error $errMsg $savedInfo $savedCode }
    
}

proc calculate_core {args} {

  set vmdexec vmd

  # Set some defaults

    # Print usage information if no arguments are given

    if {[llength $args] < 6} { calculate_usage }

    # Confirm the existence of psffile and dcdfile

    lassign $args psffile dcdfile
    set args [lrange $args 2 end]

    if {![file exists $psffile]} { error "psffile not found: $psffile" }
    if {![file exists $dcdfile]} { error "dcdfile not found: $dcdfile" }

    # Check for even number of args

    if {[llength $args] % 2 != 0} { calculate_usage }

    # Specify initialization options

    array set options {

      -first 0    -last -1     -n 4      -s 1
      -cal "type OT"

    }

    array set user_params $args

    set required_params {-topatm -lowatm}

    foreach param $required_params {

      if {![info exists user_params($param)]} { error "Missing required parameter: $param" }

    }

    array set options $args

    proc ::validate_int {val} {string is integer -strict $val}

    foreach {key validator} {

      -first  validate_int    -last   validate_int
      -n      validate_int    -s      validate_int

    } {

      if {![{*}$validator $options($key)]} { error "Invalid value for $key: '$options($key)'" }

    }

  # Start counting

    # Split DCD File

    set step [splitDCD $psffile $dcdfile $options(-first) $options(-last) $options(-n) $options(-s)]

    # Dump scripts

    dumpscripts $psffile $step $options(-n) "$options(-topatm)" "$options(-lowatm)" "$options(-cal)"

    # Parallel run

    parallelrun $options(-n) $vmdexec $step

  # Check and integrate

    # Make check and integration

    check $options(-n) "Data_permeation"
    merge_permeation_data $options(-n)
    analyze_permeation_behaviors Data_permeation_merge.dat Data_permeation.dat

  # Delete temporary files

    for {set i 0} {$i < $options(-n)} {incr i} {

      file delete $i.dcd
      file delete $i.tcl

    }

}

proc splitDCD {psffile dcdfile first last n stride} {

  if {$last != -1 && $last < $first} {error "Invalid frame range: last($last) < first($first)"}

  mol load psf $psffile
  mol addfile $dcdfile type dcd first $first last $last step $stride waitfor all
  set nframe [molinfo top get numframes]

  if {[molinfo top get numframes] == 0} {error "Couldn't load psf/pdb files!"}

  set id [molinfo top get id]

  if {$nframe % $n != 0} {error "The number of frames in each segment must be an integer"}

  set step [expr $nframe / $n]

  vmdcon -info "Splitting DCD into $n segments..."

  for { set i 0 } { $i < $n } { incr i 1 } {

    set nb [expr $i * $step]
    set ne [expr ($i + 1) * $step - 1]
    vmdcon -info "$i.dcd: from frame $nb to frame $ne"
    animate write dcd "$i.dcd" beg $nb end $ne $id

  }

  vmdcon -info "Note: above files will be deleted later."
  mol delete all
  return $step

}

proc dumpscripts {psffile step n topatm lowatm calatm} {

  vmdcon -info "writing scripts for each trajectory ..."

  for {set i 0} {$i < $n} {incr i} {

    set initframe [expr $i * $step]
    set script [open $i.tcl w]

    puts $script "mol load psf $psffile dcd $i.dcd"
    puts $script "source Data_permeation.tcl"
    puts $script "measure_permeation Data_permeation_$i \"$topatm\" \"$lowatm\" \"$calatm\""
    puts $script "mol delete all"
    puts $script "exit"
    puts $script "   "

    close $script
    vmdcon -info "script $i.tcl has been written."

  }

}

proc parallelrun {n vmdexec step} {

  puts "Running tasks ..."
  
  # First (n - 1) tasks will be run in background.
  
  for { set i 0 } { $i < [expr $n - 1] } { incr i 1 } {

    exec vmd -dispdev none -e $i.tcl > /dev/null 2>> /dev/null &

  }

  # The last task will be run in foreground. When this task
  # terminaled, all tasks should have been terminaled.

  set lasttask [expr $n - 1]
  exec $vmdexec -dispdev none -e $lasttask.tcl > /dev/null 2>> /dev/null

  # Wait

  sleep [expr $n * 1.5]
  puts "#--------------------------------------"

}

proc check {n fileid} {

    set flag 0

    for {set i 0} {$i < $n} {incr i} {

      if {[file exists $fileid\_$i.dat]} {

        set flag [expr $flag + 1]

      }

    }

    if {$flag == $n} {

      puts "$fileid tasks finished."

    } else {

      error "$fileid tasks unfinished!!!"

    }

}

proc measure_permeation {outputname topatm lowatm calatm} {

    set fileid [open $outputname\.dat w]
    set nframe [molinfo top get numframes]

    set topsel [atomselect top "$topatm"]
    set lowsel [atomselect top "$lowatm"]
    set calsel [atomselect top "$calatm"]

    array set layers {}

    for {set i 0} {$i < $nframe} {incr i} {

        $topsel frame $i
        $topsel update

        $lowsel frame $i
        $lowsel update

        set top_z [lindex [measure center $topsel] 2]
        set low_z [lindex [measure center $lowsel] 2]

        foreach index [$calsel get index] {

            set sel [atomselect top "index $index"]
            set z   [$sel get z]
            
            if {$z > $top_z} {
                lappend layers($index) 1
            } elseif {$z > $low_z} {
                lappend layers($index) 2
            } else {
                lappend layers($index) 3
            }

            $sel delete
        }
    }

    puts $fileid "$nframe"

    foreach index [$calsel get index] {

        puts $fileid "$index $layers($index)"
    }

    close $fileid

}

proc merge_permeation_data {n} {

    set atom_data {}
    set total_frames 0

    for {set i 0} {$i < $n} {incr i} {

        set filename "Data_permeation_$i.dat"
        set fid [open $filename r]
        
        gets $fid nframe
        incr total_frames $nframe

        while {[gets $fid line] >= 0} {

            set tokens [split $line]
            set index [lindex $tokens 0]
            set values [lrange $tokens 1 end]

            if {[info exists atom_data($index)]} {

                lappend atom_data($index) {*}$values

            } else {

                set atom_data($index) $values
            }
        }
        close $fid
    }

    set fileid [open Data_permeation_merge.dat w]

    puts $fileid $total_frames

    foreach index [lsort -integer [array names atom_data]] {
        puts $fileid "$index $atom_data($index)"
    }

    close $fileid
}


proc analyze_permeation_behaviors {inputfile outputfile {tcut 5}} {

    set fin [open $inputfile r]
    gets $fin total_frames

    set fout [open $outputfile w]

    while {[gets $fin line] >= 0} {
        set tokens [split $line]
        set ntokens [llength $tokens]

        if {$ntokens != [expr {$total_frames + 1}]} {
            continue
        }

        set index [lindex $tokens 0]
        set states [lrange $tokens 1 end]

        for {set i 0} {$i < [expr {$total_frames - 2}]} {incr i} {

            set s0 [lindex $states $i]
            set s1 [lindex $states [expr {$i + 1}]]
            set s2 [lindex $states [expr {$i + 2}]]

            # analyze A: 1 -> 2 -> 3

            if {$s0 == 1 && $s1 == 2} {

                for {set j [expr {$i + 2}]} {$j < $total_frames} {incr j} {

                    if {[lindex $states $j] == 3} {

                        puts $fout "A $index [expr $i + 1] $j [expr $j - $i]"
                        break

                    } elseif {[lindex $states $j] == 1} {
                        break
                    }
                }
            }

            # analyze B: 1 -> 2 -> 1

            if {$s0 == 1 && $s1 == 2} {

                for {set j [expr {$i + 2}]} {$j < $total_frames} {incr j} {

                    if {[lindex $states $j] == 1} {

                        set dwell [expr $j - $i]
                        if {$dwell >= $tcut} {
                            puts $fout "B $index [expr $i + 1] $j $dwell"
                        }
                        break
                    } elseif {[lindex $states $j] == 3} {
                        break
                    }
                }
            }

            # analyze C: 3 -> 2 -> 1
            if {$s0 == 3 && $s1 == 2} {

                for {set j [expr {$i + 2}]} {$j < $total_frames} {incr j} {

                    if {[lindex $states $j] == 1} {

                        puts $fout "C $index [expr $i + 1] $j [expr $j - $i]"
                        break
                    } elseif {[lindex $states $j] == 3} {
                        break
                    }
                }
            }

            # analyze D: 3 -> 2 -> 3
            if {$s0 == 3 && $s1 == 2} {

                for {set j [expr {$i + 2}]} {$j < $total_frames} {incr j} {

                    if {[lindex $states $j] == 3} {

                        set dwell [expr $j - $i]
                        if {$dwell >= $tcut} {
                            puts $fout "D $index [expr $i + 1] $j $dwell"
                        }
                        break
                    } elseif {[lindex $states $j] == 1} {
                        break
                    }
                }
            }
        }
    }

    close $fin
    close $fout
}