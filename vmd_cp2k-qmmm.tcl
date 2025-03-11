set auto_path [linsert $auto_path 0 */usr/share/*]
package require math::statistics

puts "1) MAKE SURE YOUR SYSTEM IS THE TOP MOLECULE!"
set topid [molinfo top get id]
puts "CONTROL STRING: TOP MOLECULE IS $topid!"

puts "\n2) INTRODUCE SELECTION FOR QM LAYER IN VMD SELECTION FORMAT:"
set qm_layer [gets stdin]

#--> set INCLUDE file for CP2K SUBSYS section
set subsys_list [list]

#--> set SYSTEM cell
lappend subsys_list "! CP2K > FORCE_EVAL > SUBSYS > CELL"
set cell_vec [measure minmax [atomselect top "noh all"]]
set cell_size [vecsub [lindex $cell_vec 1] [lindex $cell_vec $topid]]
#set cell_vec [pbc get -now]
#set cell_size [lrange [lindex $cell_vec $topid] 0 2]
set cell_shape [lrange [lindex $cell_vec $topid] 3 5]
lappend subsys_list "&CELL\n  ABC \[angstrom\] [join $cell_size "  "]     !size of cell\n  ALPHA_BETA_GAMMA [join $cell_shape " "]      !shape of cell\n&END CELL\n"

#--> set QM layer
set qm_layer_sel [atomselect top $qm_layer frame first]
puts "CONTROL STRING: [$qm_layer_sel num] atoms were selected!"

mol selection "noh $qm_layer"
mol representation DynamicBonds 1.7 0.1
mol material AOChalky
mol addrep top

mol selection "$qm_layer"
mol representation DynamicBonds 1.3 0.1
mol material AOChalky
mol addrep top

set qm_resi [lsort -unique [$qm_layer_sel get resid]]
mol selection "resid $qm_resi and not $qm_layer"
mol representation Lines 3
mol material AOChalky
mol addrep top

#--> set INCLUDE file for CP2K QMMM section
set qmmm_list [list]

#--> set QM cell
lappend qmmm_list "! CP2K > FORCE_EVAL > QMMM > CELL"
set qmcell_vec [measure minmax $qm_layer_sel]
set qmcell_size [vecadd [vecsub [lindex $qmcell_vec 1] [lindex $qmcell_vec $topid]] {10.0 10.0 10.0}]
set qmcell_shape [list 90.0 90.0 90.0]
lappend qmmm_list "&CELL\n  ABC \[angstrom\] [join $qmcell_size "  "]     !size of cell\n  ALPHA_BETA_GAMMA [join $qmcell_shape " "]      !shape of cell\n&END CELL\n"

#--> set QM atoms
lappend qmmm_list "! CP2K > FORCE_EVAL > QMMM > QM_KIND"
set qm_kind [lsort -unique [$qm_layer_sel get element]]
foreach kind $qm_kind {
    mol selection "element $kind and ($qm_layer)"
    mol representation VDW 0.2
    mol material AOChalky
    mol addrep top
    lappend qmmm_list "&QM_KIND $kind\n  MM_INDEX [[atomselect top "element $kind and ($qm_layer)" frame first] get serial]\n&END QM_KIND"
}

#--> set link atoms
lappend qmmm_list "! CP2K > FORCE_EVAL > QMMM > LINK"
set qm_index [$qm_layer_sel get index]
set heavy_qm_index [[atomselect top "noh $qm_layer" frame first] get index]
foreach id $heavy_qm_index {
    foreach conect [lindex [[atomselect top "index $id"] getbonds] $topid] {
        if { [lsearch -exact $qm_index $conect] == -1 } {
            mol selection "index $id $conect"
            mol representation CPK 0.3 0.2
            mol material AOChalky
            mol addrep top
            label add Bonds $topid/$id $topid/$conect
            set resi [[atomselect top "index $id"] get resid]
            set q_link [[atomselect top "index $conect"] get charge]
            set mean_q [::math::statistics::mean [[atomselect top "resid $resi"] get charge]]
            set stdev_q [::math::statistics::stdev [[atomselect top "resid $resi"] get charge]]
            set qmmm_sc [expr 0.83*exp(-0.5*(abs(($q_link-$mean_q)/$stdev_q)**2))]
            set alpha_imomm [expr 1.38+0.3*(0.83-$qmmm_sc)**2]
            lappend qmmm_list "&LINK\n  ALPHA_IMOMM [format "%.2f" $alpha_imomm]     ![join [[atomselect top "index $id $conect"] get type] "-"] bond (default is ALPHA_IMOMM 1.38 for C-C bond)\n  LINK_TYPE IMOMM\n  MM_INDEX  [expr $conect+1]\n  QM_INDEX  [expr $id+1]\n  QMMM_SCALE_FACTOR [format "%.1f" $qmmm_sc]\n&END LINK"
        }
    }
}

#--> write INCLUDE file for CP2K QMMM section
set cp2k_inc [open forceeval_qmmm.inc [list WRONLY CREAT EXCL]]
catch { puts $cp2k_inc [join $qmmm_list "\n"]; close $cp2k_inc }

#--> correct atom charge with parmed
set parmed_list [list]
set molname [molinfo top get name]
lappend parmed_list "parm $molname\n"
puts "\n3) DOUBLE-CHECK THAT INTEGER CHARGE OF SELECTION MATCHES ITS PARCELS!"
set qm_q [expr [join [$qm_layer_sel get charge] +]]
puts "CONTROL STRING: Charge of selection is $qm_q!"

foreach resi $qm_resi {
    set qmresi_q [expr [join [[atomselect top "resid $resi and ($qm_layer)" frame first] get charge] +]]
    set mm_atom [[atomselect top "resid $resi and not ($qm_layer)" frame first] get serial]
    if { ![llength $mm_atom] == 0} {
        set dq [expr ($qmresi_q - round($qmresi_q))/[llength $mm_atom]]
        puts "CONTROL STRING: Charge of the selection in resid $resi is $qmresi_q (applying correction $dq to make it [expr round($qmresi_q)])!"

        foreach atom $mm_atom {
            set q [[atomselect top "serial $atom" frame first] get charge]
            lappend parmed_list "change charge @$atom [expr $q + $dq]"
        }
    }
}
lappend parmed_list "changeLJSingleType :SOL@H* 0.3019 0.047"
lappend parmed_list "changeLJSingleType :WAT@H* 0.3019 0.047"
lappend parmed_list "changeLJSingleType :SER@HG 0.3019 0.047"
lappend parmed_list "changeLJSingleType :THR@HG1 0.3019 0.047"
lappend parmed_list "changeLJSingleType :ASH@HD2 0.3019 0.047"
lappend parmed_list "changeLJSingleType :GLH@HE2 0.3019 0.047"
lappend parmed_list "changeLJSingleType :*@HO 0.3019 0.047"
lappend parmed_list "\noutparm [file root $molname]_ee.prmtop"

#--> write SCRIPT for PARMED charge correction
set parmed_file [open qmmm-ee.parmed.in [list WRONLY CREAT EXCL]]
catch { puts $parmed_file [join $parmed_list "\n"]; close $parmed_file }
