#!/usr/bin/env bash

if [ $# -ne 6 ]; then

    echo "Usage: ./res_del_cp2k.sh <residue_list> <topology> <reactant_structure> <ts_structure> <cp2k_template> <qm_selection>"

else

res_list="$1"
topology="$2"
r_structure=$(echo $3 | sed 's/.pdb//g')
ts_structure=$(echo $4 | sed 's/.pdb//g')
cp2k_input="$5"
qm_selection="$6"

cat <<EOF > cpptraj_del.in
parm ../PRMTOP_TAG
trajin ../STATE_TAG 1
strip :RES_TAG 
trajout res_RES_TAG_FILE_TAG.pdb
EOF

total=$(cat $res_list | wc -l) ; printf "\rProgress: [%-50s] %d/%d" " " 0 $total
counter=0

for i in $(cat $res_list); do
        ((counter++))
        
        mkdir RES_"$i"
        cd RES_"$i"

        cp ../cpptraj_del.in cpptraj_del_"$r_structure".in
	sed -i 's/strip :RES_TAG/strip :RES_TAG parmout res_RES_TAG.prmtop/' cpptraj_del_"$r_structure".in 	
        cp ../cpptraj_del.in cpptraj_del_"$ts_structure".in
        cp ../$cp2k_input res_del_"$r_structure".inp
        cp ../$cp2k_input res_del_"$ts_structure".inp
	
	sed -i 's/PRMTOP_TAG/'"$topology"'/g' cpptraj_del_*.in
        sed -i 's/RES_TAG/'"$i"'/g' cpptraj_del_*.in
        sed -i 's/STATE_TAG/'"$r_structure"'.pdb/g' cpptraj_del_"$r_structure".in
        sed -i 's/STATE_TAG/'"$ts_structure"'.pdb/g' cpptraj_del_"$ts_structure".in
        sed -i 's/FILE_TAG/'"$r_structure"'/g' cpptraj_del_"$r_structure".in
        sed -i 's/FILE_TAG/'"$ts_structure"'/g' cpptraj_del_"$ts_structure".in

        cpptraj -i cpptraj_del_"$r_structure".in >> ../cpptraj.log 2>&1
        cpptraj -i cpptraj_del_"$ts_structure".in >> ../cpptraj.log 2>&1

        vmd res_"$i"_"$r_structure".pdb res_"$i".prmtop -e ../vmd_forceeval.tcl -dispdev none < ../$qm_selection > ../vmd.log 2>&1
	
	qm_charge=$(printf "%.0f\n" `cat qm_charge.dat`)
	sed -i 's/CHARGE .*/CHARGE '"$qm_charge"'/g' res_del_"$r_structure".inp
        sed -i 's/CHARGE .*/CHARGE '"$qm_charge"'/g' res_del_"$ts_structure".inp	

        sed -i 's/STATE_TAG/res_'"$i"_"$r_structure"'.pdb/g' res_del_"$r_structure".inp
        sed -i 's/STATE_TAG/res_'"$i"_"$ts_structure"'.pdb/g' res_del_"$ts_structure".inp
        sed -i 's/PRMTOP_TAG/res_'"$i"'.prmtop/g' res_del_*.inp

	rm cpptraj_del_"$r_structure".in cpptraj_del_"$ts_structure".in qm_charge.dat

        cd ..

        printf "\rProgress: [%-50s] %d/%d" $(printf '#%.0s' $(seq 1 $((counter * 50 / total)))) $counter $total

done

echo ""

rm cpptraj_del.in 

fi
