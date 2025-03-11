#!/usr/bin/env bash

if [ $# -ne 6 ]; then

    echo "Usage: ./res_qmmm_cp2k.sh <residue_list> <topology> <reactant_structure> <ts_structure> <cp2k_template> <qm_selection>"

else

cat <<EOF > cpptraj_pdb.in
parm ../PRMTOP_TAG
trajin ../STATE_TAG 1
strip :RES_TAG
trajout res_RES_TAG_FILE_TAG.pdb
EOF

cat <<EOF > cpptraj_prmtop.in
parm ../PRMTOP_TAG
parmstrip :RES_TAG
parmwrite out res_RES_TAG.prmtop
EOF

total=$(cat $1 | wc -l) ; printf "\rProgress: [%-50s] %d/%d" " " 0 $total
counter=0

for i in $(cat $1); do
        ((counter++))
        
        mkdir RES_"$i"
        cd RES_"$i"

        cp ../cpptraj_pdb.in cpptraj_pdb_"$i"_R.in
        cp ../cpptraj_pdb.in cpptraj_pdb_"$i"_TS.in
        cp ../cpptraj_prmtop.in cpptraj_prmtop_"$i".in
        cp ../$5 res_qmmm_"$(echo "$3" | sed 's/.pdb//g')".inp
        cp ../$5 res_qmmm_"$(echo "$4" | sed 's/.pdb//g')".inp

        sed -i 's/RES_TAG/'"$i"'/g' cpptraj_pdb_"$i"_*.in
        sed -i 's/STATE_TAG/'"$3"'/g' cpptraj_pdb_"$i"_R.in
        sed -i 's/STATE_TAG/'"$4"'/g' cpptraj_pdb_"$i"_TS.in
        sed -i 's/PRMTOP_TAG/'"$2"'/g' cpptraj_pdb_"$i"_R.in
        sed -i 's/PRMTOP_TAG/'"$2"'/g' cpptraj_pdb_"$i"_TS.in
        sed -i 's/FILE_TAG/'"$(echo "$3" | sed 's/.pdb//g')"'/g' cpptraj_pdb_"$i"_R.in
        sed -i 's/FILE_TAG/'"$(echo "$4" | sed 's/.pdb//g')"'/g' cpptraj_pdb_"$i"_TS.in
        sed -i 's/RES_TAG/'"$i"'/g' cpptraj_prmtop_"$i".in
        sed -i 's/PRMTOP_TAG/'"$2"'/g' cpptraj_prmtop_"$i".in

        cpptraj -i cpptraj_pdb_"$i"_R.in >> ../cpptraj.log 2>&1
        cpptraj -i cpptraj_pdb_"$i"_TS.in >> ../cpptraj.log 2>&1
        cpptraj -i cpptraj_prmtop_"$i".in >> ../cpptraj.log 2>&1

        rm cpptraj_pdb_"$i"_R.in cpptraj_pdb_"$i"_TS.in cpptraj_prmtop_"$i".in

        vmd res_"$i"_"$(echo "$3" | sed 's/.pdb//g')".pdb res_"$i".prmtop -e ../vmd_cp2k-qmmm.tcl -dispdev none < ../$6 >> ../vmd.log 2>&1

        sed -i 's/res_'"$i"'_'"$(echo "$3" | sed 's/.pdb//g')"'.pdb/res_'"$i"'.prmtop/g' qmmm-ee.parmed.in
        sed -i 's/res_'"$i"'_'"$(echo "$3" | sed 's/.pdb//g')"'/res_'"$i"'/g' qmmm-ee.parmed.in

        parmed -i qmmm-ee.parmed.in >> ../parmed.log 2>&1

        sed -i 's/STATE_TAG/res_'"$i"'_'"$(echo "$3" | sed 's/.pdb//g')"'.pdb/g' res_qmmm_R.inp
        sed -i 's/STATE_TAG/res_'"$i"'_'"$(echo "$4" | sed 's/.pdb//g')"'.pdb/g' res_qmmm_TS.inp
        sed -i 's/PRMTOP_TAG/res_'"$i"'_ee.prmtop/g' res_qmmm_*.inp

        cd ..

        printf "\rProgress: [%-50s] %d/%d" $(printf '#%.0s' $(seq 1 $((counter * 50 / total)))) $counter $total

done

echo ""

rm cpptraj_pdb.in cpptraj_prmtop.in

fi
