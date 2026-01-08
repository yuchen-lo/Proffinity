#!/bin/bash 
#author: Y.C. Lo

rm ./summary.txt

function one2three() {
    
    if [[ $1 == "A" ]]; then
   echo "ALA"
   elif [[ $1 == "R" ]]; then
   echo "ARG"
   elif [[ $1 == "N" ]]; then
   echo "ASN"
   elif [[ $1 == "D" ]]; then
   echo "ASP"
   elif [[ $1 == "C" ]]; then
   echo "CYS"
   elif [[ $1 == "E" ]]; then
   echo "GLU"
   elif [[ $1 == "Q" ]]; then
   echo "GLN"
   elif [[ $1 == "G" ]]; then
   echo "GLY"
   elif [[ $1 == "H" ]]; then
   echo "HIS"
   elif [[ $1 == "I" ]]; then
   echo "ILE"
   elif [[ $1 == "L" ]]; then
   echo "LEU"
   elif [[ $1 == "K" ]]; then
   echo "LYS"
   elif [[ $1 == "M" ]]; then
   echo "MET"
   elif [[ $1 == "P" ]]; then
   echo "PRO"
   elif [[ $1 == "F" ]]; then
   echo "PHE"
   elif [[ $1 == "S" ]]; then
   echo "SER"
   elif [[ $1 == "T" ]]; then
   echo "THR"
   elif [[ $1 == "W" ]]; then
   echo "TRP"
   elif [[ $1 == "Y" ]]; then
   echo "TYR"
   elif [[ $1 == "V" ]]; then
   echo "VAL"
    fi
    
}

function three2one() {
    
    if [[ $1 == "ALA" ]]; then
   echo "A"
   elif [[ $1 == "ARG" ]]; then
   echo "R"
   elif [[ $1 == "ASN" ]]; then
   echo "N"
   elif [[ $1 == "ASP" ]]; then
   echo "D"
   elif [[ $1 == "CYS" ]]; then
   echo "C"
   elif [[ $1 == "GLU" ]]; then
   echo "E"
   elif [[ $1 == "GLN" ]]; then
   echo "Q"
   elif [[ $1 == "GLY" ]]; then
   echo "G"
   elif [[ $1 == "HIS" ]]; then
   echo "H"
   elif [[ $1 == "ILE" ]]; then
   echo "I"
   elif [[ $1 == "LEU" ]]; then
   echo "L"
   elif [[ $1 == "LYS" ]]; then
   echo "K"
   elif [[ $1 == "MET" ]]; then
   echo "M"
   elif [[ $1 == "PHE" ]]; then
   echo "F"
   elif [[ $1 == "PRO" ]]; then
   echo "P"
   elif [[ $1 == "SER" ]]; then
   echo "S"
   elif [[ $1 == "THR" ]]; then
   echo "T"
   elif [[ $1 == "TRP" ]]; then
   echo "W"
   elif [[ $1 == "TYR" ]]; then
   echo "Y"
   elif [[ $1 == "VAL" ]]; then
   echo "V"
    fi
    
}

rm ./log.txt
cat ./skempi_v2.csv | sed 1d | grep -v "n.b" > ./skempi_v2_pro.csv

count=1

while read line; do

pdb_tag=$(echo "$line" | awk -F';' '{print $1}')
echo "pdb_tag: $pdb_tag"
pdb_id=$(echo "$pdb_tag" | awk -F'_' '{print $1}')
echo "pdb_id: $pdb_id"
chain_1=$(echo $pdb_tag | awk -F'_' '{print $2}')
echo "chain_1: $chain_1"
chain_2=$(echo $pdb_tag | awk -F'_' '{print $3}')
echo "chain_2: $chain_2"
kvalue=$(echo "$line" | awk -F';' '{print $8}' | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g')
k2=`echo "-1*l($kvalue*1)/2.303" | bc -l`
echo "affinity:$k2"

echo "$line" | awk -F';' '{print $3}' | tr ',' '\n' > mutant_list.txt

	
	cat ./PDBs/"$pdb_id".mapping > ./temp"$count".txt

	while read line2; do
	
		echo "mutation from $pdb_tag is $line2"
		resi1=$(echo "$line2" | head -c 1)
		resi2=$(echo "$line2" | tail -c 2)
		target_chain=$(echo "$line2" | head -c 2 | tail -c 1)
		resn=$(echo "$line2" | cut -c 3- |sed 's/.$//')
		t=$(one2three "$resi1")
		j=$(one2three "$resi2")
		echo "resi1=$resi1, resi2=$resi2, resn=$resn, target_chain=$target_chain, t=$t, j=$j" 
		cat ./temp"$count".txt | awk -v resn="$resn" -v chain="$target_chain" -v resi="$j" '$4==resn && $2==chain {$1=resi}1' >./temp.txt
		cat ./temp.txt > ./temp"$count".txt
	
	
	done < ./mutant_list.txt
	
	rm ./seq"$count".txt
	
	
	start_chain=$(cat ./temp"$count".txt | head -1 | awk '{print $2}')
	start_resi=$(cat ./temp"$count".txt | head -1 | awk '{print $4}')
	recent_chain=$start_chain
	
	while read line3; do
	
		echo "line3:$line3"
		tlc=$(echo "$line3" | awk '{print $1}')
		echo "tlc:$tlc"
		olc=$(three2one "$tlc")
		echo "olc:$olc"
		echo "$olc" >> ./seq$count.txt
		current_chain=$(echo $line3 | awk '{print $2}')
		current_resi=$(echo $line3 | awk '{print $4}')
	
		echo "pdb_tag: $pdb_tag, start_chain:$start_chain, current_chain=$current_chain" >> ./log.txt
		echo "pdb_tag: $pdb_tag, start_resi:$start_resi, current_resi=$current_resi" >> ./log.txt
	
		if [ $recent_chain == $current_chain ]; then
  			echo "same chain"
  		else
  			echo "different chain"
  			echo "/" >> ./seq$count.txt
  			recent_chain=$current_chain
		fi
	
	
	done < ./temp"$count".txt
	
	echo "*" >> ./seq$count.txt
	
	rm ./temp"$count".txt
	echo ">P1;$pdb_id" >> ./header.txt
	echo ">P1;$pdb_tag" > ./header2.txt
	echo "sequence:$pdb_id:$start_resi:$start_chain:$current_resi:$current_chain::::" >> ./header.txt
	cat ./header.txt >> ./log.txt
	cat seq$count.txt | awk '{print $1}' | tr '\n' ',' | sed 's/,//g' > ./seq"$count"_out.txt
	cat seq$count.txt | awk '{print $1}' | tr '\n' ',' | sed 's/,//g' | tr '/' ':' | sed 's/*//g' > ./seq"$count"_out2.txt
	cat ./header.txt ./seq"$count"_out.txt > ./seq"$count"_out.ali
	cat ./header2.txt ./seq"$count"_out2.txt > ./af_seq_output/seq"$count"_out2.fasta
	rm ./header.txt
	
	

#modeller

model_number=1

cat align2d-template.py | sed 's/pdb_id/'$pdb_id'/g' | sed 's/chain_1/'$chain_1'/g'|sed 's/chain_2/'$chain_2'/g' | sed 's/pdbfile/'$pdb_id'.pdb/g' | sed 's/seqX_out.ali/seq'$count'_out.ali/g' | sed 's/count.ali/model'$count'.ali/g' | sed 's/count.pap/model'$count'.pap/g' | sed 's/pdb_tag/'$pdb_tag'/g' > align2d-model"$count".py

python3 ./align2d-model"$count".py


cat model-single-template.py | sed 's/pdb_id/'$pdb_id'/g' | sed 's/pdb_tag/'$pdb_tag'/g' | sed 's/modelX.ali/model'$count'.ali/g' | sed 's/model_number/'$model_number'/g'  > model-single-model"$count".py

report=$(python3 ./model-single-model"$count".py)

echo "$report" > ./report.txt


cat ./report.txt | tail -"$model_number" > ./model"$count"_score_rank.txt

best_model_name=$(cat ./model"$count"_score_rank.txt | awk '{print $1}' | head -1)
best_model_score=$(cat ./model"$count"_score_rank.txt | awk '{print $2}' | head -1)
best_model_score2=`echo "-1*$best_model_score" | bc -l`
echo "$pdb_tag,SKEMPI_$count,$best_model_score2,$k2" >> ./summary.txt

cat $best_model_name > ./modeller_output/best_model_$count.pdb

count=$(echo "$count + 1" | bc)	



done < ./skempi_v2_binary.csv











