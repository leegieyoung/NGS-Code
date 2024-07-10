#!/bin/bash

if [ $# -ne 1 ];then
        echo "Please enter a folder name."
               exit
fi

INPUT=$1
#Dir="Change2Path"  #Replace the path 
Dir="/ichrogene/project/temp/gylee/0.GWAS/0.IAAP"
REF="/ichrogene/project/temp/gylee/0.GWAS/REFERENCE/IAAP_REF"

echo "======================================================"
echo "               IAAP...ing"
echo "======================================================"
mkdir -p ${Dir}/02.idat2ped/${INPUT}/
mkdir -p ${Dir}/03.make_csv_Input/${INPUT}/
mkdir -p ${Dir}/04.csv_Output/${INPUT}/
mkdir -p ${Dir}/temp

yes | rm ${Dir}/02.idat2ped/${INPUT}/*.ped #anti file collision

iaap-cli gencall ${REF}/ASA-24v1-0_A1.bpm \
 ${REF}/ASA-24v1-0_A1_ClusterFile.egt \
 ${Dir}/02.idat2ped/${INPUT}/ \
 --gender-estimate-call-rate-threshold -0.1 \
 -f ${Dir}/01.idat/${INPUT} -p -pt -t 30

sed h ${Dir}/02.idat2ped/${INPUT}/*.ped > ${Dir}/03.make_csv_Input/${INPUT}/${INPUT}.ped
mv ${Dir}/02.idat2ped/${INPUT}/ASA-24v1-0_A1.map ${Dir}/03.make_csv_Input/${INPUT}/


echo "======================================================"
echo "                   IAAP...Finish"
echo "======================================================"



echo "======================================================"
echo "             Update allele Top to Plus &"
echo "               Sortinng by chr-position"
echo "======================================================"
plink --ped ${Dir}/03.make_csv_Input/${INPUT}/${INPUT}.ped \
 --map ${Dir}/03.make_csv_Input/${INPUT}/ASA-24v1-0_A1.map \
 --recode \
 --keep-allele-order \
 --allow-no-sex \
 --update-alleles ${REF}/update_Top2Plus.txt \
 --out ${Dir}/temp/Noupdate_${INPUT}


echo "======================================================"
echo "                   Sorting ...Finish"
echo "                  Convert to tab file"
echo "======================================================"

awk '{$1=$3=$4=$5=$6=""; sub(/^[[:space:]]+/,""); print $0}' ${Dir}/temp/Noupdate_${INPUT}.ped > ${Dir}/temp/raw1_${INPUT}.ped
awk '{ printf "%s ", $1; for (i=2; i<=NF; i+=2) printf "%s%s ", $(i), $(i+1); printf "\n" }' ${Dir}/temp/raw1_${INPUT}.ped > ${Dir}/temp/raw2_${INPUT}.ped
awk '{ for (i=1; i<=NF; i++) a[i,NR]=$i } NF>p { p=NF } END { for(j=1; j<=p; j++) { str=a[j,1]; for(i=2; i<=NR; i++) str=str" "a[j,i]; print str } }' ${Dir}/temp/raw2_${INPUT}.ped > ${Dir}/temp/back_${INPUT}.ped

paste -d ' ' ${REF}/ASAsnpInfo.txt ${Dir}/temp/back_${INPUT}.ped > ${Dir}/04.csv_Output/${INPUT}/tab.txt

echo "======================================================"
echo "        Convert tab to csv file & Ordering"
echo "======================================================"

export OUTPUT="${Dir}/04.csv_Output/${INPUT}/"
export REF="${REF}/"
Rscript ${Dir}/Code/tab2csv.R
yes | rm ${Dir}/temp/*.bim #Remove used files
yes | rm ${Dir}/temp/*.bed #Remove used files
yes | rm ${Dir}/temp/*.fam #Remove used files
yes | rm ${Dir}/temp/*.ped #Remove used files
yes | rm ${Dir}/temp/*.map #Remove used files


echo "======================================================"
echo "                    Finish by gylee"
echo "======================================================"

