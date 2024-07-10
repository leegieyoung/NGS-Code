#!/bin/sh
#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter Sample_Name, thread(s)"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
Dir="/mnt/nas/gylee/0.GWAS"
impute_dir="${Dir}/1.Input/${Sample}/0.Impute"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
Case_pheno="${ana_dir}/temp/11.Case_pheno.txt"
Control_pheno="${ana_dir}/temp/11.Control_pheno.txt"
inversion="${Dir}/REFERENCE/inversion.txt"
Rcode="/mnt/nas/gylee/0.GWAS/Code/R"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
#=================================================================

## Add phenotype
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"

#----------------------------------

for A in $(seq 1 22)
do
${p2_Dir}/plink2 --bfile ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/1000_for_vr_${Sample} \
 --chr ${A} \
 --make-bed \
 --out ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.1000_for_vr_${Sample}
done

for A in $(seq 1 22)
do
${p2_Dir}/plink2 --bfile ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample} \
--chr ${A} \
 --make-bed \
 --out ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}
done
