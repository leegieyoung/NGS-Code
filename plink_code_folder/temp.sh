#!/bin/bash
Sing="/ichrogene/project/temp/gylee/Singularity"
Plink_Dir="/ichrogene/project/temp/gylee/0.GWAS/2.plink_result/Train6Vali4"
Plink_Code_Dir="/ichrogene/project/temp/gylee/0.GWAS/Code/plink_code_folder"
C_list=("C15" "C16" "C18_C19_C20" "C22" "C25" "C50" "C53" "C56" "C61" "C62" "C64" "C67" "C73")
for A in ${C_list[@]}
do
echo ${A}

mkdir ${Plink_Dir}/Ukb_Train_${A}/prune/SHAP
cp ${Plink_Dir}/Ukb_Train_C15/prune/SHAP/shap.C15.edit.R ${Plink_Dir}/Ukb_Train_${A}/prune/SHAP/shap.${A}.R
sed -i "s/C15/${A}/g" ${Plink_Dir}/Ukb_Train_${A}/prune/SHAP/shap.${A}.R

${Sing}/plink2/plink2 --pfile ${Plink_Dir}/Ukb_Train_${A}/imputation_g_m_maf_hwe \
 --extract ${Plink_Dir}/Ukb_Train_${A}/prune/iMAC.PRSice.p1.r10.window250kb.prune.in \
 --recode vcf \
 --out ${Plink_Dir}/Ukb_Train_${A}/prune/${A}.CT

singularity  exec --bind /ichrogene/:/ichrogene/ ${Sing}/AI.sif Rscript ${Plink_Dir}/Ukb_Train_${A}/prune/SHAP/shap.${A}.R
sleep 2

#불필요한 코드인 경우 진행하지 않음
lines_file1=$(wc -l < "${Plink_Dir}/Ukb_Train_${A}/prune/${A}.CT.SHAP.prune.in")
lines_file2=$(wc -l < "${Plink_Dir}/Ukb_Train_${A}/prune/iMAC.PRSice.p1.r10.window250kb.prune.in")

if [ "$lines_file1" -eq "$lines_file2" ]; then
    echo "두 파일의 행 수가 같습니다."
else
    echo "두 파일의 행 수가 다릅니다. 코드를 실행합니다."
singularity  exec --bind /ichrogene/:/ichrogene/ ${Sing}/PRS.sif bash ${Plink_Code_Dir}/3_prune_customCUTOFF_forSAIGE_Ukb_nosingularity_21.1.TrainVali.sh Ukb_Train_${A} 10 10000 ${A}.CT.SHAP ${A} Train6Vali4
sleep 2
fi

done
