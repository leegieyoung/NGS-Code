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
Output_dir="${Dir}/2.plink_result/${Sample}"
#=================================================================
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"
SAGE_CD="/mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder"

mkdir -p ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/result
Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG/"
 
#SV-Step2
for A in $(seq 7 22)
do
nohup Rscript ${SAGE_CD}/step2_SPAtests.R \
        --bedFile=${Output_dir}/chr${A}.imputation_g_m_maf_hwe_bfile.bed \
        --bimFile=${Output_dir}/chr${A}.imputation_g_m_maf_hwe_bfile.bim \
        --famFile=${Output_dir}/chr${A}.imputation_g_m_maf_hwe_bfile.fam \
        --sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx \
        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
        --AlleleOrder=alt-first \
       --chrom=${A} \
        --GMMATmodelFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${A}.${Sample}.LDpruned.total.rda \
        --varianceRatioFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${A}.${Sample}.LDpruned.total.varianceRatio.txt \
        --LOCO=FALSE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.125 \
        --is_output_moreDetails=TRUE \
        --SAIGEOutputFile=${Dir}/2.SAIGE_result/${Sample}/result/chr${A}.${Sample}.LDpruned.total.SAIGE.test.markers.txt \
        --is_fastTest=TRUE \
        --is_output_moreDetails=TRUE > ${Log}/chr${A}.${Sample}.LDpruned.total.SAIGE.test.markers.log 2>&1 &
done


