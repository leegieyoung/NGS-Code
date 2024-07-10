#!/bin/sh
if [ $# -ne 4 ];then
        echo "Please enter Sample_Name, thread(s)"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
Chr=$3
POS=$4
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
SAIGE_CD="/mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder"

mkdir ${Dir}/2.plink_result/${Sample}/chr${Chr}
cp ${Dir}/2.plink_result/${Sample}/chr${Chr}.imputation_g_m_maf_hwe_bfile.bim ${Dir}/2.plink_result/${Sample}/chr${Chr}/

echo "=========================="
echo "Start Step0.createSparseGRM.sh"
echo "=========================="
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/result

Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG"

export TRAIT1="${Sample}"
export CHR="${Chr}"
Rscript ${SAIGE_CD}/SAIGE.Ukb_Train_Obesity.chr11-64.R

echo "=========================="
echo "Start step2_SPAtests.R"
echo "=========================="
Rscript ${SAIGE_CD}/step2_SPAtests.R \
        --bedFile=${Output_dir}/chr${Chr}/chr${Chr}-${POS}.bed \
        --bimFile=${Output_dir}/chr${Chr}/chr${Chr}-${POS}.bim \
        --famFile=${Output_dir}/chr${Chr}/chr${Chr}-${POS}.fam \
        --sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${Chr}.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${Chr}.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
        --AlleleOrder=alt-first \
        --chrom=${Chr} \
        --GMMATmodelFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${Chr}.${Sample}.LDpruned.total.rda \
        --varianceRatioFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${Chr}.${Sample}.LDpruned.total.varianceRatio.txt \
        --LOCO=FALSE \
	--is_imputed_data=TRUE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.125 \
        --is_output_moreDetails=TRUE \
        --SAIGEOutputFile=${Dir}/2.SAIGE_result/${Sample}/result/chr${Chr}-${POS}.${Sample}.LDpruned.total.SAIGE.test.markers.txt \
        --is_fastTest=TRUE \
        --is_output_moreDetails=TRUE > ${Log}/chr${Chr}-${POS}.${Sample}.LDpruned.total.SAIGE.test.markers.log 2>&1

