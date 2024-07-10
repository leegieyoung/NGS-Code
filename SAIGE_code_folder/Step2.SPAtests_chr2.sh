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
A=2
Rscript ${SAGE_CD}/step2_SPAtests.R \
        --bedFile=${Output_dir}/chr${A}.imputation_g_m_maf_hwe_bfile.bed \
        --bimFile=${Output_dir}/chr${A}.imputation_g_m_maf_hwe_bfile.bim \
        --famFile=${Output_dir}/chr${A}.imputation_g_m_maf_hwe_bfile.fam \
        --sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
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
        --is_output_moreDetails=TRUE > ${Log}/chr${A}.${Sample}.LDpruned.total.SAIGE.test.markers.log 2>&1

mkdir ${Dir}/2.SAIGE_result/${Sample}/result/temp/
for A in $(seq 2 22)
do
sed '1,1d' ${Dir}/2.SAIGE_result/${Sample}/result/chr${A}.${Sample}.LDpruned.total.SAIGE.test.markers.txt > ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr${A}.merge
done

cat ${Dir}/2.SAIGE_result/${Sample}/result/chr1.${Sample}.LDpruned.total.SAIGE.test.markers.txt \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr2.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr3.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr4.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr5.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr6.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr7.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr8.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr9.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr10.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr11.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr12.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr13.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr14.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr15.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr16.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr17.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr18.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr19.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr20.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr21.merge \
 ${Dir}/2.SAIGE_result/${Sample}/result/temp/chr22.merge > ${Dir}/2.SAIGE_result/${Sample}/result/temp/merge.result
awk '{print $1,$2,$3,$4,$5,$8,$9,$10,$13,$15,$16,$17,$18,$19}' ${Dir}/2.SAIGE_result/${Sample}/result/temp/merge.result > ${Dir}/2.SAIGE_result/${Sample}/result/merge.result

rm -rf ${Dir}/2.SAIGE_result/${Sample}/result/temp/

export RESULT_DIR="${Dir}/2.SAIGE_result/${Sample}/result/"
#echo ${ANA_DIR}
export SAMPLE="${Sample}"
Rscript ${Dir}/Code/R/Manhattan_plot.SAIGE.R



