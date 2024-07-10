#!/bin/bash
if [ $# -ne 3 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
MEM=$3
Dir="/mnt/nas/gylee/0.GWAS"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
Inversion="/mnt/nas/gylee/0.GWAS/REFERENCE/inversion.txt"
Output_dir="${Dir}/2.plink_result/${Sample}"

mkdir -p ${QC_dir}
mkdir -p ${Output_dir} 
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG/"

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"

#Module
make_pgen() {
 local file="/mnt/nas/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/Ukb_Train_Cancer/${Sample}/${Sample}.IID"
 if [ -e "$file" ]; then
	 echo "파일이 존재합니다. 파일 이름: $file"
	 echo "make_pgen" >> ${Log}/logfile_${0}_${Sample}.txt
         date >> ${Log}/logfile_${0}_${Sample}.txt

cp /mnt/nas/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/Ukb_Train_Cancer/${Sample}/${Sample}_pheno.txt ${Output_dir}/${Sample}_pheno.txt
${p2_Dir}/plink2 --pfile /mnt/nas/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/Ukb_Cancers_64949_R2_0.3_MAF0.005/imputation \
 --keep /mnt/nas/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/Ukb_Train_Cancer/${Sample}/${Sample}.IID \
 --make-pgen \
 --memory ${MEM} \
 --out ${Output_dir}/imputation

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation \
 --geno 0.02 \
 --make-pgen \
 --memory ${MEM} \
 --update-sex /mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/1.raw/subset_of_ukb45411_for_sex.txt \
 --out ${Output_dir}/imputation_g_m

rm -rf ${Output_dir}/imputation.p*

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m \
 --maf 0.01 \
 --make-pgen \
 --memory ${MEM} \
 --out ${Output_dir}/imputation_g_m_maf

rm -rf ${Output_dir}/imputation_g_m.p*

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf \
 --hwe 1e-50 \
 --rm-dup force-first \
 --memory ${MEM} \
 --make-pgen \
 --out ${Output_dir}/imputation_g_m_maf_hwe

rm -rf ${Output_dir}/imputation_g_m_maf.p*
 else
 	echo "파일이 존재하지 않습니다: $file, make_pgen 단계"
	echo "Error : make_pgen, file name : $file"  >> ${Log}/logfile_${0}_${Sample}.txt 
         date >> ${Log}/logfile_${0}_${Sample}.txt
fi
}

make_indep() {
 local file="${Output_dir}/imputation_g_m_maf_hwe.pvar"
 if [ -e "$file" ]; then
	 echo "파일이 존재합니다. 파일 이름: $file"
         echo "make_indep" >> ${Log}/logfile_${0}_${Sample}.txt
         date >> ${Log}/logfile_${0}_${Sample}.txt

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --exclude ${Inversion} \
 --indep-pairwise 50 5 0.02 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --threads ${Thread} \
 --memory ${MEM} \
 --out ${Output_dir}/${Sample}_indepSNP

 else
 	echo "파일이 존재하지 않습니다: $file, indep 단계"
	echo "make_indep, file name : $file" >> ${Log}/logfile_${0}_${Sample}.txt
         date >> ${Log}/logfile_${0}_${Sample}.txt
fi
}

make_pca() {
 local file="${Output_dir}/imputation_g_m_maf_hwe.pvar"
 if [ -e "$file" ]; then
	 echo "파일이 존재합니다. 파일 이름: $file"
         echo "make_pca" >> ${Log}/logfile_${0}_${Sample}.txt
         date >> ${Log}/logfile_${0}_${Sample}.txt

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --double-id \
 --extract ${Output_dir}/${Sample}_indepSNP.prune.in \
 --pca approx 10 \
 --memory ${MEM} \
 --set-missing-var-ids @:# \
 --threads ${Thread} \
 --out ${Output_dir}/${Sample}_PCA

#Make_COV
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/covariate_SI_forSAIGE.R
cp ${Output_dir}/covariate.txt ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/${Sample}.covariate.txt

#Divide hwe file to chr
nohup bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/2_divide.chr.sh ${Sample} ${Thread} ${MEM} > ${Log}/plink.log 2>&1 &

#Divide Prune to chr
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC
${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --extract ${Output_dir}/${Sample}_indepSNP.prune.in \
 --memory ${MEM} \
 --make-bed \
 --out ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/LD_pruned_QC_${Sample}

for A in $(seq 1 22)
do
${p2_Dir}/plink2 --bfile ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/LD_pruned_QC_${Sample} \
 --extract ${Output_dir}/${Sample}_indepSNP.prune.in \
 --memory ${MEM} \
 --make-bed \
 --chr ${A} \
 --out ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/chr${A}.LD_pruned_QC_${Sample}
done

 else
 	echo "파일이 존재하지 않습니다: PCA 단계"
	echo "Error : make_pca, file name : $file" >> ${Log}/logfile_${0}_${Sample}.txt
         date >> ${Log}/logfile_${0}_${Sample}.txt
fi
}

#CreateSparseGRM
make_CreateSparseGRM() {
 local file="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/LD_pruned_QC_${Sample}.bim"
if [ -e "$file" ]; then
 	 echo "파일이 존재합니다. 파일 이름: $file"
         echo "make_CreateSparseGRM" >> ${Log}/logfile_${0}_${Sample}.txt
         date >> ${Log}/logfile_${0}_${Sample}.txt

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/SAIGE_indepGRM_CreateGRM_231017.sh ${Sample} ${Thread}

 else
 	echo "파일이 존재하지 않습니다: CreateSparseGRM 단계"
	echo "Error : make_CreateSparseGRM, file name : $file" >> ${Log}/logfile_${0}_${Sample}.txt
         date >> ${Log}/logfile_${0}_${Sample}.txt
fi
}

#SAIGE
do_SAIGE() {
 local file="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr1.sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
if [ -e "$file" ]; then
	 echo "파일이 존재합니다. 파일 이름: $file"
	 echo "do_SAIGE" >> ${Log}/logfile_${0}_${Sample}.txt
         date >> ${Log}/logfile_${0}_${Sample}.txt

for A in $(seq 1 20)
do
nohup singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/SAIGE_indepGRM2000_2310167.sh ${Sample} ${Thread} ${A} > /dev/null 2>&1 &
done

for A in $(seq 21 22)
do
nohup singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/SAIGE_indepGRM1000_2310167.sh ${Sample} ${Thread} ${A} > /dev/null 2>&1 &
done

wait
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/Step3.merge.sh ${Sample} ${Thread}

 else
 	echo "파일이 존재하지 않습니다: SAIGE 단계"
	echo "Error : do_SAIGE, file name : $file" >> ${Log}/logfile_${0}_${Sample}.txt
         date >> ${Log}/logfile_${0}_${Sample}.txt
fi
}

###########Start###############
rm -rf ${Log}/logfile_${0}_${Sample}.txt
echo "Used Main Code : $0" >> ${Log}/logfile_${0}_${Sample}.txt
date >> ${Log}/logfile_${0}_${Sample}.txt
echo "input file : ${Sample}" >> ${Log}/logfile_${0}_${Sample}.txt
echo "Threads : ${Thread}" >> ${Log}/logfile_${0}_${Sample}.txt
echo "=======================================" >> ${Log}/logfile_${0}_${Sample}.txt
echo "Executed  module" >> ${Log}/logfile_${0}_${Sample}.txt

#make_pgen
make_indep
make_pca

make_CreateSparseGRM
do_SAIGE

