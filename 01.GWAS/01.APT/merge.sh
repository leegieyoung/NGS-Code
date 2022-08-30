#!/bin/sh
APT="/scratch/hpc46a05/raw_data/APT/apt_2.11.4_linux_64_bit_x86_binaries/bin"
analysis_files_path="/scratch/hpc46a05/raw_data/02.Korean_chip/CODA_DATA"
KORV1_1_Lib="/scratch/hpc46a05/raw_data/APT/Library/Axiom_KORV1_1_Analysis.r1"
output="/scratch/hpc46a05/raw_data/APT"

#For Singularity
image_path="/scratch/hpc46a05/image"
Code_path="/scratch/hpc46a05/raw_data/APT/Code"

${APT}/apt-geno-qc \
 --analysis-files-path ${analysis_files_path} \
 --xml-file ${KORV1_1_Lib}/Axiom_KORV1_1.r1.apt-geno-qc.AxiomQC1_nurion.xml \
 --cel-files ${output}/01.geno-qc/CODA_list.txt \
 --out-file ${output}/01.geno-qc/qc.txt \
 --log-file ${output}/01.geno-qc/apt-geno-qc.log

#02.apt-genotype-axiom/CODA_list.txt must remove DQC < 0.82 samples.
${APT}/apt-genotype-axiom \
 --arg-file ${KORV1_1_Lib}/Axiom_KORV1_1_96orMore_Step1.r1.apt-genotype-axiom.AxiomGT1.apt2_nurion.xml \
 --analysis-files-path ${analysis_files_path} \
 --out-dir ${output}/02.apt-genotype-axiom/step1 \
 --dual-channel-normalization true \
 --cel-files ${output}/02.apt-genotype-axiom/CODA_list.txt \
 --log-file ${output}/02.apt-genotype-axiom/02.apt-genotype-axiom.log

#02.apt-genotype-axiom/step1/AxiomGT1.report.txt must remove call_rate < 0.97 samples and write divide1/DQC_CR_pCR_CODA_list.txt.
#Guidelines for passing plates in:
#Plate Pass Rate = of samples passing DQC and 97% QC call rate / Total samples X 100
#High-quality Plates that samples derived from tissue, blood opr cell line Plate pass Rate >= 95%, and >93% if sample source is saliva.
#Calculate the average QC call rate of passing samples on the plate.
#If non-passing plates are identified, all samples from these plates also must be removed in thr process of creating DQC_CR_pCR_CODA_list.txt.
#Average QC call rate of passing samples < 98.5% are of low quality and can result in lower genotyping performance of samples on plates with hihg average QC call rates. So Samples on such plates should be removed from the final genotyping run and considered for reprocessing.
${APT}/apt-genotype-axiom \
 --arg-file ${KORV1_1_Lib}/Axiom_KORV1_1_96orMore_Step2.r1.apt-genotype-axiom.AxiomGT1.apt2_nurion.xml \
 --analysis-files-path ${analysis_files_path} \
 --out-dir ${output}/divide1/step2 \
 --dual-channel-normalization true \
 --summaries \
 --genotyping-node:snp-posteriors-output true \
 --batch-folder ${output}/divide1/suitefiles \
 --multi-genotyping-node:multi-posteriors-output true \
 --cel-files ${output}/divide1/DQC_CR_pCR_CODA_list.txt \
 --log-file ${output}/divide1/03.apt-genotype-axiom.log

#Run ps-metrics
mkdir ${output}/divide1/SNPolisher

${APT}/ps-metrics \
 --posterior-file ${output}/divide1/step2/AxiomGT1.snp-posteriors.txt \
 --call-file ${output}/divide1/step2/AxiomGT1.calls.txt \
 --summary-file ${output}/divide1/step2/AxiomGT1.summary.txt \
 --metrics-file ${output}/divide1/SNPolisher/metrics.txt \
 --log-file ${output}/divide1/snpQC.log

#Run ps-classification
${APT}/ps-classification \
 --species-type human \
 --metrics-file ${output}/divide1/SNPolisher/metrics.txt \
 --output-dir ${output}/divide1/SNPolisher \
 --ps2snp-file ${KORV1_1_Lib}/Axiom_KORV1_1.r1.ps2snp_map.ps \
 --log-file ${output}/divide1/test.snp_classification.log

singularity shell ${image_path}/RNAseq.sif << EOJ
Rscript ${Code_path}/05.snp_classification_divide1_finish.R 
EOJ

#APT-format-result
mkdir ${output}/divide1/plink
pedi="/scratch/hpc46a05/raw_data/APT/01.geno-qc/pedigree.txt"

#Make pedigree.txt
#초기화
rm ${output}/divide1/pedigree.txt
#New pedigree
sed 's/\//\t/g' ${output}/divide1/DQC_CR_pCR_CODA_list.txt | awk '{print $6}' > ${output}/divide1/pedigree.list
for A in $(cat ${output}/divide1/pedigree.list)
do
grep "^${A}" ${pedi} >> ${output}/divide1/pedigree.txt

done

sed -i '1iSample Filename\tFamily_ID\tIndividual_ID\tFather_ID\tMother_ID\tAffection Status\tSex\tPhenotype' ${output}/divide1/pedigree.txt

${APT}/apt-format-result \
 --calls-file ${output}/divide1/step2/snpQC_somatic.txt \
 --annotation-file ${KORV1_1_Lib}/Axiom_KORV1_1.NA35.20210318.annot.db \
 --export-plink-file ${output}/divide1/plink/divide1 \
 --pedigree-file ${output}/divide1/pedigree.txt \
 --log-file ${output}/divide1/apt-format-result.log


