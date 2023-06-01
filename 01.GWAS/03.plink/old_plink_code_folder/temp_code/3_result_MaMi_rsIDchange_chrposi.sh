if [ $# -ne 1 ];then
	echo "please input Imputed_##_##"
		exit
fi

Sample=$1
result_folder="/scratch/x1997a11/GWAS/pdxen_AD/result_folder/${Sample}"

#chr:posi_Major/Minor SNP
awk '{$1=$2="";print $0}' ${result_folder}/MaMi/${Sample}.bim > ${result_folder}/MaMi/raw_other_col_merge_case.bim
awk '{print $1}' ${result_folder}/MaMi/${Sample}.bim > ${result_folder}/MaMi/raw_col_chr_merge_case.bim
awk '{print $4}' ${result_folder}/MaMi/${Sample}.bim > ${result_folder}/MaMi/raw_col_posi_merge_case.bim
awk '{print $6}' ${result_folder}/MaMi/${Sample}.bim > ${result_folder}/MaMi/raw_Major_merge_case.bim
awk '{print $5}' ${result_folder}/MaMi/${Sample}.bim > ${result_folder}/MaMi/raw_Minor_merge_case.bim
mv ${result_folder}/MaMi/${Sample}.bim ${result_folder}/MaMi/original_${Sample}.bim
#chr:posi
paste -d : ${result_folder}/MaMi/raw_col_chr_merge_case.bim ${result_folder}/MaMi/raw_col_posi_merge_case.bim > ${result_folder}/MaMi/raw_12_case.bim

#make bim
paste -d '\t' ${result_folder}/MaMi/raw_col_chr_merge_case.bim ${result_folder}/MaMi/raw_12_case.bim  > ${result_folder}/MaMi/raw_12bim.bim
awk '{print($0"\t'0'")}' ${result_folder}/MaMi/raw_12bim.bim > ${result_folder}/MaMi/raw_123bim.bim
paste -d '\t' ${result_folder}/MaMi/raw_123bim.bim ${result_folder}/MaMi/raw_col_posi_merge_case.bim > ${result_folder}/MaMi/raw_1234bim.bim
paste -d '\t' ${result_folder}/MaMi/raw_1234bim.bim ${result_folder}/MaMi/raw_Minor_merge_case.bim > ${result_folder}/MaMi/raw_12345bim.bim
paste -d '\t' ${result_folder}/MaMi/raw_12345bim.bim  ${result_folder}/MaMi/raw_Major_merge_case.bim > ${result_folder}/MaMi/${Sample}.bim

rm ${result_folder}/MaMi/raw_*
