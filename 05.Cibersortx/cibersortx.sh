singularity shell --bind /data/keeyoung/scRNA/cibersortx/input:/src/data --bind /data/keeyoung/scRNA/cibersortx/output:/src/outdir /data/sskimb/Singularity/cibersortx_fractions.sif
#cd /src
#CRC
##Make Signature Matrix
#iCD_cpm_220923_Only_CodingGene_addBcell
singularity shell --bind /data/keeyoung/scRNA/cibersortx/input/iCD/Only_CodingGene_220914/cpm:/src/data --bind /data/keeyoung/scRNA/cibersortx/output/01.SignatureMatrix/iCD/cpm_OnlyCodingGene_Bcell_220923:/src/outdir /data/sskimb/Singularity/cibersortx_fractions.sif

/src/CIBERSORTxFractions --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --refsample /src/data/iCD_cpm_Only_CodingGene_addBcell_220923.txt --single_cell TRUE

#iCD_rpkm_221201_Only_CodingGene_addBcell
singularity shell --bind /data/keeyoung/scRNA/cibersortx/input/iCD/Only_CodingGene_220914/rpkm:/src/data --bind /data/keeyoung/scRNA/cibersortx/output/01.SignatureMatrix/iCD/rpkm_OnlyCodingGene_Bcell_221201:/src/outdir /data/sskimb/Singularity/cibersortx_fractions.sif

/src/CIBERSORTxFractions --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --refsample /src/data/iCD_rpkm_Only_CodingGene_addBcell_rmUnknown_221201.txt --single_cell TRUE

#CRC_GSE132465
singularity shell --bind /data/keeyoung/scRNA/cibersortx/input/CRC/OnlyCodingGene_220919:/src/data --bind /data/keeyoung/scRNA/cibersortx/output/01.SignatureMatrix/GSE132465/OnlyCodingGene_220919:/src/outdir /data/sskimb/Singularity/cibersortx_fractions.sif

/src/CIBERSORTxFractions --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --refsample /src/data/GSE132465_cpm_220919.txt --single_cell TRUE

#Make Cell fraction by S-mode

#IBD #Norm, Protein_Coding_Gene, addBcells, Consider platform
singularity shell --bind /data/keeyoung/scRNA/cibersortx/input/iCD/Only_CodingGene_220914/cpm:/src/data --bind /data/keeyoung/scRNA/cibersortx/output/02.Imputed_Cell_Fraction/IBD/IBD_cpm_221102:/src/outdir /data/sskimb/Singularity/cibersortx_fractions.sif

/src/CIBERSORTxFractions --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --mixture /src/data/IBD_cpm_221102.txt --sigmatrix /src/data/SM_cpm_Only_CodingGene_addBcell_220923.txt --perm 100 --refsample /src/data/iCD_cpm_Only_CodingGene_addBcell_220923.txt --rmbatchSmode TRUE

#IBD #Norm, Protein_Coding_Gene, addBcells, Totalseq
singularity shell --bind /data/keeyoung/scRNA/cibersortx/input/iCD/Only_CodingGene_220914/cpm:/src/data --bind /data/keeyoung/scRNA/cibersortx/output/02.Imputed_Cell_Fraction/IBD/IBD_Total_cpm_221103:/src/outdir /data/sskimb/Singularity/cibersortx_fractions.sif

/src/CIBERSORTxFractions --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --mixture /src/data/IBD_Total_cpm_221103.txt --sigmatrix /src/data/SM_cpm_Only_CodingGene_addBcell_220923.txt --perm 100 --refsample /src/data/iCD_cpm_Only_CodingGene_addBcell_220923.txt --rmbatchSmode TRUE

#IBD #Norm, Protein_Coding_Gene, addBcells, Truseq
singularity shell --bind /data/keeyoung/scRNA/cibersortx/input/iCD/Only_CodingGene_220914/cpm:/src/data --bind /data/keeyoung/scRNA/cibersortx/output/02.Imputed_Cell_Fraction/IBD/IBD_Tru_cpm_221103:/src/outdir /data/sskimb/Singularity/cibersortx_fractions.sif
/src/CIBERSORTxFractions --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --mixture /src/data/IBD_Tru_cpm_221103.txt --sigmatrix /src/data/SM_cpm_Only_CodingGene_addBcell_220923.txt --perm 100 --refsample /src/data/iCD_cpm_Only_CodingGene_addBcell_220923.txt --rmbatchSmode TRUE

#RISK_199 #Nrom Protein_Coding_gene, addBcells, cpm
singularity shell --bind /data/keeyoung/scRNA/cibersortx/input/iCD/Only_CodingGene_220914/cpm:/src/data --bind /data/keeyoung/scRNA/cibersortx/output/02.Imputed_Cell_Fraction/IBD/RISK_199_cpm_221206:/src/outdir /data/sskimb/Singularity/cibersortx_fractions.sif
/src/CIBERSORTxFractions --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --mixture /src/data/RISK_199_cpm_Norm_221206.txt --sigmatrix /src/data/SM_cpm_Only_CodingGene_addBcell_220923.txt --perm 100 --refsample /src/data/iCD_cpm_Only_CodingGene_addBcell_220923.txt --rmbatchSmode TRUE

#RISK #Norm, Protein_Coding_Gene, addBcells, fpkm(rpkm) #wrong
singularity shell --bind /data/keeyoung/scRNA/cibersortx/input/iCD/Only_CodingGene_220914/rpkm:/src/data --bind /data/keeyoung/scRNA/cibersortx/output/02.Imputed_Cell_Fraction/IBD/RISK_rpkm_221201:/src/outdir /data/sskimb/Singularity/cibersortx_fractions.sif

/src/CIBERSORTxFractions --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --mixture /src/data/GSE134881_rmdupli_LOC_RISKfpkmAll_Coding_rm_min200.txt --sigmatrix /src/data/SM_iCD_rpkm_Only_CodingGene_addBcell_rmUnknown_221201.txt --perm 100 --refsample /src/data/iCD_rpkm_Only_CodingGene_addBcell_rmUnknown_221201.txt --rmbatchSmode TRUE

/src/CIBERSORTxFractions --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --mixture /src/data/GSE134881_rmdupli_LOC_RISKfpkmAll_Coding_rm_min200_change1e-6.txt --sigmatrix /src/data/SM_iCD_rpkm_Only_CodingGene_addBcell_rmUnknown_221201.txt --perm 100 --refsample /src/data/iCD_rpkm_Only_CodingGene_addBcell_rmUnknown_221201.txt --rmbatchSmode TRUE

##TEST
RISK_cpm_Norm_221130.txt
/src/CIBERSORTxFractions --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --mixture /src/data/RISK_cpm_Norm_221130.txt --sigmatrix /src/data/SM_iCD_rpkm_Only_CodingGene_addBcell_rmUnknown_221201.txt --perm 100 --refsample /src/data/iCD_rpkm_Only_CodingGene_addBcell_rmUnknown_221201.txt --rmbatchSmode TRUE

#03.GEP
#Groud Truth 221108_IBD
singularity shell --bind /data/keeyoung/scRNA/cibersortx/input/iCD/Only_CodingGene_220914/cpm:/src/data --bind /data/keeyoung/scRNA/cibersortx/output/03.GEP/CD_infla/cpm_OnlyCodingGene_Bcell_Smode_GEP_221108:/src/outdir /data/sskimb/Singularity/cibersortx_gep.sif

Rscript /src/R_modules/CIBERSORTxGEP.R --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --mixture /src/data/IBD_cpm_221102.txt --sigmatrix /src/data/SM_cpm_Only_CodingGene_addBcell_220923.txt --refsample /src/data/iCD_cpm_Only_CodingGene_addBcell_220923.txt --rmbatchSmode TRUE --groundtruth /src/data/iCD_cpm_Only_CodingGene_addBcell_220923.txt --threads 20

#Groud Truth 221115 Total
singularity shell --bind /data/keeyoung/scRNA/cibersortx/input/iCD/Only_CodingGene_220914/cpm:/src/data --bind /data/keeyoung/scRNA/cibersortx/output/03.GEP/CD_infla/cpm_OnlyCodingGene_Bcell_Smode_GEP_Total_221115:/src/outdir /data/sskimb/Singularity/cibersortx_gep.sif

Rscript /src/R_modules/CIBERSORTxGEP.R --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --mixture /src/data/IBD_Total_cpm_221103.txt --sigmatrix /src/data/SM_cpm_Only_CodingGene_addBcell_220923.txt --refsample /src/data/iCD_cpm_Only_CodingGene_addBcell_220923.txt --rmbatchSmode TRUE --groundtruth /src/data/iCD_cpm_Only_CodingGene_addBcell_220923.txt --threads 20

#Groud Truth 221115 Tru
singularity shell --bind /data/keeyoung/scRNA/cibersortx/input/iCD/Only_CodingGene_220914/cpm:/src/data --bind /data/keeyoung/scRNA/cibersortx/output/03.GEP/CD_infla/cpm_OnlyCodingGene_Bcell_Smode_GEP_Tru_221115:/src/outdir /data/sskimb/Singularity/cibersortx_gep.sif

Rscript /src/R_modules/CIBERSORTxGEP.R --username gylee@soongsil.ac.kr --token 1f766aa1888a03dc415d7d023548af8b --mixture /src/data/IBD_cpm_221102.txt --sigmatrix /src/data/SM_cpm_Only_CodingGene_addBcell_220923.txt --refsample /src/data/iCD_cpm_Only_CodingGene_addBcell_220923.txt --rmbatchSmode TRUE --groundtruth /src/data/iCD_cpm_Only_CodingGene_addBcell_220923.txt --threads 20

