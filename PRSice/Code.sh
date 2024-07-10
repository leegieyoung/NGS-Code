singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript PRSice.R \
 --prsice ./PRSice_linux \
 --base TOY_BASE_GWAS.assoc \
 --target TOY_TARGET_DATA \
 --thread 4 \
 --stat OR --binary-target T
