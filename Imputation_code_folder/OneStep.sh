#!/bin/bash
Input=$1
Image="/ichrogene/project/temp/gylee/Singularity"

bash Step1.QC.plink2.sh ${Input}
singularity exec --bind /ichrogene/:/ichrogene/ ${Image}/PRS.sif bash Step2.bcftools_231228.sh ${Input}

singularity exec --bind /ichrogene/:/ichrogene/ ${Image}/shapeit4.2.sif bash Step3.phasing_231228.sh ${Input}

singularity exec --bind /ichrogene/:/ichrogene/ ${Image}/minimac4_1.0.3.sif bash Step4.imputation_minimac4.1.0_231228.sh ${Input}
