#/bin/bash
for A in $(seq 1 29)
do
nohup singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash Step2.SAIGE.Ukb_Train_Obesity.chr.32_231011.sh Ukb_Train_cliT2D 1 22 ${A} & > /dev/null
done
