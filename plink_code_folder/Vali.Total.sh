#!/bin/bash
Sample=$1
bash 3_prune_customCUTOFF_forSAIGE_KOR.sh KoGES_Train_T2D 32  123123123 ${Sample}
bash 3_prune_customCUTOFF_forSAIGE_KOR_Total.sh KoGES_Train_T2D 32 123412312 ${Sample}
