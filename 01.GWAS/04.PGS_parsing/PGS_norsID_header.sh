#/!/bin/sh
for A in $(cat /mnt/nas/gylee/0.GWAS/9.etc/PGS/0.Input/PGS/1343nonrsID.list)
do
#zcat /mnt/nas/PGScatalogScores/${A}.txt.gz | grep '^chr_name' > /mnt/nas/gylee/0.GWAS/9.etc/PGS/1.Output/norsID/${A}_head.txt
zcat /mnt/nas/PGScatalogScores/${A}.txt.gz | grep '^#genome_build' > /mnt/nas/gylee/0.GWAS/9.etc/PGS/1.Output/norsID/${A}_build.txt
done
