#/bin/sh
A=$1
mkdir -p SMK084.chr${A}/workDir
samtools view -b SMK-084-1-TR.STARAlignedAligned.sortedByCoord.out.markduplicates.Split.withRG.BQSR.bam chr${A} > SMK084.chr${A}/SMK084.chr${A}.bam
samtools index SMK084.chr${A}/SMK084.chr${A}.bam
singularity exec /data/sskimb/Singularity/weCall.sif  weCall --inputs SMK084.chr${A}/SMK084.chr${A}.bam \
 --output SMK084.chr${A}/SMK084.chr${A}.vcf \
 --refFile Homo_sapiens_assembly38_reducted.fa \
# --numberOfJobs=16 \
 --workDir SMK084.chr${A}/workDir/ \
 --logFilename SMK084.chr${A}/SMK084.chr${A}.log 

