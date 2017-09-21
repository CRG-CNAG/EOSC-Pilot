<@begin mem="10gb" walltime="49:59:00"/>

sleep 40

#Realign#

inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al07.mark_duplicates.ftl.${index_file}.dedup.bam"
inputs "${resdir}/${genome}/intervals/realign_intervals_${genome}_${index_file}.intervals"
inputs "${resdir}/${genome}/dbsnp/dbsnp_129_b37_${index_file}.rod"
inputs "${resdir}/${genome}/indels/1kg_pilot_release_merged_indels_sites_${genome}_${index_file}.vcf"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al08.realign.ftl.${index_file}.realigned.bam"

java -Djava.io.tmpdir=${tempdir} -Xmx10g -jar \
${tooldir}/GATK-1.0.5069/Sting/dist/GenomeAnalysisTK.jar -l INFO \
-T IndelRealigner \
-U ALLOW_UNINDEXED_BAM \
-I ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al07.mark_duplicates.ftl.${index_file}.dedup.bam \
-targetIntervals ${resdir}/${genome}/intervals/realign_intervals_${genome}_${index_file}.intervals \
-R ${resdir}/${genome}/indices/${index_file}.fa \
-D ${resdir}/${genome}/dbsnp/dbsnp_129_b37_${index_file}.rod \
-B:indels,VCF ${resdir}/${genome}/indels/1kg_pilot_release_merged_indels_sites_${genome}_${index_file}.vcf \
-o ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al08.realign.ftl.${index_file}.realigned.bam \
-knownsOnly \
-LOD 0.4 \
-compress 0 \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al08.realign.ftl.${index_file}.log

#Fix mates#

inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al08.realign.ftl.${index_file}.realigned.bam"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al09.fixmates.ftl.${index_file}.matefixed.bam"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al09.fixmates.ftl.${index_file}.matefixed.bam.bai"

java -jar -Xmx6g ${tooldir}/picard-tools-1.32/FixMateInformation.jar \
INPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al08.realign.ftl.${index_file}.realigned.bam \
OUTPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al09.fixmates.ftl.${index_file}.matefixed.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al09.fixmates.ftl.${index_file}.log

java -jar -Xmx3g ${tooldir}/picard-tools-1.32/BuildBamIndex.jar \
INPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al09.fixmates.ftl.${index_file}.matefixed.bam \
OUTPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al09.fixmates.ftl.${index_file}.matefixed.bam.bai \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM=1000000 \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al09.fixmates.ftl.${index_file}.log
<@end />