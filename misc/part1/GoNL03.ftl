<@begin walltime="89:59:00" mem="4gb"/>

sleep 40

#BWA SAMpe#

inputs "${resdir}/${genome}/indices/${index_file}.fa"
inputs "${datadir}/${file}_1.fq.gz"
inputs "${datadir}/${file}_2.fq.gz"
inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al01.bwa_align_pair1.ftl.${index_file}.1.sai"
inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al02.bwa_align_pair2.ftl.${index_file}.2.sai"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al03.bwa_sampe.ftl.${index_file}.sam"


${tooldir}/bwa_45_patched/bwa sampe -P \
-p illumina \
-i ${lane} \
-m ${sample} \
-l ${lib} \
${resdir}/${genome}/indices/${index_file}.fa \
${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al01.bwa_align_pair1.ftl.${index_file}.1.sai \
${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al02.bwa_align_pair2.ftl.${index_file}.2.sai \
${datadir}/${file}_1.fq.gz \
${datadir}/${file}_2.fq.gz \
-f ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al03.bwa_sampe.ftl.${index_file}.sam \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al03.bwa_sampe.ftl.${index_file}.log

#SAM to BAM#

inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al03.bwa_sampe.ftl.${index_file}.sam"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al04.sam_to_bam.ftl.${index_file}.bam"

java -jar -Xmx4g ${tooldir}/picard-tools-1.32/SamFormatConverter.jar \
INPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al03.bwa_sampe.ftl.${index_file}.sam \
OUTPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al04.sam_to_bam.ftl.${index_file}.bam \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM=2000000 \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al04.sam_to_bam.ftl.${index_file}.log

#SAM sort#

inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al04.sam_to_bam.ftl.${index_file}.bam"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam.bai"

java -jar -Xmx4g ${tooldir}/picard-tools-1.32/SortSam.jar \
INPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al04.sam_to_bam.ftl.${index_file}.bam \
OUTPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM=1000000 \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_sort.ftl.${index_file}.log

java -jar -Xmx4g ${tooldir}/picard-tools-1.32/BuildBamIndex.jar \
INPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam \
OUTPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam.bai \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM=1000000 \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_sort.ftl.${index_file}.log

#Picard QC#

inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam"
inputs "${resdir}/${genome}/indices/${index_file}.fa"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.AlignmentSummaryMetrics"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.GcBiasMetrics"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.GcBiasMetrics.pdf"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.CollectInsertSizeMetrics"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.CollectInsertSizeMetrics.pdf"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.MeanQualityByCycle"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.MeanQualityByCycle.pdf"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.QualityScoreDistribution"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.QualityScoreDistribution.pdf"

java -jar -Xmx4g ${tooldir}/picard-tools-1.32/CollectAlignmentSummaryMetrics.jar \
I=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam \
O=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.AlignmentSummaryMetrics \
R=${resdir}/${genome}/indices/${index_file}.fa \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.log

java -jar ${tooldir}/picard-tools-1.32/CollectGcBiasMetrics.jar \
R=${resdir}/${genome}/indices/${index_file}.fa \
I=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam \
O=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.GcBiasMetrics \
CHART=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.GcBiasMetrics.pdf \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.log

java -jar ${tooldir}/picard-tools-1.32/CollectInsertSizeMetrics.jar \
I=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam \
O=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.CollectInsertSizeMetrics \
H=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.CollectInsertSizeMetrics.pdf \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.log

java -jar ${tooldir}/picard-tools-1.32/MeanQualityByCycle.jar \
I=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam \
O=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.MeanQualityByCycle \
CHART=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.MeanQualityByCycle.pdf \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.log

java -jar ${tooldir}/picard-tools-1.32/QualityScoreDistribution.jar \
I=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam \
O=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.QualityScoreDistribution \
CHART=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.QualityScoreDistribution.pdf \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.log

java -jar ${tooldir}/picard-tools-1.32/BamIndexStats.jar \
INPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${tempdir} \
> ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.BamIndexStats \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al06.picardQC.ftl.${index_file}.log

#Mark duplicates#

inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam"
inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam.bai"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al07.mark_duplicates.ftl.${index_file}.dedup.bam"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al07.mark_duplicates.ftl.${index_file}.dedup.metrics"

java -Xmx4g -jar ${tooldir}/picard-tools-1.32/MarkDuplicates.jar \
INPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al05.sam_to_bam.ftl.${index_file}.sorted.bam \
OUTPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al07.mark_duplicates.ftl.${index_file}.dedup.bam \
METRICS_FILE=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al07.mark_duplicates.ftl.${index_file}.dedup.metrics \
REMOVE_DUPLICATES=false \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al07.mark_duplicates.ftl.${index_file}.log

java -jar -Xmx4g ${tooldir}/picard-tools-1.32/BuildBamIndex.jar \
INPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al07.mark_duplicates.ftl.${index_file}.dedup.bam \
OUTPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al07.mark_duplicates.ftl.${index_file}.dedup.bam.bai \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM=1000000 \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al07.mark_duplicates.ftl.${index_file}.log
<@end />