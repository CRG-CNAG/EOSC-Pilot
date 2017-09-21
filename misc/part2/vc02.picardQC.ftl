<@begin mem="4gb" walltime="45:00:00"/>

sleep 40

inputs "${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam"
inputs "${resdir}/${genome}/indices/${index_file}.fa"
outputs "${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.AlignmentSummaryMetrics"
outputs "${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.GcBiasMetrics"
outputs "${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.GcBiasMetrics.pdf"
outputs "${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.CollectInsertSizeMetrics"
outputs "${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.CollectInsertSizeMetrics.pdf"
outputs "${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.MeanQualityByCycle"
outputs "${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.MeanQualityByCycle.pdf"
outputs "${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.QualityScoreDistribution"
outputs "${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.QualityScoreDistribution.pdf"

java -jar -Xmx4g ${tooldir}/picard-tools-1.32/CollectAlignmentSummaryMetrics.jar \
I=${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam \
O=${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.AlignmentSummaryMetrics \
R=${resdir}/${genome}/indices/${index_file}.fa \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=/local

java -jar ${tooldir}/picard-tools-1.32/CollectGcBiasMetrics.jar \
R=${resdir}/${genome}/indices/${index_file}.fa \
I=${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam \
O=${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.GcBiasMetrics \
CHART=${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.GcBiasMetrics.pdf \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=/local

java -jar ${tooldir}/picard-tools-1.32/CollectInsertSizeMetrics.jar \
I=${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam \
O=${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.CollectInsertSizeMetrics \
H=${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.CollectInsertSizeMetrics.pdf \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=/local

java -jar ${tooldir}/picard-tools-1.32/MeanQualityByCycle.jar \
I=${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam \
O=${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.MeanQualityByCycle \
CHART=${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.MeanQualityByCycle.pdf \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=/local

java -jar ${tooldir}/picard-tools-1.32/QualityScoreDistribution.jar \
I=${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam \
O=${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.QualityScoreDistribution \
CHART=${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.QualityScoreDistribution.pdf \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=/local

java -jar ${tooldir}/picard-tools-1.32/BamIndexStats.jar \
INPUT=${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=/local \
> ${outputdir}/${sample}/${sample}.vc02.picardQC.ftl.${index_file}.BamIndexStats
<@end />