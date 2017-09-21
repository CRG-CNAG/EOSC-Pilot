<@begin mem="6gb" walltime="75:59:00" queue="gcc" ppn="2"/>

sleep 40

<#list lanes as lane>inputs "${outputdir}/${lane.getString("sample")}/${lane.getString("sample")}.${lane.getString("flowcell")}_${lane.getString("lane")}.al11.recalibrate.ftl.${index_file}.recal.sorted.bam"
</#list>
outputs "${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam"
outputs "${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam.bai"

java -jar -Xmx6g ${tooldir}/picard-tools-1.32/MergeSamFiles.jar \
<#list lanes as lane>INPUT=${outputdir}/${lane.getString("sample")}/${lane.getString("sample")}.${lane.getString("flowcell")}_${lane.getString("lane")}.al11.recalibrate.ftl.${index_file}.recal.sorted.bam \
</#list>
ASSUME_SORTED=true USE_THREADING=true \
TMP_DIR=${tempdir} MAX_RECORDS_IN_RAM=6000000 \
OUTPUT=${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=SILENT

java -jar -Xmx3g ${tooldir}/picard-tools-1.32/BuildBamIndex.jar \
INPUT=${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam \
OUTPUT=${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam.bai \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM=1000000 \
TMP_DIR=${tempdir}
<@end />