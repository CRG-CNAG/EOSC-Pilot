<@begin walltime="15:00:00" mem="512mb"/>

sleep 40

#FastQC#

inputs "${datadir}/${file}_1.fq.gz"
inputs "${datadir}/${file}_2.fq.gz"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}_1.fastqcsummary.txt"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}_2.fastqcsummary.txt"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}_1.fastqcsummary.log"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}_2.fastqcsummary.log"
outputs "${outputdir}/${sample}/${file}_1.fq_fastqc.zip"
outputs "${outputdir}/${sample}/${file}_2.fq_fastqc.zip"

mkdir -p "${outputdir}/${sample}/"

${tooldir}/fastqc-v0.7.0/fastqc ${datadir}/${file}_1.fq.gz \
-Dfastqc.output_dir=${outputdir}/${sample} \
-Dfastqc.unzip=false \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}.log

perl ${tooldir}/scripts/fastqc_report_v1.pl \
-r ${datadir}/${file}_1.fq.gz \
-p ${outputdir}/${sample} \
-o ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}.1.fastqcsummary.txt \
-l ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}.1.fastqcsummary.log \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}.log

${tooldir}/fastqc-v0.7.0/fastqc ${datadir}/${file}_2.fq.gz \
-Dfastqc.output_dir=${outputdir}/${sample} \
-Dfastqc.unzip=false \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}.log

perl ${tooldir}/scripts/fastqc_report_v1.pl \
-r ${datadir}/${file}_2.fq.gz \
-p ${outputdir}/${sample} \
-o ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}.2.fastqcsummary.txt \
-l ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}.2.fastqcsummary.log \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al00.fastqc.ftl.${index_file}.log
<@end />