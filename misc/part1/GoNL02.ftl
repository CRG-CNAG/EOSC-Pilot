<@begin mem="6gb" walltime="30:00:00" ppn="4"/>

sleep 40

#Align pair1#

inputs "${resdir}/${genome}/indices/${index_file}.fa" 
inputs "${datadir}/${file}_1.fq.gz"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al01.bwa_align_pair1.ftl.${index_file}.1.sai"

mkdir -p "${outputdir}/${sample}/"

${tooldir}/bwa-0.5.8c_patched/bwa aln \
${resdir}/${genome}/indices/${index_file}.fa \
${datadir}/${file}_1.fq.gz -t ${cores} \
-f ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al01.bwa_align_pair1.ftl.${index_file}.1.sai \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al01.bwa_align_pair1.ftl.${index_file}.log


#Align pair2#

inputs "${resdir}/${genome}/indices/${index_file}.fa" 
inputs "${datadir}/${file}_2.fq.gz"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al02.bwa_align_pair2.ftl.${index_file}.2.sai"

${tooldir}/bwa-0.5.8c_patched/bwa aln \
${resdir}/${genome}/indices/${index_file}.fa \
${datadir}/${file}_2.fq.gz -t ${cores} \
-f ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al02.bwa_align_pair2.ftl.${index_file}.2.sai \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al02.bwa_align_pair2.ftl.${index_file}.log
<@end />