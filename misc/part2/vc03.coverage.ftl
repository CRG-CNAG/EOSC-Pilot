<@begin walltime="65:59:00" ppn="4" mem="12gb"/>

sleep 40

inputs "${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam"
inputs "${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam.bai"
inputs "${resdir}/${genome}/indices/${index_file}.fa"
outputs "${outputdir}/${sample}/${sample}.vc03.coverage.ftl.${index_file}.bam.coverage.sample_cumulative_coverage_counts"
outputs "${outputdir}/${sample}/${sample}.vc03.coverage.ftl.${index_file}.bam.coverage.sample_cumulative_coverage_proportions"
outputs "${outputdir}/${sample}/${sample}.vc03.coverage.ftl.${index_file}.coverage.pdf"
outputs "${outputdir}/${sample}/${sample}.vc03.coverage.ftl.${index_file}.cumulative_coverage.pdf"

#java -Xmx4g -jar ${tooldir}/Sting/dist/GenomeAnalysisTK.jar \
java -Djava.io.tmpdir=${tempdir} -Xmx12g -jar ${tooldir}/GATK-1.0.5069/Sting/dist/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R ${resdir}/${genome}/indices/${index_file}.fa \
-I ${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam \
-o ${outputdir}/${sample}/${sample}.vc03.coverage.ftl.${index_file}.bam.coverage \
--omitDepthOutputAtEachBase -omitIntervals -omitSampleSummary \
-nt 4


#Create coverage graphs for sample
${rscript} ${tooldir}/scripts/plot_coverage-1.1.R \
--in ${outputdir}/${sample}/${sample}.vc03.coverage.ftl.${index_file}.bam.coverage.sample_cumulative_coverage_counts \
--out ${outputdir}/${sample}/${sample}.vc03.coverage.ftl.${index_file}.coverage.pdf \
--expected-coverage 12 \
--max-depth 40 \
--title "Coverage ${sample}"

${rscript} ${tooldir}/scripts/plot_cumulative_coverage-1.1.R \
--in ${outputdir}/${sample}/${sample}.vc03.coverage.ftl.${index_file}.bam.coverage.sample_cumulative_coverage_proportions \
--out ${outputdir}/${sample}/${sample}.vc03.coverage.ftl.${index_file}.cumulative_coverage.pdf \
--expected-coverage 12 \
--max-depth 40 \
--title "Cumulative coverage ${sample}"



#-geneList ${resdir}/${genome}/indices/genelist.${index_file}.bed \
#mmq 50 and -mbq 20
# http://www.broadinstitute.org/gsa/wiki/index.php/Callable_Loci_Walker
# genelist http://getsatisfaction.com/gsa/topics/refseq_transcript_annotation_file
<@end />