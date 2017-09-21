<@begin walltime="229:59:00" mem="4gb"/>

sleep 40

#Covariates before#

inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al09.fixmates.ftl.${index_file}.matefixed.bam"
inputs "${resdir}/${genome}/indices/${index_file}.fa"
inputs "${resdir}/${genome}/dbsnp/dbsnp_129_b37_${index_file}.rod"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al10.covariates_before.ftl.${index_file}.matefixed.covariate_table.csv"

java -jar -Xmx4g \
${tooldir}/GATK-1.0.5069/Sting/dist/GenomeAnalysisTK.jar -l INFO \
-T CountCovariates \
-U ALLOW_UNINDEXED_BAM \
-R ${resdir}/${genome}/indices/${index_file}.fa \
--DBSNP ${resdir}/${genome}/dbsnp/dbsnp_129_b37_${index_file}.rod \
-I ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al09.fixmates.ftl.${index_file}.matefixed.bam \
-cov ReadGroupcovariate \
-cov QualityScoreCovariate \
-cov CycleCovariate \
-cov DinucCovariate \
-recalFile ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al10.covariates_before.ftl.${index_file}.matefixed.covariate_table.csv \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al10.covariates_before.ftl.${index_file}.log

#Recalibrate#

inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al09.fixmates.ftl.${index_file}.matefixed.bam"
inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al10.covariates_before.ftl.${index_file}.matefixed.covariate_table.csv"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.recal.bam"

java -jar -Xmx4g ${tooldir}/GATK-1.0.5069/Sting/dist/GenomeAnalysisTK.jar -l INFO \
-T TableRecalibration \
-U ALLOW_UNINDEXED_BAM \
-R ${resdir}/${genome}/indices/${index_file}.fa \
-I ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al09.fixmates.ftl.${index_file}.matefixed.bam \
--recal_file ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al10.covariates_before.ftl.${index_file}.matefixed.covariate_table.csv \
--out ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.recal.bam \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.log

#SAM sort

inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.recal.bam"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.recal.sorted.bam"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.recal.sorted.bam.bai"

java -jar -Xmx3g ${tooldir}/picard-tools-1.32/SortSam.jar \
INPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.recal.bam \
OUTPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.recal.sorted.bam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al12.sam_sort.ftl.${index_file}.log

java -jar -Xmx3g ${tooldir}/picard-tools-1.32/BuildBamIndex.jar \
INPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.recal.sorted.bam \
OUTPUT=${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.recal.sorted.bam.bai \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${tempdir} \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al12.sam_sort.ftl.${index_file}.log

#Covariates after#

inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.recal.sorted.bam"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al13.covariates_after.ftl.${index_file}.recal.covariate_table.csv"

java -jar -Xmx4g \
${tooldir}/GATK-1.0.5069/Sting/dist/GenomeAnalysisTK.jar -l INFO \
-T CountCovariates \
-U ALLOW_UNINDEXED_BAM \
-R ${resdir}/${genome}/indices/${index_file}.fa \
--DBSNP ${resdir}/${genome}/dbsnp/dbsnp_129_b37_${index_file}.rod \
-I ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al11.recalibrate.ftl.${index_file}.recal.sorted.bam \
-cov ReadGroupcovariate \
-cov QualityScoreCovariate \
-cov CycleCovariate \
-cov DinucCovariate \
-recalFile ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al13.covariates_after.ftl.${index_file}.recal.covariate_table.csv \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al13.covariates_after.ftl.${index_file}.log

#Analyze covariates#

inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al10.covariates_before.ftl.${index_file}.matefixed.covariate_table.csv"
inputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al13.covariates_after.ftl.${index_file}.recal.covariate_table.csv"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al14.analyze_covariates.ftl.${index_file}.recal.stats_before/${lane}.CycleCovariate.dat"
outputs "${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al14.analyze_covariates.ftl.${index_file}.recal.stats_after/${lane}.CycleCovariate.dat"

java -jar -Xmx4g ${tooldir}/GATK-1.0.5069/Sting/dist/AnalyzeCovariates.jar -l INFO \
-resources ${resdir}/${genome}/indices/${index_file}.fa \
--recal_file ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al10.covariates_before.ftl.${index_file}.matefixed.covariate_table.csv \
-outputDir ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al14.analyze_covariates.ftl.${index_file}.recal.stats_before/ \
-Rscript ${rscript} \
-ignoreQ 5 \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al14.analyze_covariates.ftl.${index_file}.log

java -jar -Xmx4g ${tooldir}/GATK-1.0.5069/Sting/dist/AnalyzeCovariates.jar -l INFO \
-resources ${resdir}/${genome}/indices/${index_file}.fa \
--recal_file ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al13.covariates_after.ftl.${index_file}.recal.covariate_table.csv \
-outputDir ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al14.analyze_covariates.ftl.${index_file}.recal.stats_after/ \
-Rscript ${rscript} \
-ignoreQ 5 \
2>&1 | tee -a ${outputdir}/${sample}/${sample}.${flowcell}_${lane}.al14.analyze_covariates.ftl.${index_file}.log
<@end />