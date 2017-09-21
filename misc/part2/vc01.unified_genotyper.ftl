<@begin walltime="28:59:00" ppn="1" mem="10gb" />

sleep 40

inputs "${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam"
inputs "${resdir}/${genome}/indices/${index_file}.fa"
inputs "${resdir}/${genome}/dbsnp/dbsnp_129_b37.rod"
outputs "${outputdir}/${sample}/${sample}.vc01.unified_genotyper.ftl.${index_file}.qc_check_snps.vcf"
outputs "${outputdir}/${sample}/${sample}.vc01.unified_genotyper.ftl.${index_file}.qc_check_snps.vcf.metrics"
outputs "${outputdir}/${sample}/${sample}.vc01.unified_genotyper.ftl.${index_file}.qc_check_snps.concordance.eval"
outputs "${outputdir}/${sample}/${sample}.vc01.unified_genotyper.ftl.${index_file}.qc_check_snps.concordance.txt"

#
# Append path to R lib required for reading GATK reports in R to R_LIB path 
# as it is not installed in the default location for R libs inside the R dir, 
# but in a GATK dir instead.
#
export R_LIBS="${tooldir}/GATK-1.0.5069/Sting/R/"

#Call SNPs on all positions in GoNL ICHIP data
java -Xmx10g -Djava.io.tmpdir=${tempdir} -jar ${tooldir}/GATK-1.0.5069/Sting/dist/GenomeAnalysisTK.jar \
-l INFO \
-T UnifiedGenotyper \
-I ${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam \
--out ${outputdir}/${sample}/${sample}.vc01.unified_genotyper.ftl.${index_file}.qc_check_snps.vcf \
-R ${resdir}/${genome}/indices/${index_file}.fa \
-D ${resdir}/${genome}/dbsnp/dbsnp_129_b37.rod \
-stand_call_conf 40.0 \
-stand_emit_conf 10.0 \
-dcov 200 \
--metrics_file ${outputdir}/${sample}/${sample}.vc01.unified_genotyper.ftl.${index_file}.qc_check_snps.vcf.metrics \
-L ${resdir}/${genome}/intervals/GoNL_SNP_concordance_regions_20110927.interval_list

#-nt 4 \

#Calculate concordance between iChip and VCF data
java -Xmx2g -Djava.io.tmpdir=${tempdir} -jar ${tooldir}/GATK-1.0.5506/Sting/dist/GenomeAnalysisTK.jar \
-T VariantEval \
-B:eval,VCF ${outputdir}/${sample}/${sample}.vc01.unified_genotyper.ftl.${index_file}.qc_check_snps.vcf \
-B:comp_immuno,VCF /target/gpfs2/gcc/home/fvandijk/GoNL_ichip_noAmbiguous/b37/GoNL_Immuno.noAmbiguous.vcf \
-o ${outputdir}/${sample}/${sample}.vc01.unified_genotyper.ftl.${index_file}.qc_check_snps.concordance.eval \
-R ${resdir}/${genome}/indices/${index_file}.fa \
-D ${resdir}/${genome}/dbsnp/dbsnp_129_b37.rod \
-EV GenotypeConcordance

#Retrieve concordance stats and write to file
Rscript ${tooldir}/scripts/extract_info_GATK_variantEval_V2.R \
--in ${outputdir}/${sample}/${sample}.vc01.unified_genotyper.ftl.${index_file}.qc_check_snps.concordance.eval \
--step concordance \
--name ${sample} \
--comp comp_immuno \
--header >> ${outputdir}/${sample}/${sample}.vc01.unified_genotyper.ftl.${index_file}.qc_check_snps.concordance.txt
<@end />