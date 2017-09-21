<@begin walltime="45:59:00" ppn="1" mem="4gb"/>

sleep 40

inputs "${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam" 
inputs "${resdir}/${genome}/indices/${index_file}.fa"
###Change when chipconcordance etc. is automated###
#inputs "/target/gpfs2/gcc/home/fvandijk/GoNL_ichip_noAmbiguous/b37/GoNL_Immuno.noAmbiguous.bed"
#inputs "/target/gpfs2/gcc/home/fvandijk/GoNL_ichip_noAmbiguous/b37/GoNL_Immuno.noAmbiguous.bim"
#inputs "/target/gpfs2/gcc/home/fvandijk/GoNL_ichip_noAmbiguous/b37/GoNL_Immuno.noAmbiguous.fam"

inputs "/target/gpfs2/gcc/home/fvandijk/GoNL_ichip_noAmbiguous/b37/GoNL_Immuno.noAmbiguous.bed"
inputs "/target/gpfs2/gcc/home/fvandijk/GoNL_ichip_noAmbiguous/b37/GoNL_Immuno.noAmbiguous.bim"
inputs "/target/gpfs2/gcc/home/fvandijk/GoNL_ichip_noAmbiguous/b37/GoNL_Immuno.noAmbiguous.fam"

outputs "${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.best.merged"
outputs "${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.bestRG"
outputs "${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.bestSM"
outputs "${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.self.merged"
outputs "${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.selfRG"
outputs "${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.selfSM"

##Tool doesn't work correctly. Make symlink to reference first and index the symlink file
ln -s ${resdir}/${genome}/indices/${index_file}.fa ${outputdir}/${sample}/verifyBamID.fa


${tooldir}/verifyBamID-0.0.5/verifyBamID/verifyBamID \
--reference ${outputdir}/${sample}/verifyBamID.fa \
--in ${outputdir}/${sample}/${sample}.vc00.merge.ftl.${index_file}.bam \
--bfile /target/gpfs2/gcc/home/fvandijk/GoNL_ichip_noAmbiguous/b37/GoNL_Immuno.noAmbiguous \
--out ${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file} \
--maxDepth 30 \
--verbose

unlink ${outputdir}/${sample}/verifyBamID.fa
rm  ${outputdir}/${sample}/verifyBamID-bs.umfa

cat ${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.bestRG \
${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.bestSM > \
${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.best.merged

cat ${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.selfRG \
${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.selfSM > \
${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.self.merged

#Check for contamination
#perl ${tooldir}/scripts/check_contamination.pl \
#${sample} \
#${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.self.merged \
#${outputdir}/${sample}/${sample}.vc03b.contamination_checker.${index_file}.best.merged
<@end />