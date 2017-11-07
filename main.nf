params.gatk = '/gatk-1.1'
params.R_resources = "${params.gatk}/public/R" 
params.picard = '/picard-tools-1.32'
params.platform = 'illumina'


/* 
 * download human genome reference file, reference genome indexing and downloading of 1000Genomes ancillary files
 */
process '0_download' {
  output:
  file 'human_g1k_v37.fasta' into gen_fasta_ch
  file '000G_phase1.indels.b37.vcf' into indels_ch 
  file 'dbsnp_138.b37.excluding_sites_after_129.vcf' into snp_ch
 
  script:
    """
   	wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.gz
	wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.dict.gz
	wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.fai.gz
	wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.excluding_sites_after_129.vcf.gz                
	wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz
	wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_07/low_coverage/snps/CEU.low_coverage.2010_07.genotypes.vcf.gz
	
	gunzip 1000G_phase1.indels.b37.vcf.gz
    gunzip CEU.low_coverage.2010_07.genotypes.vcf.gz
    gunzip dbsnp_138.b37.excluding_sites_after_129.vcf.gz
    gunzip human_g1k_v37.dict.gz
    gunzip human_g1k_v37.fasta.fai.gz
    gunzip human_g1k_v37.fasta.gz
    
    bwa index -a bwtsw human_g1k_v37.fasta #bwt file is <fastaFile.bwt>
    """ 
}

/* 
 * for each sample lane the two fastq files are processed to get the the BAM recalibrated files. All the following tasks are performed in the script "FromFastqToBam.pl".
 * 
 * Quality control of the fastq files
 */
process '1_quality_control' {
  input: 
  set file(sample_lane_1), file(sample_lane_2) from xx
 
  output:
  file 'fqc{1,2}summary.{txt,log}
  file 'sample_out'
  
  script:
  """
  mkdir sample_out
  fastqc $sample_lane_1 -Dfastqc.output_dir=sample_out -Dfastqc.unzip=false 
  fastqc_report_v1.pl -r $sample_lane_1-p sample_out -o fqc1summary.txt -l fqc1summary.log

  fastqc $sample_lane_2 -Dfastqc.output_dir=sample_out -Dfastqc.unzip=false 
  fastqc_report_v1.pl -r $sample_lane_2 -p sample_out -o fqc2summary.txt -l fqc2summary.log 
  """
}


process '2_create_sai_files' {
  input: 
  file gen_fasta from gen_fasta_ch
  set file(sample_lane_1), file(sample_lane_2) from xx

  output:
  set file('*_1.sai'), file('*_2.sai') into sai_ch
  
  script:
  """
  bwa aln $gen_fasta $sample_lane_1 -t $task.cpus -f "sample_1.sai"
  bwa aln $gen_fasta $sample_lane_2 -t $task.cpus -f "sample_2.sai" 
  """
}

process '3_align_to_genome' {
  input:
  set lane, sampleID, library, file(sample_lane_1), file(sample_lane_2) from xx
  set file(sai1), file(sai2) from sai_ch 
  file gen_fasta from gen_fasta_ch
  
  output: 
  file '*.bam' into bam_ch 
  
  // lane: L1
  // sampleID: SRR062634
  // library: HUMjxbRPPDIABPE
  script:
  """
  bwa sampe -P -p $params.platform -i $lane -m $sampleID -l $library $gen_fasta $sai1 $sai2 $sample_lane_1 $sample_lane_2 | \\
  java -Xmx4g -jar ${params.picard}/SamFormatConverter.jar INPUT=/dev/stdin OUTPUT="${sampleID}_${lane}.bam" VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=2000000 TMP_DIR=\$TMPDIR 
  """
}

/* 
 * sorting and indexing of the bam file generated in step 3 
 */
process '4_' {
  input: 
  file bam_file from bam_ch 
  
  output: 
  set file ('*.sorted.bam'), file ('*.sorted.bam.bai') into sorted_ch
  
  script:
  """
  java -Xmx4g -jar ${params.picard}/SortSam.jar INPUT=$bam_file OUTPUT=${bam_file.baseName}.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=${bam_file.baseName}.sorted.bam OUTPUT=${bam_file.baseName}.sorted.bam.bai VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR 
  """
}

/* 
 * removing of optical duplicated from the sorted bam file and subsequent indexing of the duplicates free bam file
 */
process '5_dedup_and_index {
  input: 
  set file(sorted_bam), file(sorted_bai) from sorted_ch
  
  output: 
  set file ('*.dedup.bam'), file ('*.dedup.bam.bai') into dedup_ch
  
  script:
  """
  java -Xmx4g -jar ${params.picard}/MarkDuplicates.jar INPUT=$sorted_bam OUTPUT=${sorted_bam.baseName}.dedup.bam METRICS_FILE=${sorted_bam.baseName}.dedup.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=\$TMPDIR
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=${sorted_bam.baseName}.dedup.bam OUTPUT=${sorted_bam.baseName}.dedup.bam.bai VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR 
  """

} 

/* 
 *  realign reads around indels
 */
process '6_realign_indels' {
  input:
  file gen_file from gen_file_ch
  file indels from indels_ch
  file snp_file from snp_ch
  set file(dedup_bam), file(dedup_bai) from dedup_ch

  output: 
  file '*.realigned.bam' into realigned_ch 
  
  script:
  """
  java -Xmx10g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T RealignerTargetCreator -R $gen_file -I $dedup_bam -o ${dedup_bam.baseName}.intervals
  java -Xmx10g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T IndelRealigner -U ALLOW_UNINDEXED_BAM -I $dedup_bam -targetIntervals ${dedup_bam.baseName}.intervals -known $indels -known $snp_file -o ${dedup_bam.baseName}.realigned.bam -LOD 0.4 -compress 0
  """ 	

} 

/* 
 * fixing mate reads and indexing of the fixed mate bam file
 */
process '7_fixing_and_indexing {
  input:
  file realigned_bam from realigned_ch 
  
  output: 
  set file('*.matefixed.bam'), file('*.matefixed.bam.bai') into matefixed_ch
  
  script:
  """
  java -Xmx4g -jar ${params.picard}/FixMateInformation.jar INPUT=$realigned_bam OUTPUT=${realigned_bam.baseName}.matefixed.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=\$TMPDIR 
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=*.matefixed.bam OUTPUT=${realigned_bam.baseName}.matefixed.bam.bai VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR 
  """
}

/* 
 * bam recalibration and subsequent sorting and indexing
 */
process '8_recalibrate_and_sort' {
  input:
  file gen_file from gen_file_ch
  file snp_file from snp_ch
  set file(matefixed_bam), file(matefixed_bai) from matefixed_ch
  
  output:
  file '*.recal.sorted.bam' into recal_ch
  file '*.matefixed.covariate_table.csv' into matefixed_cov_ch
  script:
  """
  java -Xmx4g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T CountCovariates -U ALLOW_UNINDEXED_BAM -R $gen_file -knownSites $snp_file -I $matefixed_bam -cov ReadGroupcovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile ${matefixed_bam.baseName}.matefixed.covariate_table.csv 
  java -Xmx4g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T TableRecalibration -U ALLOW_UNINDEXED_BAM -R $gen_file -I $matefixed_bam --recal_file *.matefixed.covariate_table.csv --out  ${matefixed_bam.baseName}.recal.bam
  java -Xmx4g -jar ${params.picard}/SortSam.jar INPUT=*.recal.bam OUTPUT=${matefixed_bam.baseName}.recal.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=\$TMPDIR 
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=*.recal.sorted.bam OUTPUT=${matefixed_bam.baseName}.recal.sorted.bam.bai VALIDATION_STRINGENCY=LENIENT TMP_DIR=\$TMPDIR 
  """
}

/* 
 * 	analysis of the recalibration process through comparison of some metrics between MateFixedBam (before recalibration) and recalSortedBam (after recalibration)
 */  
process '9_recalibrate_and_compare' {
  input:
  file gen_file from gen_file_ch
  file snp_file from snp_ch
  file recal_bam from recal_ch
  file matefixed_cov from matefixed_cov_ch

  script:
  """
  mkdir Before
  mkdir After 
  java -Xmx4g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T CountCovariates -U ALLOW_UNINDEXED_BAM -R $gen_file -knownSites $snp_file -I $recal_bam -cov ReadGroupcovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile ${recal_bam.baseName}.recal.covariate_table.csv
  java -Xmx4g -jar ${params.gatk}/AnalyzeCovariates.jar -l INFO -resources ${params.R_resources} --recal_file $matefixed_cov -outputDir Before -Rscript `which R` -ignoreQ 5 
  java -Xmx4g -jar ${params.gatk}/AnalyzeCovariates.jar -l INFO -resources ${params.R_resources} --recal_file *.recal.covariate_table.csv -outputDir After -Rscript `which R` -ignoreQ 5
  """  
} 
 