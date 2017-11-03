params.gatk = '/gatk-1.1'
params.picard = '/picard-tools-1.32'


/* 
 * download human genome reference file, reference genome indexing and downloading of 1000Genomes ancillary files
 */
process '0_download' {
 
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


  script:
  """
  fastqc <Sample_Lane_1.fq.gz> -Dfastqc.output_dir=<SampleLaneOutputFolder> -Dfastqc.unzip=false 
  fastqc_report_v1.pl -r <Sample_Lane_1.fq.gz> -p <SampleLaneOutputFolder> -o <fqc1summary> -l <fqc1summaryLog>

  fastqc <Sample_Lane_2.fq.gz> -Dfastqc.output_dir=<SampleLaneOutputFolder> -Dfastqc.unzip=false 
  fastqc_report_v1.pl -r <Sample_Lane_2.fq.gz> -p <SampleLaneOutputFolder> -o <fqc2summary> -l <fqc2summaryLog> 
  """
}


process '2_create_sai_files' {

  script:
  """
  bwa aln <fastaFile> <Sample_Lane_1.fq.gz> -t $task.cpus -f <sai1> 
  bwa aln <fastaFile> <Sample_Lane_2.fq.gz> -t $task.cpus -f <sai2> 

  """
}

process '3_align_to_genome' {

  script:
  """
  bwa sampe -P -p <platform> -i <lane> -m <sampleID> -l <library> <fastaFile> <sai1> <sai2> <fq1> <fq2> | \\
  java -Xmx4g -jar ${params.picard}/SamFormatConverter.jar INPUT=/dev/stdin OUTPUT=<bam> VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=2000000 TMP_DIR=\$TMPDIR 
  """
}

/* 
 * sorting and indexing of the bam file generated in step 3 
 */
process '4_' {
  
  script:
  """
  java -Xmx4g -jar ${params.picard}/SortSam.jar INPUT=<bam> OUTPUT=<sortedBam> SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=<sortedBam> OUTPUT=<sortedBamIndex> VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR 
  """
}

/* 
 * removing of optical duplicated from the sorted bam file and subsequent indexing of the duplicates free bam file
 */
process '5_dedup_and_index {

  script:
  """
  java -Xmx4g -jar ${params.picard}/MarkDuplicates.jar INPUT=<sortedBam> OUTPUT=<NodupBam> METRICS_FILE=<NodupMetrics> REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=\$TMPDIR
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=<NodupBam> OUTPUT=<NodupBamIndex> VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR 
  """

} 

/* 
 *  realign reads around indels
 */
process '6_realign_indels' {

  script:
  """
  java -Xmx10g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T RealignerTargetCreator -R <fastaFile> -I <NodupBam> -o <intervalFile>
  java -Xmx10g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T IndelRealigner -U ALLOW_UNINDEXED_BAM -I <NodupBam> -targetIntervals <intervalFile> -known <IndelsFile> -known dbSNPfile -o <realignedBam> -LOD 0.4 -compress 0
  """ 	

} 

/* 
 * fixing mate reads and indexing of the fixed mate bam file
 */
process '7_fixing_and_indexing {

  script:
  """
  java -Xmx4g -jar ${params.picard}/FixMateInformation.jar INPUT=<realignedBam> OUTPUT=<MateFixedBam> SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=\$TMPDIR 
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=<MateFixedBam> OUTPUT=<MateFixedIndex> VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR 
  """
}

/* 
 * bam recalibration and subsequent sorting and indexing
 */
process '8_recalibrate_and_sort' {

  script:
  """
  java -Xmx4g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T CountCovariates -U ALLOW_UNINDEXED_BAM -R <fastaFile> -knownSites <dbSNPfile> -I <MateFixedBam> -cov ReadGroupcovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile <MateFixedCovariate> 
  java -Xmx4g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T TableRecalibration -U ALLOW_UNINDEXED_BAM -R <fastaFile> -I <MateFixedBam> --recal_file <MateFixedCovariate> --out <recalBam> 
  java -Xmx4g -jar ${params.picard}/SortSam.jar INPUT=<recalBam> OUTPUT=<recalSortedBam> SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=\$TMPDIR 
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=<recalSortedBam> OUTPUT=<recalSortedIndex> VALIDATION_STRINGENCY=LENIENT TMP_DIR=\$TMPDIR 
  """
}

/* 
 * 	analysis of the recalibration process through comparison of some metrics between MateFixedBam (before recalibration) and recalSortedBam (after recalibration)
 */  
process '9_recalibrate_and_compare' {

  script:
  """
  java -Xmx4g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T CountCovariates -U ALLOW_UNINDEXED_BAM -R <fastaFile> -knownSites <dbSNPfile> -I <recalSortedBam> -cov ReadGroupcovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile <RecalCovariate> 
  java -Xmx4g -jar ${params.gatk}/AnalyzeCovariates.jar -l INFO -resources <Rresources> --recal_file <MateFixedCovariate> -outputDir <StatsBeforeFolder> -Rscript <RscriptPath> -ignoreQ 5 
  java -Xmx4g -jar ${params.gatk}/AnalyzeCovariates.jar -l INFO -resources <Rresources> --recal_file <RecalCovariate> -outputDir $StatsAfterFolder> -Rscript <RscriptPath> -ignoreQ 5
  """  
} 
 