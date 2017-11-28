/*
 * Copyright (c) 2017, Centre for Genomic Regulation (CRG).
 *
 *   This file is part of 'EOSC-Pilot': 
 *   A Nextflow pipeline for Variant Calling with NGS data
 *
 *   EOSC-Pilot is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   EOSC-Pilot is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with EOSC-Pilot.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * EOSC-Pilot project 
 * 
 * Authors:
 * - Nino Spataro <nino.spataro@crg.eu>
 * - Paolo Di Tommaso <paolo.ditommaso@crg.eu>
 */

params.index = 'data/test/fastq/GonlSamplesToFilesTest.txt'
params.test = false
params.platform = 'illumina'
params.output = 'results'
params.gatk = '/gatk-1.2'
params.R_resources = "/gatk-protected-1.2/public/R" 
params.picard = '/picard-tools-1.32'
params.genome = "${params.output}/genome"


/* 
 * download human genome reference file, reference genome indexing and downloading of 1000Genomes ancillary files
 */
process '0_download' {
  storeDir params.genome
  
  output:
  file 'human_g1k_v37.fasta' into gen_fasta_ch
  file '1000G_phase1.indels.b37.vcf' into indels_ch 
  file 'dbsnp_138.b37.excluding_sites_after_129.vcf' into snp_ch
  file 'human_g1k_v37.fasta.{bwt,amb,ann,pac,rbwt,rpac,rsa,sa}' into gen_files_ch
  file 'human_g1k_v37.dict' into dict_ch
  file 'human_g1k_v37.fasta.fai' into gen_fai_ch 
 
  script:
    """
    ${( !params.test ? 
      '''
    wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.gz
    wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.fai.gz
    wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.dict.gz
    wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.excluding_sites_after_129.vcf.gz                
    wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz
    wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_07/low_coverage/snps/CEU.low_coverage.2010_07.genotypes.vcf.gz
      '''
      : 
      """
    cp $baseDir/data/test/human_g1k_v37.fasta.gz human_g1k_v37.fasta.gz
    cp $baseDir/data/test/human_g1k_v37.fasta.fai.gz human_g1k_v37.fasta.fai.gz 
    cp $baseDir/data/test/human_g1k_v37.dict.gz human_g1k_v37.dict.gz
    cp $baseDir/data/test/1000G_phase1.indels.b37.vcf.gz 1000G_phase1.indels.b37.vcf.gz
    cp $baseDir/data/test/dbsnp_138.b37.excluding_sites_after_129.vcf.gz dbsnp_138.b37.excluding_sites_after_129.vcf.gz 
    cp $baseDir/data/test/CEU.low_coverage.2010_07.genotypes.vcf.gz CEU.low_coverage.2010_07.genotypes.vcf.gz
      """	
    )}

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
 * parse the index file and extract relevant metadata such as: 
 * - lane 
 * - sampleID 
 * - library
 * - prefixID
 * - read pair files
 */

READ_FILENAME_REGEX = /^([\w]+_[\w]+_[\w]+_([\w]+)_([\w]+))_[12]$/

Channel.fromPath(params.index)
       .splitCsv(sep:'\t', skip:1)
       .map{ sampleId, fileId, fullPath -> 
            def readFile = file(fullPath)
            def name = readFile.simpleName 
            def regex = (name =~ READ_FILENAME_REGEX) 
            if( !regex.matches() ) 
                error "Invalid read pair file name format: $fullPath" 
          
            def baseName = regex.group(1)
            def meta = [:]
            meta.lane = regex.group(2)
            meta.library = regex.group(3)
            meta.sampleId = sampleId
            meta.prefixId = "${sampleId}_${meta.lane}"
  
            tuple(baseName, meta, readFile)
            }
       .groupTuple(size:2)
       .map { base, metas, files -> 
            assert metas[0]==metas[1]; 
            files.sort(); 
            tuple(metas[0], files[0], files[1]) 
            }
       .into { reads_ch1; reads_ch2; reads_ch3 }   

/* 
 * for each sample lane the two fastq files are processed to get the the BAM recalibrated files. All the following tasks are performed in the script "FromFastqToBam.pl".
 * 
 * Quality control of the fastq files
 */
process '1_quality_control' {
  tag "${meta.prefixId}"
  
  input: 
  set meta, file(read_1), file(read_2) from reads_ch1
 
  output:
  file 'fqc{1,2}summary.{txt,log}'
  file 'sample_out'
  
  script:
  """
  mkdir sample_out
  fastqc $read_1 -Dfastqc.output_dir=sample_out -Dfastqc.unzip=false 
  fastqc_report_v1.pl -r $read_1 -p sample_out -o fqc1summary.txt -l fqc1summary.log

  fastqc $read_2 -Dfastqc.output_dir=sample_out -Dfastqc.unzip=false 
  fastqc_report_v1.pl -r $read_2 -p sample_out -o fqc2summary.txt -l fqc2summary.log 
  """
}


process '2_create_sai_files' {
  tag "${meta.prefixId}"
  
  input: 
  file gen_fasta from gen_fasta_ch
  file bwt_file from gen_files_ch 
  set meta, file(read_1), file(read_2) from reads_ch2

  output:
  set val(meta), file('*_1.sai'), file('*_2.sai') into sai_ch
  
  script:
  """
  bwa aln $gen_fasta $read_1 -t $task.cpus -f "${meta.prefixId}_1.sai"
  bwa aln $gen_fasta $read_2 -t $task.cpus -f "${meta.prefixId}_2.sai" 
  """
}

process '3_align_to_genome' {
  tag "${meta.prefixId}"
  
  input:
  file gen_fasta from gen_fasta_ch
  file bwt_files from gen_files_ch
  set meta, file(read_1), file(read_2), \
            file(sai1), file(sai2) from reads_ch3.join(sai_ch)
  
  output: 
  set val(meta), file('*.bam') into bam_ch 
  
  script:
  """
  bwa sampe -P -p ${params.platform} -i ${meta.lane} -m ${meta.sampleId} -l ${meta.library} $gen_fasta $sai1 $sai2 $read_1 $read_2 | \\
  java -Xmx4g -jar ${params.picard}/SamFormatConverter.jar INPUT=/dev/stdin OUTPUT="${meta.prefixId}.bam" VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=2000000 TMP_DIR=\$TMPDIR 
  """
}

/* 
 * sorting and indexing of the bam file generated in step 3 
 */
process '4_sort_and_index' {
  tag "${meta.prefixId}"
  
  input: 
  set meta, file(bam_file) from bam_ch 
  
  output: 
  set meta, file('*.sorted.bam'), file('*.sorted.bam.bai') into sorted_ch
  
  script:
  """
  java -Xmx4g -jar ${params.picard}/SortSam.jar INPUT=$bam_file OUTPUT=${meta.prefixId}.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=${meta.prefixId}.sorted.bam OUTPUT=${meta.prefixId}.sorted.bam.bai VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR 
  """
}

/* 
 * removing of optical duplicated from the sorted bam file and subsequent indexing of the duplicates free bam file
 */
process '5_dedup_and_index' {
  tag "${meta.prefixId}"
  
  input: 
  set meta, file(sorted_bam), file(sorted_bai) from sorted_ch
  
  output: 
  set meta, file ('*.dedup.bam'), file ('*.dedup.bam.bai') into dedup_ch
  
  script:
  """
  java -Xmx4g -jar ${params.picard}/MarkDuplicates.jar INPUT=$sorted_bam OUTPUT=${meta.prefixId}.dedup.bam METRICS_FILE=${meta.prefixId}.dedup.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=\$TMPDIR
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=${meta.prefixId}.dedup.bam OUTPUT=${meta.prefixId}.dedup.bam.bai VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR 
  """

} 

/* 
 *  realign reads around indels
 */
process '6_realign_indels' {
  tag "${meta.prefixId}"
  
  input:
  file gen_file from gen_fasta_ch
  file gen_fai from gen_fai_ch
  file indels from indels_ch
  file snp_file from snp_ch
  file dict_file from dict_ch
  set meta, file(dedup_bam), file(dedup_bai) from dedup_ch

  output: 
  set meta, file('*.realigned.bam') into realigned_ch 
  
  script:
  """
  java -Xmx10g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T RealignerTargetCreator -R $gen_file -I $dedup_bam -o ${meta.prefixId}.intervals
  java -Xmx10g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T IndelRealigner -U ALLOW_UNINDEXED_BAM -I $dedup_bam -targetIntervals ${meta.prefixId}.intervals -R $gen_file -known $indels -known $snp_file -o ${meta.prefixId}.realigned.bam -LOD 0.4 -compress 0
  """ 	

} 

/* 
 * fixing mate reads and indexing of the fixed mate bam file
 */
process '7_fixing_and_indexing' {
  tag "${meta.prefixId}"
  
  input:
  set meta, file(realigned_bam) from realigned_ch 
  
  output: 
  set meta, file('*.matefixed.bam'), file('*.matefixed.bam.bai') into matefixed_ch
  
  script:
  """
  java -Xmx4g -jar ${params.picard}/FixMateInformation.jar INPUT=$realigned_bam OUTPUT=${meta.prefixId}.matefixed.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=\$TMPDIR 
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=${meta.prefixId}.matefixed.bam OUTPUT=${meta.prefixId}.matefixed.bam.bai VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=\$TMPDIR 
  """
}

/* 
 * bam recalibration and subsequent sorting and indexing
 */
process '8_recalibrate_and_sort' {
  tag "${meta.prefixId}"
  publishDir params.output
  
  input:
  file gen_file from gen_fasta_ch
  file snp_file from snp_ch
  file dict_file from dict_ch
  file gen_fai from gen_fai_ch
  set meta, file(matefixed_bam), file(matefixed_bai) from matefixed_ch
  
  output:
  set meta, file('*.recal.sorted.bam') into recal_ch
  set meta, file('*.matefixed.covariate_table.csv') into matefixed_cov_ch
  file '*.recal.sorted.bam'
  script:
  """
  java -Xmx4g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T CountCovariates -U ALLOW_UNINDEXED_BAM -R $gen_file -knownSites $snp_file -I $matefixed_bam -cov ReadGroupcovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile ${meta.prefixId}.matefixed.covariate_table.csv 
  java -Xmx4g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T TableRecalibration -U ALLOW_UNINDEXED_BAM -R $gen_file -I $matefixed_bam --recal_file *.matefixed.covariate_table.csv --out  ${meta.prefixId}.recal.bam
  java -Xmx4g -jar ${params.picard}/SortSam.jar INPUT=${meta.prefixId}.recal.bam OUTPUT=${meta.prefixId}.recal.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=\$TMPDIR 
  java -Xmx4g -jar ${params.picard}/BuildBamIndex.jar INPUT=${meta.prefixId}.recal.sorted.bam OUTPUT=${meta.prefixId}.recal.sorted.bam.bai VALIDATION_STRINGENCY=LENIENT TMP_DIR=\$TMPDIR 
  """
}

/* 
 * 	analysis of the recalibration process through comparison of some metrics between MateFixedBam (before recalibration) and recalSortedBam (after recalibration)
 */  
process '9_recalibrate_and_compare' {
  tag "${meta.prefixId}"
  
  input:
  file gen_file from gen_fasta_ch
  file snp_file from snp_ch
  file dict_file from dict_ch
  file gen_fai from gen_fai_ch
  set meta, file(recal_bam) from recal_ch
  set meta, file(matefixed_cov) from matefixed_cov_ch

  script:
  """
  mkdir Before
  mkdir After 
  java -Xmx4g -jar ${params.gatk}/GenomeAnalysisTK.jar -l INFO -T CountCovariates -U ALLOW_UNINDEXED_BAM -R $gen_file -knownSites $snp_file -I $recal_bam -cov ReadGroupcovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile ${meta.prefixId}.recal.covariate_table.csv
  java -Xmx4g -jar ${params.gatk}/AnalyzeCovariates.jar -l INFO -resources ${params.R_resources} --recal_file $matefixed_cov -outputDir Before -Rscript `which R` -ignoreQ 5 
  java -Xmx4g -jar ${params.gatk}/AnalyzeCovariates.jar -l INFO -resources ${params.R_resources} --recal_file ${meta.prefixId}.recal.covariate_table.csv -outputDir After -Rscript `which R` -ignoreQ 5
  """  
} 
 