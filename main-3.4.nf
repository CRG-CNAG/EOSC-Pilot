/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 * 
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */

/* 
 * EOSC-Pilot project 
 * 
 * Authors:
 * - Nino Spataro <nino.spataro@crg.eu>
 * - Paolo Di Tommaso <paolo.ditommaso@crg.eu>
 */

params.index = 'data/test2/GonlSamplesToFilesTest.txt'
params.intervals = 'data/test2/Intervals.bed'
params.test = false
params.platform = 'illumina'
params.output = 'results'
params.gatk = '/gatk-1.2'
params.R_resources = "/gatk-protected-1.2/public/R" 
params.picard = '/picard-tools-1.32'
params.genome = "${params.output}/genome"

intervals_file = file(params.intervals)

/* 
 * download human genome reference file, reference genome indexing and downloading of 1000Genomes ancillary files
 */
process '0_download' {
  storeDir params.genome
  
  output:
  file 'human_g1k_v37.fasta' into gen_fasta_ch
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

READ_FILENAME_REGEX = /^(([\w]+_[\w]+_[\w]+)_([\w]+)_([\w]+))_[12]$/

Channel.fromPath(params.index)
       .splitCsv(sep:'\t', skip:1)
       .tap { records_ch }
       .map{ sampleId, fileId, fullPath -> 
            def readFile = file(fullPath)
            def name = readFile.simpleName 
            def regex = (name =~ READ_FILENAME_REGEX) 
            if( !regex.matches() ) 
                error "Invalid read pair file name format: $fullPath" 
          
            def baseName = regex.group(1)
            def meta = [:]
            meta.lane = regex.group(3)
            meta.library = regex.group(4)
            meta.sampleId = sampleId
            meta.prefixId = "${sampleId}_${regex.group(2)}_${regex.group(3)}"
            meta.familyId = sampleId[0..-2]
  
            tuple(baseName, meta, readFile)
            }
       .groupTuple(size:2)
       .map { base, metas, files -> 
            assert metas[0]==metas[1]; 
            files.sort(); 
            tuple(metas[0], files[0], files[1]) 
            }
       .into { reads_ch1; reads_ch2 }   

/* 
 * create a PED file for each family ID
 */
records_ch
  .map { sampleId, fileId, fullPath -> sampleId }
  .unique()
  .map { sampleId -> 
      def familyId = sampleId[0..-2]
      if( sampleId.endsWith('a') ) 
        return "$familyId\t$sampleId\t0\t0\t1\t-9"
      if( sampleId.endsWith('b') ) 
        return "$familyId\t$sampleId\t0\t0\t2\t-9"
      else 
        return "$familyId\t$sampleId\t${familyId}a\t${familyId}b\t-9\t-9"   
  }
  .collectFile(newLine:true) { line -> 
    def items = line.tokenize('\t');
    tuple("${items[0]}.ped", line) 
  }
  .map { file -> tuple(file.simpleName, file) }
  .into { ped_files_ch; ped_files_ch2 }

/* 
 * for each sample lane the two fastq files are processed to get the the BAM recalibrated files. All the following tasks are performed in the script "FromFastqToBam.pl".
 * 
 * Quality control of the fastq files
 */
process '1_quality_control' {
  tag "${meta.prefixId}"
  
  input: 
  set meta, file(fq1), file(fq2) from reads_ch1
 
  output:
  file 'sample_out'
  
  script:
  """
  mkdir sample_out
  fastqc $fq1 $fq2 -o sample_out --noextract 
  """
}

process '2_align_and_sort' {
  tag "${meta.prefixId}"
  
  input: 
  file gen_fasta from gen_fasta_ch
  file gen_index from gen_files_ch 
  set meta, file(fq1), file(fq2) from reads_ch2

  output:
  set val(meta.sampleId), file('*.sorted.bam') into sorted_bam_ch

  """
  bwa mem -M -R "@RG\\tID:${meta.prefixId}\\tPL:illumina\\tLB:${meta.lane}\\tSM:${meta.sampleId}" -t ${task.cpus} $gen_fasta $fq1 $fq2 > temp.sam  
  grep -v '@PG' temp.sam > aligned.sam
  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -jar ${params.picard}/picard.jar \
    SortSam INPUT=aligned.sam \
    OUTPUT=${meta.prefixId}.sorted.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true
  """
}

process '3_merge_bam' {
  tag "${prefixId}"

  input:
    set val(prefixId), file('*') from sorted_bam_ch.groupTuple()
  output: 
    set val(prefixId), file('*.bam') into merged_bam_ch

  script:
  """
  sambamba merge ${prefixId}.merged.bam *.bam
  """
}

process '4_dedup_merge' {
  tag "${prefixId}"

  input:
    file gen_fasta from gen_fasta_ch
    file gen_fai from gen_fai_ch
    file snp_file from snp_ch
    file dict_file from dict_ch
    set val(prefixId), file(merged_bam) from merged_bam_ch  
  output:
    set val(prefixId), file('*.recalibration.csv'), file('*.merged.dedup.bam') into recalibrate_dedup_bam_ch
  
  script:
  """
  sambamba markdup --nthreads=${task.cpus} --overflow-list-size 1000000 --hash-table-size 1000000 -p $merged_bam ${prefixId}.merged.dedup.bam
  java \
    -XX:ParallelGCThreads=${task.cpus} -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R $gen_fasta \
    -nct 8 \
    -knownSites $snp_file \
    -I $merged_bam \
    -o ${prefixId}.recalibration.csv \
    -allowPotentiallyMisencodedQuals -fixMisencodedQuals
  """
}


process '5_varcall' {
  tag "${prefixId}"

  input:
    file gen_fasta from gen_fasta_ch
    file gen_fai from gen_fai_ch
    file snp_file from snp_ch
    file dict_file from dict_ch 
    set val(prefixId), file(recalibrated_table), file(dedup_bam) from recalibrate_dedup_bam_ch

  output:
    set prefixId, file('*.g.vcf') into vcf_file_ch

  script:
  """
  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R $gen_fasta \
    -I $dedup_bam \
    -newQual \
    --BQSR $recalibrated_table \
    --dbsnp $snp_file \
    -o ${prefixId}.g.vcf \
    --emitRefConfidence GVCF \
    -ploidy 2 \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000
  """
}

/* 
 * group together all gvcf files within the same family 
 * the family-id is given by the prefix-id from which is remved the last 
 * character (eg. `gonl-25b` -> `gonl-25`)
 */ 

vcf_file_ch
  .map { prefixId, files -> tuple(prefixId[0..-2], files) }
  .groupTuple()
  .tap { family_vcf_ch }
  .join( ped_files_ch )
  .set { vcf_and_ped_files_ch }   


process '6_combine_vcf' {
  tag "${familyId}"

  input: 
    file gen_fasta from gen_fasta_ch
    file gen_fai from gen_fai_ch
    file snp_file from snp_ch
    file dict_file from dict_ch  
    file intervals_file 
    set familyId, file(vcf_files) from family_vcf_ch

  output: 
    set familyId, file('*.g.vcf') into combined_vcf_ch

  script:
  def gvcfs = vcf_files.collect { "--variant $it" }.join(' ')

  """
  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T CombineGVCFs \
    -R $gen_fasta \
    --dbsnp $snp_file $gvcfs \
    -o ${familyId}.g.vcf
  """
}

process '7_merge_family' {
  tag "${familyId}"

  input:
    file gen_fasta from gen_fasta_ch
    file gen_fai from gen_fai_ch
    file snp_file from snp_ch
    file dict_file from dict_ch   
    set familyId, file(vcf_files), file(ped_file) from vcf_and_ped_files_ch

  output:
    set familyId, file('*.vcf') into family_merge_vcf_ch

  // TODO 
  // 1. add the cli option `-L $intervals_file` to the GenotypeGVCFs command 
  script:
  def names = vcf_files.collect { "--variant $it" }.join(' ')
  """
  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R $gen_fasta \
    --dbsnp $snp_file \
    -o ${familyId}.vcf \
    $names
  """
}

process '8_phase_by_transmission' {
  tag "${familyId}"

  input:
    file gen_fasta from gen_fasta_ch
    file gen_fai from gen_fai_ch
    file dict_file from dict_ch
    set familyId, file(family_vcf), file(ped_file) from family_merge_vcf_ch.join(ped_files_ch2)
  output: 
    file('*.PBT.vcf') into pbt_files_ch
  script:
  """
  awk '\$2!="${familyId}d"' $ped_file > abc.ped
  awk '\$2!="${familyId}c"' $ped_file > abd.ped
  
  # process `abc.ped` file 
  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T PhaseByTransmission \
    -R $gen_fasta \
    -V $family_vcf \
    --pedigreeValidationType SILENT \
    -mvf ${familyId}.ABC.mvf \
    -ped abc.ped \
    -o ${familyId}.abc.vcf 

  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $gen_fasta \
    --variant ${familyId}.abc.vcf \
    -o ${familyId}.abc.filtered.vcf \
    -xl_se \\'gonl-.+d\\'  
 
  if [[ `grep -c ${familyId} abd.ped ` != 3 ]]; then
      cp ${familyId}.abc.filtered.vcf ${familyId}.PBT.vcf
      exit 0 
  fi     

  # process `abd.ped` file 
  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T PhaseByTransmission \
    -R $gen_fasta \
    -V $family_vcf \
    --pedigreeValidationType SILENT \
    -mvf ${familyId}.ABD.mvf \
    -ped abd.ped \
  -o ${familyId}.abd.vcf

  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $gen_fasta \
    --variant ${familyId}.abd.vcf \
    -o ${familyId}.abd.filtered.vcf

  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R $gen_fasta \
    --variant:'' *.abc.filtered.vcf --variant:'' *.abd.filtered.vcf \
    -o ${familyId}.PBT.vcf \
    -genotypeMergeOptions UNIQUIFY

    correction.pl ${familyId}.PBT.vcf 
  """
}

process '9_merge_pbt' {
  input:
    file gen_fasta from gen_fasta_ch
    file gen_fai from gen_fai_ch
    file dict_file from dict_ch
    file (vcf_files) from pbt_files_ch.collate(50).collect()

  output:
    file 'sorted.vcf' into (merged_vcf_ch1, merged_vcf_ch2)

  script:
  def names = vcf_files.collect { "--variant $it" }.join(' ')
  """
  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R $gen_fasta \
    $names \
    -o merged.vcf \
    -genotypeMergeOptions UNIQUIFY

    correction.pl merged.vcf  

    sortVCFbyFai.pl \
      -fastaIndexFile $gen_fai \
      -inputVCF merged.vcf \
      -outputVcf sorted.vcf
  """
}

process '10a_snps_filtering' {
  input:
    file gen_fasta from gen_fasta_ch
    file gen_fai from gen_fai_ch
    file dict_file from dict_ch
    file joined_vcf_file from merged_vcf_ch1
  output:
    file 'snps.filtered.vcf' into filtered_snps_ch
    
  script:
  """
  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $gen_fasta \
    --variant $joined_vcf_file \
    -o snp.vcf \
    --selectTypeToExclude INDEL

  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $gen_fasta \
    --variant snp.vcf \
    -o snps.filtered.vcf \
    --filterExpression 'QD < 2.0' \
    --filterName 'filterQD' \
    --filterExpression 'MQ < 25.0' \
    --filterName 'filterMQ' \
    --filterExpression 'FS > 60.0' \
    --filterName 'filterFS' \
    --filterExpression 'MQRankSum < -12.5' \
    --filterName 'filterMQRankSum' \
    --filterExpression 'ReadPosRankSum < -8.0' \
    --filterName 'filterReadPosRankSum'
  """
}

process '10b_indels_filtering' {
  input:
    file gen_fasta from gen_fasta_ch
    file gen_fai from gen_fai_ch
    file dict_file from dict_ch
    file joined_vcf_file from merged_vcf_ch2
  output:
    file 'indels.filtered.vcf' into filtered_indels_ch
  script:
  """
  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $gen_fasta \
    --variant $joined_vcf_file \
    -o indels.vcf \
    --selectTypeToInclude INDEL

  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $gen_fasta \
    --variant indels.vcf \
    -o indels.filtered.vcf \
    --filterExpression 'QD < 2.0' \
    --filterName 'filterQD' \
    --filterExpression 'FS > 200.0' \
    --filterName 'filterFS' \
    --filterExpression 'ReadPosRankSum < -20.0' \
    --filterName 'filterReadPosRankSum'
  """
}

process '11_combine_filtered' {
  publishDir params.output
  input:
    file gen_fasta from gen_fasta_ch
    file gen_fai from gen_fai_ch
    file dict_file from dict_ch
    file indels_file from filtered_indels_ch
    file snps_file from filtered_snps_ch  
  output:
    file 'joined.filtered.vcf.gz.*'

  script:
  """
    java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R $gen_fasta \
    --variant $snps_file \
    --variant $indels_file \
    --genotypemergeoption UNSORTED \
    -o joined.filtered.vcf

  bgzip -c joined.filtered.vcf > joined.filtered.vcf.gz
  tabix -p vcf joined.filtered.vcf.gz
  """
}