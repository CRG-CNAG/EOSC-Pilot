/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG).
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

params.index = 'data/test2/GonlSamplesToFilesTest.txt'
params.intervals = 'data/test2/Intervals.bed'
params.genref = 'test'  
params.platform = 'illumina'
params.output = 'results'
params.gatk = '/gatk-1.2'
params.R_resources = "/gatk-protected-1.2/public/R" 
params.picard = '/picard-tools-1.32'
params.genome = "${params.output}/genome/$params.genref"

intervals_file = file(params.intervals)

assert params.genref in ['test','v37','v38'], "Invalid `genref` parameter: $params.genref"

/* 
 * download human genome reference file, reference genome indexing and downloading of 1000Genomes ancillary files
 */
process '0_download' {
  tag params.genref
  storeDir params.genome
  
  // make sure there's just file matching the declared extension 
  output:
  file '*.fasta' into gen_fasta_ch
  file '*.vcf' into snp_ch
  file '*.fasta.{bwt,amb,ann,pac,rbwt,rpac,rsa,sa}' into gen_files_ch
  file '*.dict' into dict_ch
  file '*.fasta.fai' into gen_fai_ch 
 
  script:
    if(params.genref=='v38')
    """
    wget -q ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
    wget -q ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi 
    wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai

    gunzip All_20180418.vcf.gz
    < GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai gunzip -d - > GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai
    < GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz gunzip -d - > GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta

    java \
    -XX:ParallelGCThreads=${task.cpus} \
    -jar ${params.picard}/picard.jar CreateSequenceDictionary R=GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta O=GCA_000001405.15_GRCh38_no_alt_analysis_set

    bwa index -a bwtsw GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
    """

    else if(params.genref=='v37') 
    '''
    wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.gz
    wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.fai.gz
    wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.dict.gz
    wget -q ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.excluding_sites_after_129.vcf.gz

    gunzip dbsnp_138.b37.excluding_sites_after_129.vcf.gz
    gunzip human_g1k_v37.dict.gz
    gunzip human_g1k_v37.fasta.fai.gz
    gunzip human_g1k_v37.fasta.gz
    
    bwa index -a bwtsw human_g1k_v37.fasta 
    '''

    else
    """
    cp $baseDir/data/test/human_g1k_v37.fasta.gz human_g1k_v37.fasta.gz
    cp $baseDir/data/test/human_g1k_v37.fasta.fai.gz human_g1k_v37.fasta.fai.gz 
    cp $baseDir/data/test/human_g1k_v37.dict.gz human_g1k_v37.dict.gz
    cp $baseDir/data/test/dbsnp_138.b37.excluding_sites_after_129.vcf.gz dbsnp_138.b37.excluding_sites_after_129.vcf.gz 

    gunzip dbsnp_138.b37.excluding_sites_after_129.vcf.gz
    gunzip human_g1k_v37.dict.gz
    gunzip human_g1k_v37.fasta.fai.gz
    gunzip human_g1k_v37.fasta.gz
    
    bwa index -a bwtsw human_g1k_v37.fasta 
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


process '5_bam_recalibration' {
  tag "${prefixId}"
  publishDir params.output 
  input:
    file gen_fasta from gen_fasta_ch
    file gen_fai from gen_fai_ch
    file snp_file from snp_ch
    file dict_file from dict_ch 
    set val(prefixId), file(recalibrated_table), file(dedup_bam) from recalibrate_dedup_bam_ch

  output:
  file '*.recalibrated.bam'
  """
  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R $gen_fasta \
    -I $dedup_bam \
    -knownSites $snp_file \
    -fixMisencodedQuals \
    -o recal_data.table
 
  java \
    -XX:ParallelGCThreads=${task.cpus} \
    -Xmx${task.memory?.giga?:1}g \
    -jar ${params.gatk}/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R $gen_fasta \
    -I $dedup_bam \
    -BQSR recal_data.table \
    -o ${prefixId}.recalibrated.bam
   """
}