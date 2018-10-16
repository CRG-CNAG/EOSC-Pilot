#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Cwd;

sub FASTQCreport{
	my $sampleID=${$_[0]}; my $sampleOutputFolder=${$_[1]}; my $fq1=${$_[2]}; my $fq2=${$_[3]}; my $lane=${$_[4]}; my $fastqcCommand=${$_[5]}; my $fastqcReportScript=${$_[6]};
	my $FastqFilesFolder=${$_[7]};

	my $samplePrefix=$sampleOutputFolder."/".$sampleID."_".$lane;	
	my $log1=$samplePrefix."_1.log";	my $log2=$samplePrefix."_2.log";

	system("$fastqcCommand $fq1 -Dfastqc.output_dir=$sampleOutputFolder -Dfastqc.unzip=false 2>&1 | tee -a $log1");
	my $base1=basename($fq1);   $base1 =~ s{\.[^.]+$}{};	my $Fq1=$FastqFilesFolder."/".$base1."_fastqc";
	my $fqc1summary=$samplePrefix."_1.fastqcsummary.txt";
	my $fqc1summaryLog=$samplePrefix."_1.fastqcsummary.log";
	system("perl $fastqcReportScript -r $Fq1 -p $sampleOutputFolder -o $fqc1summary -l $fqc1summaryLog 2>&1 | tee -a $log1");

	system("$fastqcCommand $fq2 -Dfastqc.output_dir=$sampleOutputFolder -Dfastqc.unzip=false 2>&1 | tee -a $log2");
	my $base2=basename($fq2);   $base2 =~ s{\.[^.]+$}{};	my $Fq2=$FastqFilesFolder."/".$base2."_fastqc";
	my $fqc2summary=$samplePrefix."_2.fastqcsummary.txt";
	my $fqc2summaryLog=$samplePrefix."_2.fastqcsummary.log";
	system("perl $fastqcReportScript -r $Fq2 -p $sampleOutputFolder -o $fqc2summary -l $fqc2summaryLog 2>&1 | tee -a $log2");
}

sub BAMstats{
	my $samplePrefix=${$_[0]}; my $picardToolsCommand=${$_[1]}; my $tmpDir=${$_[2]}; my $logSample=${$_[3]}; my $sortedBam=${$_[4]}; my $fastaFile=${$_[5]};

	my $AlignMetrics=$samplePrefix.".AlignmentSummaryMetrics";
	my $GCBiasMetrics=$samplePrefix.".GcBiasMetrics";                       my $GCBiasMetricsPDF=$GCBiasMetrics.".pdf";
	my $InsertSizeMetrics=$samplePrefix.".CollectInsertSizeMetrics";        my $InsertSizeMetricsPDF=$InsertSizeMetrics.".pdf";
	my $QualityByCycle=$samplePrefix.".MeanQualityByCycle";                 my $QualityByCyclePDF=$QualityByCycle.".pdf";
	my $QualityScoreDist=$samplePrefix.".QualityScoreDistribution";         my $QualityScoreDistPDF=$QualityScoreDist.".pdf";
	my $IndexStat=$samplePrefix.".BamIndexStats";

	system("java -jar -Xmx4g $picardToolsCommand/CollectAlignmentSummaryMetrics.jar I=$sortedBam O=$AlignMetrics R=fastaFile VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");
	system("java -jar -Xmx4g $picardToolsCommand/CollectGcBiasMetrics.jar R=$fastaFile I=$sortedBam O=$GCBiasMetrics CHART=$GCBiasMetricsPDF VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");
	system("java -jar -Xmx4g $picardToolsCommand/CollectInsertSizeMetrics.jar I=$sortedBam O=$InsertSizeMetrics H=$InsertSizeMetricsPDF VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");
	system("java -jar -Xmx4g $picardToolsCommand/MeanQualityByCycle.jar I=$sortedBam O=$QualityByCycle CHART=$QualityByCyclePDF VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");
	system("java -jar -Xmx4g $picardToolsCommand/QualityScoreDistribution.jar I=$sortedBam O=$QualityScoreDist CHART=$QualityScoreDistPDF VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");
	system("java -jar -Xmx4g $picardToolsCommand/BamIndexStats.jar INPUT=$sortedBam VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir > $IndexStat 2>&1 | tee -a $logSample");
}
################################################################################################################
################################################################################################################
my $sampleID=$ARGV[0]; my $lane=$ARGV[1]; my $library=$ARGV[2]; my $platform=$ARGV[3]; my $fq1=$ARGV[4]; my $fq2=$ARGV[5]; my $BAMFolder=$ARGV[6]; my $FastqFilesFolder=$ARGV[7];

my $currentDir=getcwd;
my $ToolsDir=$currentDir."/Tools";
my $BundleDir=$currentDir."/BundleFiles";
my $GenomeFolder=$BundleDir."/GenomeAssembly";
my $fastaFile=$GenomeFolder."/human_g1k_v37.fasta";

my $javaPath=$ToolsDir."/jre1.7.0_80/bin/java";
my $Rscript="/usr/bin/Rscript";
my $bwaCommand=$ToolsDir."/bwa-0.5.8c_patched/bwa";
my $picardToolsCommand=$ToolsDir."/picard-tools-1.32";
my $GATKcommand=$ToolsDir."/dist/GenomeAnalysisTK.jar";
my $AnalyzeCovariateCommand=$ToolsDir."/dist/AnalyzeCovariates.jar";
my $resources=$ToolsDir."/gatk-1.1/public/R";	
my $fastqcCommand=$ToolsDir."/FastQC/fastqc";		
my $fastqcReportScript=$ToolsDir."/fastqc_report_v1.pl";

my $sampleOutputFolder=$BAMFolder."/".$sampleID."_".$lane;     	if(!-d $sampleOutputFolder){    system("mkdir $sampleOutputFolder");}
my $samplePrefix=$sampleOutputFolder."/".$sampleID."_".$lane;

&FASTQCreport(\$sampleID,\$sampleOutputFolder,\$fq1,\$fq2,\$lane,\$fastqcCommand,\$fastqcReportScript,\$FastqFilesFolder);

my $sai1=$samplePrefix."_1.sai";	my $sai2=$samplePrefix."_2.sai";
my $log1=$samplePrefix."_1.log";        my $log2=$samplePrefix."_2.log";
my $cores=8;
system("$bwaCommand aln $fastaFile $fq1 -t $cores -f $sai1 2>&1 | tee -a $log1");
system("$bwaCommand aln $fastaFile $fq2 -t $cores -f $sai2 2>&1 | tee -a $log2");

my $logSample=$samplePrefix.".log";
my $bam=$samplePrefix.".bam";
my $tmpDir=$sampleOutputFolder."/tmp";	if(!-d $tmpDir){	system("mkdir $tmpDir");}
system("$bwaCommand sampe -P -p $platform -i $lane -m $sampleID -l $library $fastaFile $sai1 $sai2 $fq1 $fq2 | java -jar -Xmx4g $picardToolsCommand/SamFormatConverter.jar INPUT=/dev/stdin OUTPUT=$bam VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=2000000 TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");

my $sortedBam=$samplePrefix.".sorted.bam";
my $sortedBamIndex=$sortedBam.".bai";
system("java -jar -Xmx4g $picardToolsCommand/SortSam.jar INPUT=$bam OUTPUT=$sortedBam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");
system("java -jar -Xmx4g $picardToolsCommand/BuildBamIndex.jar INPUT=$sortedBam OUTPUT=$sortedBamIndex VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");
&BAMstats(\$samplePrefix,\$picardToolsCommand,\$tmpDir,\$logSample,\$sortedBam,\$fastaFile);

my $NodupBam=$samplePrefix.".dedup.bam";	my $NodupBamIndex=$NodupBam.".bai";
my $NodupMetrics=$samplePrefix.".dedup.metrics";
system("java -jar -Xmx4g $picardToolsCommand/MarkDuplicates.jar INPUT=$sortedBam OUTPUT=$NodupBam METRICS_FILE=$NodupMetrics REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");
system("java -jar -Xmx4g $picardToolsCommand/BuildBamIndex.jar INPUT=$NodupBam OUTPUT=$NodupBamIndex VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");

my $intervals=$samplePrefix.".intervals";	
my $realignedBam=$samplePrefix.".realigned.bam";
my $Genome1Kindels=$BundleDir."/1000G_phase1.indels.b37.vcf";
my $dbSNP=$BundleDir."/dbsnp_138.b37.excluding_sites_after_129.vcf";
system("$javaPath -Djava.io.tmpdir=$tmpDir -Xmx10g -jar $GATKcommand -l INFO -T RealignerTargetCreator -R $fastaFile -I $NodupBam -o $intervals");
system("$javaPath -Djava.io.tmpdir=$tmpDir -Xmx10g -jar $GATKcommand -l INFO -T IndelRealigner -U ALLOW_UNINDEXED_BAM -I $NodupBam -targetIntervals $intervals -R $fastaFile -known $Genome1Kindels -known $dbSNP -o $realignedBam -LOD 0.4 -compress 0");

my $MateFixedBam=$samplePrefix.".matefixed.bam";
my $MateFixedIndex=$MateFixedBam.".bai";
system("java -jar -Xmx4g $picardToolsCommand/FixMateInformation.jar INPUT=$realignedBam OUTPUT=$MateFixedBam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");
system("java -jar -Xmx4g $picardToolsCommand/BuildBamIndex.jar INPUT=$MateFixedBam OUTPUT=$MateFixedIndex VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");

my $MateFixedCovariate=$samplePrefix.".matefixed.covariate_table.csv";
system("$javaPath -jar -Xmx4g $GATKcommand -l INFO -T CountCovariates -U ALLOW_UNINDEXED_BAM -R $fastaFile -knownSites $dbSNP -I $MateFixedBam -cov ReadGroupcovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile $MateFixedCovariate 2>&1 | tee -a $logSample");

my $recalBam=$samplePrefix.".recal.bam";
system("$javaPath -jar -Xmx4g $GATKcommand -l INFO -T TableRecalibration -U ALLOW_UNINDEXED_BAM -R $fastaFile -I $MateFixedBam --recal_file $MateFixedCovariate --out $recalBam 2>&1 | tee -a $logSample");

my $recalSortedBam=$samplePrefix.".recal.sorted.bam";
my $recalSortedIndex=$recalSortedBam.".bai";
system("java -jar -Xmx4g $picardToolsCommand/SortSam.jar INPUT=$recalBam OUTPUT=$recalSortedBam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");
system("java -jar -Xmx4g $picardToolsCommand/BuildBamIndex.jar INPUT=$recalSortedBam OUTPUT=$recalSortedIndex VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir 2>&1 | tee -a $logSample");

my $RecalCovariate=$samplePrefix.".recal.covariate_table.csv";
system("$javaPath -jar -Xmx4g $GATKcommand -l INFO -T CountCovariates -U ALLOW_UNINDEXED_BAM -R $fastaFile -knownSites $dbSNP -I $recalSortedBam -cov ReadGroupcovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile $RecalCovariate 2>&1 | tee -a $logSample");

my $StatsBeforeFolder=$sampleOutputFolder."/Before";	if(!-d $StatsBeforeFolder){	system("mkdir $StatsBeforeFolder");}
my $StatsAfterFolder=$sampleOutputFolder."/After";	if(!-d $StatsAfterFolder){	system("mkdir $StatsAfterFolder");}
system("$javaPath -jar -Xmx4g $AnalyzeCovariateCommand -l INFO -resources $resources --recal_file $MateFixedCovariate -outputDir $StatsBeforeFolder -Rscript $Rscript -ignoreQ 5 2>&1 | tee -a $logSample");
system("$javaPath -jar -Xmx4g $AnalyzeCovariateCommand -l INFO -resources $resources --recal_file $RecalCovariate -outputDir $StatsAfterFolder -Rscript $Rscript -ignoreQ 5 2>&1 | tee -a $logSample");
