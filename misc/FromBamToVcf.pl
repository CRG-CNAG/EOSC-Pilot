#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;    
use Cwd;

sub BAMstats{
        my $samplePrefix=${$_[0]}; my $picardToolsCommand=${$_[1]}; my $bam=${$_[2]}; my $tmpDir=${$_[3]}; my $fastaFile=${$_[4]};

        my $AlignMetrics=$samplePrefix.".AlignmentSummaryMetrics";
        my $GCBiasMetrics=$samplePrefix.".GcBiasMetrics";                       my $GCBiasMetricsPDF=$GCBiasMetrics.".pdf";
        my $InsertSizeMetrics=$samplePrefix.".CollectInsertSizeMetrics";        my $InsertSizeMetricsPDF=$InsertSizeMetrics.".pdf";
        my $QualityByCycle=$samplePrefix.".MeanQualityByCycle";                 my $QualityByCyclePDF=$QualityByCycle.".pdf";
        my $QualityScoreDist=$samplePrefix.".QualityScoreDistribution";         my $QualityScoreDistPDF=$QualityScoreDist.".pdf";
        my $IndexStat=$samplePrefix.".BamIndexStats";

        system("java -jar -Xmx4g $picardToolsCommand/CollectAlignmentSummaryMetrics.jar I=$bam O=$AlignMetrics R=$fastaFile VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir");
        system("java -jar -Xmx4g $picardToolsCommand/CollectGcBiasMetrics.jar R=$fastaFile I=$bam O=$GCBiasMetrics CHART=$GCBiasMetricsPDF VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir");
        system("java -jar -Xmx4g $picardToolsCommand/CollectInsertSizeMetrics.jar I=$bam O=$InsertSizeMetrics H=$InsertSizeMetricsPDF VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir");
        system("java -jar -Xmx4g $picardToolsCommand/MeanQualityByCycle.jar I=$bam O=$QualityByCycle CHART=$QualityByCyclePDF VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir");
        system("java -jar -Xmx4g $picardToolsCommand/QualityScoreDistribution.jar I=$bam O=$QualityScoreDist CHART=$QualityScoreDistPDF VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir");
        system("java -jar -Xmx4g $picardToolsCommand/BamIndexStats.jar INPUT=$bam VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmpDir");
}


################################################################################################################
################################################################################################################    
my $sample=$ARGV[0]; my $VariantCallingFolder=$ARGV[1]; my @bams=@ARGV[2..$#ARGV];

my $currentDir=getcwd;

my $tmpDir=$VariantCallingFolder."/tmp";        if(!-d $tmpDir){        system("mkdir $tmpDir");}
my $ToolsDir=$currentDir."/Tools";
my $BundleDir=$currentDir."/BundleFiles";
my $GenomeFolder=$BundleDir."/GenomeAssembly";
my $fastaFile=$GenomeFolder."/human_g1k_v37.fasta";

my $javaPath=$ToolsDir."/jre1.7.0_80/bin/java";
my $Rscript="/usr/bin/Rscript";
my $picardToolsCommand=$ToolsDir."/picard-tools-1.32";
my $GATKcommand=$ToolsDir."/dist/GenomeAnalysisTK.jar";
my $evalRscript=$ToolsDir."/extract_info_GATK_variantEval.R";            
my $plot_coverageScript=$ToolsDir."/plot_coverage-1.1.R";
my $plot_cumulativeScript=$ToolsDir."/plot_cumulative_coverage-1.1.R";
#my $verifyBAMidScript=$ToolsDir."/"; ##this software is still not installed, we are facing some problems but it should just be related to a final quality control and should not affect data production

my $dbSNP=$GenomeFolder."/dbsnp_138.b37.excluding_sites_after_129.vcf";
my $NoAmbiguous=$BundleDir."/CEU.low_coverage.2010_07.genotypes.vcf";	##my $NoAmbiguous=$BundleDir."/GoNL_Immuno.noAmbiguous.vcf"; ###
#my $NoAmbiguousBed=$BundleDir."/GoNL_Immuno.noAmbiguous.bed";
#my $NoAmbiguousBim=""; ###
#my $NoAmbiguousFam=""; ###
#my $IntervalList=""; ###
#####All the commented lines are referring to ancillary files required for the variant calling. When possible these files were substituted by 1000Genomes Project files,
#####but to produce the same results the original files are required. Still waiting for these files from the GONL team

my $inputBAMstring=join(" ",@bams);
my $SampleVariantCallingFolder=$VariantCallingFolder."/".$sample;	if(!-d $SampleVariantCallingFolder){	system("mkdir $SampleVariantCallingFolder");}
my $samplePrefix=$SampleVariantCallingFolder."/".$sample;

my $mergedBam=$samplePrefix.".merged.bam";
my $mergedIndex=$mergedBam.".bai";        
system("java -jar -Xmx4g $picardToolsCommand/MergeSamFiles.jar $inputBAMstring ASSUME_SORTED=true USE_THREADING=true TMP_DIR=$tmpDir OUTPUT=$mergedBam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT");        
system("java -jar -Xmx4g $picardToolsCommand/BuildBamIndex.jar INPUT=$mergedBam OUTPUT=$mergedIndex VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 TMP_DIR=$tmpDir");
        
my $CheckSnps=$samplePrefix.".merged.qc_check_snps.vcf";        
my $CheckMetrics=$samplePrefix.".merged.qc_check_snps.vcf.metrics";        
system("$javaPath -Xmx10g -Djava.io.tmpdir=$tmpDir -jar $GATKcommand -l INFO -T UnifiedGenotyper -I $mergedBam --out $CheckSnps -R $fastaFile -D $dbSNP -stand_call_conf 40.0 -stand_emit_conf 10.0 -dcov 200 --metrics_file $CheckMetrics -L $IntervalList");
        
my $Eval=$samplePrefix.".merged.qc_check_snps.concordance.eval";        
system("$javaPath -Xmx10g -Djava.io.tmpdir=$tmpDir -jar $GATKcommand -T VariantEval --eval $CheckSnps --comp $NoAmbiguous -o $Eval -R $fastaFile -D $dbSNP -EV GenotypeConcordance");
        
my $Concordance=$samplePrefix.".merged.qc_check_snps.concordance.txt";        
system("Rscript $evalRscript --in $Eval --step concordance --name $sample --comp comp_immuno --header >> $Concordance");
        
&BAMstats(\$samplePrefix,\$picardToolsCommand,\$mergedBam,\$tmpDir,\$fastaFile);        
my $coverage=$samplePrefix.".coverage";        
system("$javaPath -Djava.io.tmpdir=$tmpDir -Xmx12g -jar $GATKcommand -T DepthOfCoverage -R $fastaFile -I $mergedBam -o $coverage --omitDepthOutputAtEachBase -omitIntervals -omitSampleSummary -nt 4");
        
my $coverageCumulativeCount=$samplePrefix.".coverage.sample_cumulative_coverage_counts";        
my $coveragePDF=$samplePrefix.".coverage.pdf";        
system("Rscript $plot_coverageScript --in $coverageCumulativeCount --out $coveragePDF --expected-coverage 12 --max-depth 40 --title \"Coverage $sample\"");
        
my $coverageCumulativeProportion=$samplePrefix.".coverage.sample_cumulative_coverage_proportions";        
my $coverageCumulativePDF=$samplePrefix.".cumulative_coverage.pdf";        
system("Rscript $plot_cumulativeScript --in $coverageCumulativeProportion --out $coverageCumulativePDF --expected-coverage 12 --max-depth 40 --title \"Coverage $sample\"");
        
#system("$verifyBAMidScript --reference $fastaFile --in $mergedBam --bfile GoNL_Immuno.noAmbiguous --out $samplePrefix --maxDepth 30 --verbose");
my $bestMerged=$samplePrefix.".best.merged";    my $bestRG=$samplePrefix.".bestRG";     my $bestSM=$samplePrefix.".bestSM";
#system("cat $bestRG $bestSM $bestMerged");
        
my $selfMerged=$samplePrefix.".self.merged";    my $selfRG=$samplePrefix.".selfRG";     my $selfSM=$samplePrefix.".selfSM";        
#system("cat $selfRG $selfSM $selfMerged");

