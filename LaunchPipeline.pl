#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Cwd;

sub retrieveGenomeData{
        my $GenomeFolder=${$_[0]}; my $ReferenceToRetrieve=${$_[1]}; my $bwaCommand=${$_[2]}; my $dbSNPToRetrieve=${$_[3]}; my $Genome1KindelsToRetrieve=${$_[4]};

        my $fastaToRetrieve=$ReferenceToRetrieve.".fasta.gz";      my $fastaFile=$GenomeFolder."/".$fastaToRetrieve;	my $checkFasta=$GenomeFolder."/".$ReferenceToRetrieve.".fasta";
        my $faiToRetrieve=$ReferenceToRetrieve.".fasta.fai.gz";    my $faiFile=$GenomeFolder."/".$faiToRetrieve;	my $checkFai=$GenomeFolder."/".$ReferenceToRetrieve.".fasta.fai";
        my $dictToRetrieve=$ReferenceToRetrieve.".dict.gz";        my $dictFile=$GenomeFolder."/".$dictToRetrieve;	my $checkDict=$GenomeFolder."/".$ReferenceToRetrieve.".dict";

        if(!-e $checkFasta){    system("wget ftp://gsapubftp-anonymous\@ftp.broadinstitute.org/bundle/b37/$fastaToRetrieve -O $fastaFile"); 	system("gunzip $fastaFile");}
        if(!-e $checkFai){      system("wget ftp://gsapubftp-anonymous\@ftp.broadinstitute.org/bundle/b37/$faiToRetrieve -O $faiFile");		system("gunzip $faiFile");}
        if(!-e $checkDict){     system("wget ftp://gsapubftp-anonymous\@ftp.broadinstitute.org/bundle/b37/$dictToRetrieve -O $dictFile");	system("gunzip $dictFile");}

        my $bwtFile=$checkFasta.".bwt";
        if(!-e $bwtFile){      system("$bwaCommand index -a bwtsw $checkFasta");}

	$dbSNPToRetrieve.=".gz";
	$Genome1KindelsToRetrieve.=".gz";
        my $dbSNP=$GenomeFolder."/".$dbSNPToRetrieve;				my $checkDBsnp=$GenomeFolder."/".$dbSNPToRetrieve;
        my $Genome1Kindels=$GenomeFolder."/".$Genome1KindelsToRetrieve;		my $checkGenome1Kindels=$GenomeFolder."/".$Genome1KindelsToRetrieve;
        if(!-e $checkDBsnp){                system("wget ftp://gsapubftp-anonymous\@ftp.broadinstitute.org/bundle/b37/$dbSNPToRetrieve -O $dbSNP"); 			system("gunzip $dbSNP");}
        if(!-e $checkGenome1Kindels){       system("wget ftp://gsapubftp-anonymous\@ftp.broadinstitute.org/bundle/b37/$Genome1KindelsToRetrieve -O $Genome1Kindels"); 	system("gunzip $Genome1Kindels");}
}

sub getSamplesFiles{
	my $samplesFile=${$_[0]};

	my %SamplesIDToFiles;
	open(SAMPLES,"<$samplesFile");
	while(<SAMPLES>){
        	chop($_);
        	my @split=split("\t",$_);
        	if($split[0] eq "sample_alias"){ next;}
		my $sampleName=$split[0];
		my $base=basename($split[2]);
		my @splitBase=split("_",$base);
		my $sampleID=$sampleName."_".$splitBase[3];
		push(@{$SamplesIDToFiles{$sampleID}},$split[2]);
	}
	close(SAMPLES);
	return %SamplesIDToFiles;
}

sub unencryptFastq{
        my %SamplesToFile=%{$_[0]}; my $FastqFolder=${$_[1]}; my $currentDir=${$_[2]};

	my $UnecryptionScript=$currentDir."/Unencrypt.sh";
        foreach my $sample(keys(%SamplesToFile)){
                my @splitSample=split("_",$sample);
                my $sampleName=$splitSample[0];	my $sampleLane=$splitSample[1];
		my @fqs=@{$SamplesToFile{$sample}};
		system("$UnecryptionScript $sampleName $sampleLane $fqs[0] $fqs[1] $FastqFolder $currentDir");
        }
}

sub getSamplesUnencryptedFiles{
	my $samplesFile=${$_[0]}; my $FastqFolder=${$_[1]};

        my %SamplesIDToFiles;
        open(SAMPLES,"<$samplesFile");
        while(<SAMPLES>){
                chop($_);
                my @split=split("\t",$_);
                if($split[0] eq "sample_alias"){ next;}
                my $sampleName=$split[0];
                my $base=basename($split[2]);	$base =~ s{\.[^.]+$}{};
		my $fastqFile=$FastqFolder."/".$base;
		my @splitBase=split("_",$base);
		my $sampleID=$sampleName."_".$splitBase[3];
		push(@{$SamplesIDToFiles{$sampleID}},$fastqFile);
	}
	close(SAMPLES);
	return %SamplesIDToFiles;
}
###################################################################################################################
###################################################################################################################
###################################################################################################################

my $currentDir=getcwd;
my $ToolsFolder=$currentDir."/Tools";				if(!-d $ToolsFolder){	system("mkdir $ToolsFolder");}
my $BundleFolder=$currentDir."/BundleFiles";			if(!-d $BundleFolder){	system("mkdir $BundleFolder");}
my $GenomeFolder=$BundleFolder."/GenomeAssembly";		if(!-d $GenomeFolder){	system("mkdir $GenomeFolder");}

####Set output folder to decide where to store produced data#####
my $outputFolder="/nfs/crg/new_tool/spataro/EOSCpipeline";	if(!-d $outputFolder){	system("mkdir $outputFolder");}
################################################################

my $FastqFolder=$outputFolder."/FastqFiles";			if(!-d $FastqFolder){	system("mkdir $FastqFolder");}
my $BAMFolder=$outputFolder."/BAMFiles";			if(!-d $BAMFolder){	system("mkdir $BAMFolder");}
my $VCFFolder=$outputFolder."/VCFFiles";			if(!-d $VCFFolder){	system("mkdir $VCFFolder");}

my $bwaCommand=$ToolsFolder."/bwa-0.5.8c_patched/bwa";
my $ReferenceToRetrieve="human_g1k_v37";
my $dbSNPToRetrieve="dbsnp_138.b37.excluding_sites_after_129.vcf";
my $Genome1KindelsToRetrieve="1000G_phase1.indels.b37.vcf";
####In the following function human reference genome files are downloaded and indexes################################################
####In addition some ancillary files produced by 1000Genomes Project and useful for BAM processing and variant calling are downloaded
&retrieveGenomeData(\$GenomeFolder,\$ReferenceToRetrieve,\$bwaCommand,\$dbSNPToRetrieve,\$Genome1KindelsToRetrieve);
#####################################################################################################################################

my $FromFastqToBamScript=$currentDir."/FromFastqToBam.pl";	####script to process each single pair of fastq files
my $FromBamToVcfScript=$currentDir."/FromBamToVcf.pl";		####script to merge BAM files related to the same individual and for the subsequent variant calling

my $samplesFile=$currentDir."/GonlSamplesToFiles.txt";		####File containing the relationship between sample names and fastq file path
my %SamplesIDToFiles=&getSamplesFiles(\$samplesFile);

#####The following function is required to un-encrypt data stored at EGA#######
#&unencryptFastq(\%SamplesIDToFiles,\$FastqFolder,\$currentDir);
##############################################################################

my %SamplesIDToUncryptedFiles=&getSamplesUnencryptedFiles(\$samplesFile,\$FastqFolder);
my $platform="illumina";
foreach my $sampleID(keys(%SamplesIDToUncryptedFiles)){
	####the sample ID is represented by the sample name and the lane number. The BAM processing will be carried out on fastq files pair related to the same individual and lane
	my @fqFiles=@{$SamplesIDToUncryptedFiles{$sampleID}};
	my @splitID=split("_",$sampleID);	
	my $base1=basename($fqFiles[0]);	my @splitBase1=split("_",$base1);
	my $sampleName=$splitID[0]; my $lane=$splitID[1]; my $library=$splitBase1[4];
	system("perl $FromFastqToBamScript $sampleName $lane $library $platform $fqFiles[0] $fqFiles[1] $BAMFolder $FastqFolder");
}

my %SamplesJoinedLanes;
foreach my $sampleID(keys(%SamplesIDToUncryptedFiles)){
	####in this loop all the recalibrated BAM related to the sample are collected
	my @splitID=split("_",$sampleID);
	my $recalibratedBAM=$BAMFolder."/".$sampleID."/".$sampleID.".recal.sorted.bam";
	push(@{$SamplesJoinedLanes{$splitID[0]}},$recalibratedBAM);
}

foreach my $sample(keys(%SamplesJoinedLanes)){
	#### in this loop the script for BAM merging and variant calling is run on each sample
	if($sample ne "gonl-1a"){ next;}
	my @BAMs=@{$SamplesJoinedLanes{$sample}};
	my $bamString=join(" ",@BAMs);
	system("perl $FromBamToVcfScript $sample $VCFFolder $bamString");
}



