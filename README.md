# EOSC pilot project

To run the pipeline run `perl LaunchPipeline.pl`. This script will:

* download the human genome reference file and index them
* download some ancillary files produced by 1000Genomes Project from the Bundle repository of GATK
* run-encrypt EGA deposited fastq files
* run the script "FromFastqToBam.pl" (this script is run on each pair of fastq files and it is able to perform the various steps to obtain the recalibrated BAM)
* run the script "FromBamToVcf.pl" (this script is run on each sample and firstly the various BAM recalibrated files related to the same individual are merged). 
* Subsequently variant calling is performed
