# Original pipeline 

The original pipeline is described in the following paper:

Francioli et al, *Whole-genome sequence variation, population structure and demographic history of the Dutch population*, 
[doi:10.1038/ng.3021](http://www.nature.com/ng/journal/v46/n8/full/ng.3021.html) 

### Part 1 

Script describing all the alignment analysis steps we performed for the GoNL project. 
Please note that most software is quite outdated, which in turn means all the commandline 
arguments are not 1 on 1 reusable. We have updated all our protocols, pipelines, etc. 
to current standards, but I don't know if you're interested in these. 

See [part1](part1) folder.


### Part 2 


Attached three of the statististics scripts we used to generate QC reports. Furthermore, the old Picard, 
GATK and verifyBamID we haven't archived, I assume these can be found online at their source page. 
The patched version of BWA we archived here: https://github.com/molgenis/ngs-utils/blob/master/bwa-0.5.8c_patched.tar.gz
This repository might also contain other "old" scripts we have used.

The scripts I previously attached to the email only describe the part of the single flowcell analysis. 
We afterwards merged all flowcells per sample and conducted further downstream analysis. I attached the scripts 
for some coverage QC and genotype calling using GATK UG too.

We currently have much more sophisticated (using Broad/GATK best practices) pipelines in which we automated 
almost everything (from data staging, QCing to GATK haplotypecaller and SV analysis). My colleagues Roan and 
Gerben (in cc.) use these pipelines for diagnostic data analysis and GoNL re-analysis on genome build 38. 
All these pipelines are versioned in github, but I don't know which version etc. they are using at the moment.


See [part2](part2) folder.


### Part 3 

For the GoNL b38 data analysis we've used [this version of the pipeline](https://github.com/molgenis/NGS_DNA/releases/tag/3.4.0). 

As you can see on our github we are now at 3.4.2, but those are minor updates according to the GoNL version.
