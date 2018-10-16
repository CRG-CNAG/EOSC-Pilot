FROM debian:wheezy
MAINTAINER paolo.ditommaso@gmail.com 

RUN apt-get update --fix-missing \
 && apt-get install -y openjdk-6-jdk \
 && apt-get install -y wget git

#
# compile GATK 1.1 (compilation of version 1.1 works only *after* having compiled version 1.2)
# 
RUN wget -q -O- https://archive.apache.org/dist/ant/binaries/apache-ant-1.8.1-bin.tar.gz | tar xz  \
 && wget -q -O- https://github.com/broadgsa/gatk-protected/archive/1.2.tar.gz | tar xz 
 
ENV PATH=/apache-ant-1.8.1/bin:$PATH

RUN cd /gatk-protected-1.2 && ant 

#
# main image build 
#  
FROM debian:jessie

#
# copies gatk-1.1 distribution 
#
COPY --from=0 /gatk-protected-1.2/dist /gatk-1.2 
 
#
# install Java 7 
# 
RUN apt-get update --fix-missing \
 && apt-get install -y openjdk-7-jre-headless \
 && apt-get install -y unzip wget curl 

#
# Install R
#
RUN apt-get update --fix-missing \
 && apt-get install -y r-base 
 
RUN R -e 'install.packages("gsalib", repos="http://cloud.r-project.org/")'   
 
#
# Install required Perl modules
# 
RUN wget -q cpanmin.us -O /usr/local/bin/cpanm \
  && chmod +x /usr/local/bin/cpanm \
  && cpanm File::Basename Cwd Getopt::Long Switch 
 
#
# Download and install BWA 0.5.8c_patched
#
RUN wget -q https://github.com/molgenis/ngs-utils/raw/master/bwa-0.5.8c_patched.tar.gz \
 && tar xf bwa-0.5.8c_patched.tar.gz \
 && cd bwa-* \
 && mv bwa /usr/local/bin

#
# Install picard 1.32
#
RUN wget -q 'https://downloads.sourceforge.net/project/picard/picard-tools/1.32/picard-tools-1.32.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fpicard%2Ffiles%2Fpicard-tools%2F1.32%2F&ts=1509117698&use_mirror=ayera' -O picard-tools-1.32.zip \
  && unzip picard-tools-1.32.zip
 
#
# Install FastQC 0.7.0
# 
RUN wget -q --no-check-certificate https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.7.0.zip \
 && unzip fastqc_* \
 && chmod +x /FastQC/fastqc \
 && ln -s /FastQC/fastqc /usr/local/bin/fastqc
 
RUN wget -q https://github.com/broadgsa/gatk-protected/archive/1.2.tar.gz \
 && tar xf 1.2.tar.gz 
