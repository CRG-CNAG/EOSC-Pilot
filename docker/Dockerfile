FROM r-base:3.2.3

MAINTAINER paolo.ditommaso@gmail.com

RUN apt-get update --fix-missing \
  && apt-get install -q -y bc wget curl vim nano unzip make gcc g++ gfortran procps \
  && apt-get install -q -y libnet-ssleay-perl libcrypt-ssleay-perl libssl-dev libexpat1 libexpat1-dev liblwp-protocol-https-perl \
  && apt-get install -q -y libglib2.0-0  \
  && apt-get install -q -y openjdk-8-jdk-headless \
  && apt-get clean 
  
RUN wget -q cpanmin.us -O /usr/local/bin/cpanm \
  && chmod +x /usr/local/bin/cpanm   
  
RUN cpanm File::Basename Cwd Getopt::Long Switch

RUN  R -e 'install.packages("gsalib")' 
  
RUN wget -q https://github.com/molgenis/ngs-utils/raw/master/bwa-0.5.8c_patched.tar.gz \
 && tar xf bwa-0.5.8c_patched.tar.gz

RUN wget 'https://downloads.sourceforge.net/project/picard/picard-tools/1.32/picard-tools-1.32.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fpicard%2Ffiles%2Fpicard-tools%2F1.32%2F&ts=1509117698&use_mirror=ayera' -O picard-tools-1.32.zip \
  && unzip picard-tools-1.32.zip
  
 

