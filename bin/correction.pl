#!/usr/bin/env perl

sub correctSamplesNames{
    my $vcf=${$_[0]};

    my $tmpVCF="tmp.vcf";
    open(TMP,">$tmpVCF");
    open(VCF,"<$vcf");
    while(<VCF>){
            if(substr($_,0,2) eq "##"){ print TMP "$_"; next;}
            if(substr($_,0,1) eq "#"){
                    chop($_);
                    my @split=split("\t",$_);
                    for(my $i=0;$i<9;$i++){ print TMP "$split[$i]\t";}

                    my @samplesFields=splice(@split,9);
                    for(my $i=0;$i<$#samplesFields;$i++){
                            my $sample=$samplesFields[$i];
                            my @splitSample=split(/\./,$sample);
                            print TMP "$splitSample[0]\t";
                    }
                    my $finalSample=$samplesFields[-1];
                    my @splitSample=split(/\./,$finalSample);
                    print TMP "$splitSample[0]\n";
            }
            else{   print TMP "$_";}
    }
    close(VCF); close(TMP);
    system("mv $tmpVCF $vcf");
}

my $PBTvcf=$ARGV[0];
&correctSamplesNames(\$PBTvcf);
